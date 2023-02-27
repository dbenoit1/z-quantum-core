#DMB changes 13/06/2022
################################################################################
# © Copyright 2020-2022 Zapata Computing Inc.
################################################################################
import time
from typing import Union
import numpy as np
import json

from zquantum.core.openfermion import (
    SymbolicOperator,
    bravyi_kitaev,
    get_fermion_operator,
    jordan_wigner,
    symmetry_conserving_bravyi_kitaev,
    load_interaction_operator,
    save_qubit_operator,
    QubitOperator,
)
from zquantum.core.utils import save_timing
from qeqiskit.conversions import qiskitpauli_to_qubitop, qubitop_to_qiskitpauli
from qiskit.chemistry import FermionicOperator
from zquantum.core.openfermion.ops.representations import InteractionOperator

def rearrange_both(array,array2,occ):
    #rearrange a block array of the type AAAABBBB to ABABABAB, assuming that we have occ A elements
    #Performs both one-body and two-body intergral matrices rearrengments - checked on H2/STO-3G 
    #with original qiskit code
    idx=[]
    base=0
    for i in range(len(array)//2):
        idx.append(base)
        idx.append(base+occ)
        base+=1
    print("new state ordering index array",idx)
    new_array=array[:, idx][idx]
    new_array2=array2[:, :, :, idx][:, :, idx][:, idx][idx]
    return(new_array,new_array2)

def transform_interaction_operator(
    transformation: str, input_operator: Union[str, SymbolicOperator], active_orbitals=0,active_fermions=0,hfresults=""
):
    """Transform an interaction operator through the Bravyi-Kitaev, 2qbitrer-BK or
    Jordan-Wigner transformations. The results are serialized into a JSON under the
    files: "transformed-operator.json" and "timing.json"

    ARGS:
        transformation: The transformation to use.  "Jordan-Wigner" or
            "Bravyi-Kitaev" or "BK-2qbr"
        input_operator: The interaction operator to transform
    """
    if isinstance(input_operator, str):
        input_operator = load_interaction_operator(input_operator)

    if transformation == "Jordan-Wigner":
        transformation_function = jordan_wigner
    elif transformation == "Bravyi-Kitaev":
        input_operator = get_fermion_operator(input_operator)
        transformation_function = bravyi_kitaev
    elif transformation == "BK-2qbr":
        input_operator = get_fermion_operator(input_operator)
        transformation_function = symmetry_conserving_bravyi_kitaev
    elif transformation == "qiskit":
     #   input_operator = get_fermion_operator(input_operator)
        transformation_function = jordan_wigner
    elif transformation == "qiskitBK":
        print("COMPUTING PH Hamiltonian")
        # extract h1 and h2 from the interaction operator to generate one in qiskit
        h1=input_operator.one_body_tensor 
        # there is a sign difference between OF and qiskit so changing sign.
        h2=-input_operator.two_body_tensor
        #note that the order is diffrent between the interleaved format that qiskit expects and the block format that is provided
        #re-organising the data to fit an interleaved scheme rather than a block scheme
        (newh1,newh2)=rearrange_both(h1,h2,active_fermions)
        #QISKIT CODE \/
        #then use qiskit to make a fermionic operator operator
        ferOp = FermionicOperator(h1=newh1, h2=newh2)
        #perform particle/hole tranformation (assumes number active_fermions = 2*alpha =2*beta) 
        newferOp, energy_shift = ferOp.particle_hole_transformation([active_fermions//2, active_fermions//2])
        print('Energy shift is: {}'.format(energy_shift))
        #QISKIT CODE /\
        #extract H1 and H2 from PH hamiltonian
        (phh1,phh2)=rearrange_both(newferOp.h1,newferOp.h2,active_fermions)
        #Flip sign for H2
        phh2*=-1
        input_operator = InteractionOperator(energy_shift, phh1, phh2)
        print(input_operator)
        input_operator = get_fermion_operator(input_operator)
        print(input_operator)
        #setting transformation to BK-2qbr
        #transformation_function = symmetry_conserving_bravyi_kitaev
        transformation_function = jordan_wigner #test
    else:
        raise RuntimeError("Unrecognized transformation ", transformation)

    start_time = time.time()
    #if transformation == "BK-2qbr" or "qiskitBK":
    if transformation == "BK-2qbr":
            transformed_operator = transformation_function(input_operator,active_orbitals,active_fermions)
    elif transformation == "qiskit":
        # extract h1 and h2 from the interaction operator to generate one in qiskit
        h1=input_operator.one_body_tensor 
        # there is a sign difference between OF and qiskit so changing sign.
        h2=-input_operator.two_body_tensor
        #note that the order is diffrent between the interleaved format that qiskit expects and the block format that is provided
        #re-organising the data to fit an interleaved scheme rather than a block scheme
        (newh1,newh2)=rearrange_both(h1,h2,active_fermions)
        #print(newh1)
        #print(".......")
        #print(newh2)

        print(input_operator.constant)
        #QISKIT CODE \/
        #then use qiskit to make a fermionic operator operator
        ferOp = FermionicOperator(h1=newh1, h2=newh2)
        #perform particle/hole tranformation (assumes number active_fermions = 2*alpha =2*beta) 
        newferOp, energy_shift = ferOp.particle_hole_transformation([active_fermions//2, active_fermions//2])
        print('Energy shift is: {}'.format(energy_shift))
        #map to JW in qiskit
        newqubitOp_jw = newferOp.mapping(map_type='JORDAN_WIGNER', threshold=0.00000001)
        newqubitOp_jw.chop(10**-10)
        #operator is not in the correct format 
        print(newqubitOp_jw)
        #change this to a SummedOp (even though the qeqiskit conversion routine doesnt recognise it as such 
        #(this is current fixed by removing the safety check in "qiskitpauli_to_qubitop" but fix in a cleaner way later 
        newqubitOp_jw_P=newqubitOp_jw.to_opflow()
        #checking that this is correct
        print("Qiskit version:")
        print(newqubitOp_jw_P)
        #QISKIT CODE /\
        #transform qiskit back to OF representation
        transformed_operator=qiskitpauli_to_qubitop(newqubitOp_jw_P)
        #checking the resulting operator
        print("OF version:")
        print(transformed_operator)
        #IMPORTANT NOTE: the OF tranformation removes the IIII operations so the resulting H needs shifing by the IIII coefficent
        #coeff is stored in (first element of the operator, usually ordered in qiskit, with II..I being first):
        print(newqubitOp_jw_P._oplist[0].coeff)
        ##Check if a HF ground state energy is provided and if so use it to "zero" the hamiltonian -> correlation hamiltonian
        #hf_ground_state_energy=0
        ##try to read the HF energy from psi4 calculation
        #if (hfresults != ""):
        #    #If the file is there, then read in the HF energy
        #    with open(hfresults, "r") as f:
        #        psi4results = json.load(f)
        #    hf_ground_state_energy=float(psi4results["energy"])
        #    print("HF energy:",hf_ground_state_energy)
        ##remove HF energy from constant shift (if this was not read in, then simply remove 0)  
        #constant_coefficient=newqubitOp_jw_P._oplist[0].coeff-energy_shift-hf_ground_state_energy
        constant_coefficient=newqubitOp_jw_P._oplist[0].coeff
        #Print the correction value
        print("total correction = ",constant_coefficient)
        #Add constant term to the operator
        transformed_operator+= constant_coefficient * QubitOperator(())
    else:
        transformed_operator = transformation_function(input_operator)
        print("FINAL T-OP")
        print(transformed_operator)
        
    walltime = time.time() - start_time

    save_qubit_operator(transformed_operator, "transformed-operator.json")
    save_timing(walltime, "timing.json")
