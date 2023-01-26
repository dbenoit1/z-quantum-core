#DMB changes 13/06/2022
################################################################################
# © Copyright 2020-2022 Zapata Computing Inc.
################################################################################
import time
from typing import Union
import numpy as np

from zquantum.core.openfermion import (
    SymbolicOperator,
    bravyi_kitaev,
    get_fermion_operator,
    jordan_wigner,
    symmetry_conserving_bravyi_kitaev,
    load_interaction_operator,
    save_qubit_operator,
)
from zquantum.core.utils import save_timing
from qeqiskit.conversions import qiskitpauli_to_qubitop, qubitop_to_qiskitpauli
from qiskit.chemistry import FermionicOperator


def transform_interaction_operator(
    transformation: str, input_operator: Union[str, SymbolicOperator], active_orbitals=0,active_fermions=0
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
    else:
        raise RuntimeError("Unrecognized transformation ", transformation)

    start_time = time.time()
    if transformation == "BK-2qbr":
            transformed_operator = transformation_function(input_operator,active_orbitals,active_fermions)
    elif transformation == "qiskit":
#    if transformation == "qiskit":
#    #transform orquestra qubit operator (openfermion) to qiskit
#        qiskitop=qubitop_to_qiskitpauli(transformed_operator)
#        print("OK")
#    #transform qiskit back to OF representation
#        transformed_operator=qiskitpauli_to_qubitop(qiskitop)
          
#    if transformation == "qiskit":
#    #transform orquestra qubit operator (openfermion) to qiskit
#        # This shoudl work, but doesn't convert -> qiskitop=qubitop_to_qiskitpauli(input_operator)
#        qiskitop=qubitop_to_qiskitpauli(transformed_operator)
#    #perform particle/hole tranformation (assumes number active_fermions = 2*alpha =2*beta) 
#    #QISKIT CODE \/
#        newferOp, energy_shift = qiskitop.to_legacy_op().particle_hole_transformation([active_fermions//2, active_fermions//2])
#        print('Energy shift is: {}'.format(energy_shift))
#      #  newqubitOp_jw = newferOp.mapping(map_type='JORDAN_WIGNER', threshold=0.00000001)
#      #  newqubitOp_jw.chop(10**-10)
#    #QISKIT CODE /\
#    #transform qiskit back to OF representation
#        transformed_operator=qiskitpauli_to_qubitop(newferOp)
#
        #instead extract h1 and h2 from the interaction operator
        h1=input_operator.one_body_tensor       
        h2=input_operator.two_body_tensor
        print(h1)
        print("++++++++")
        print(h2)
        #print("++++++++")
        #print(np.einsum(‘ikmj->ijkm’, h2))
        #DEBUG qiskitop=qubitop_to_qiskitpauli(jordan_wigner(get_fermion_operator(input_operator)))
        #DEBUG print(qiskitop)
        #then use qiskit to make the operator
        #QISKIT CODE \/
        ferOp = FermionicOperator(h1=h1, h2=h2)
        #perform particle/hole tranformation (assumes number active_fermions = 2*alpha =2*beta) 
        newferOp, energy_shift = ferOp.particle_hole_transformation([active_fermions//2, active_fermions//2])
        print('Energy shift is: {}'.format(energy_shift))
        #map to JW in qiskit
        newqubitOp_jw = newferOp.mapping(map_type='JORDAN_WIGNER', threshold=0.00000001)
        newqubitOp_jw.chop(10**-10)
        #QISKIT CODE /\
        #transform qiskit back to OF representation
        transformed_operator=qiskitpauli_to_qubitop(newqubitOp_jw.to_opflow().to_pauli_op())
    else:
        transformed_operator = transformation_function(input_operator)

    walltime = time.time() - start_time

    save_qubit_operator(transformed_operator, "transformed-operator.json")
    save_timing(walltime, "timing.json")
