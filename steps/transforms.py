#DMB changes 13/06/2022
################################################################################
# © Copyright 2020-2022 Zapata Computing Inc.
################################################################################
import time
from typing import Union

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
        input_operator = get_fermion_operator(input_operator)
        transformation_function = jordan_wigner
    else:
        raise RuntimeError("Unrecognized transformation ", transformation)

    start_time = time.time()
    if transformation == "BK-2qbr":
            transformed_operator = transformation_function(input_operator,active_orbitals,active_fermions)
    else:
        transformed_operator = transformation_function(input_operator)
    walltime = time.time() - start_time
    
   if transformation == "qiskit":
    #transform orquestra qubit operator (openfermion) to qiskit
        qiskitop=qubitop_to_qiskitpauli(input_operator)
    #perform particle/hole tranformation
    #transform qiskit back to OF representation
        transformed_operator=qiskitpauli_to_qubitop(qiskitop)

    save_qubit_operator(transformed_operator, "transformed-operator.json")
    save_timing(walltime, "timing.json")
