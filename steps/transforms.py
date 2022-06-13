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


def transform_interaction_operator(
    transformation: str, input_operator: Union[str, SymbolicOperator]
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
    else:
        raise RuntimeError("Unrecognized transformation ", transformation)

    start_time = time.time()
    transformed_operator = transformation_function(input_operator)
    walltime = time.time() - start_time

    save_qubit_operator(transformed_operator, "transformed-operator.json")
    save_timing(walltime, "timing.json")
