from typing import Any, Optional

import numpy as np
from zquantum.core.circuits import Circuit
from zquantum.core.circuits.layouts import CircuitConnectivity
from zquantum.core.interfaces.backend import QuantumSimulator, flip_wavefunction
from zquantum.core.measurement import Measurements, sample_from_wavefunction
from zquantum.core.wavefunction import Wavefunction


class SymbolicSimulator(QuantumSimulator):
    """A simulator computing wavefunction by consecutive gate matrix multiplication.

    Args:
        seed: the seed of the sampler
    """

    def __init__(
        self,
        noise_model: Optional[Any] = None,
        device_connectivity: Optional[CircuitConnectivity] = None,
        seed: Optional[int] = None,
    ):
        super().__init__(noise_model, device_connectivity)
        self._seed = seed

    def run_circuit_and_measure(self, circuit: Circuit, n_samples: int) -> Measurements:
        """Run a circuit and measure a certain number of bitstrings

        Args:
            circuit: the circuit to prepare the state
            n_samples: the number of bitstrings to sample
        """
        if circuit.free_symbols:
            raise ValueError("Cannot sample from circuit with symbolic parameters.")

        wavefunction = self.get_wavefunction(circuit)
        bitstrings = sample_from_wavefunction(wavefunction, n_samples, self._seed)
        return Measurements(bitstrings)

    def get_wavefunction(self, circuit: Circuit) -> Wavefunction:
        if circuit.free_symbols:
            raise ValueError("Currently circuits with free symbols are not supported")
        super().get_wavefunction(circuit)
        state = np.zeros(2 ** circuit.n_qubits)
        state[0] = 1

        for operation in circuit.operations:
            state = operation.apply(state)

        return flip_wavefunction(Wavefunction(state))
