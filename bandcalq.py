
from __future__ import annotations

import numpy as np
import cmath as cmt

from qiskit.opflow import Z, I, X, Y, PauliOp
from qiskit.algorithms.eigensolvers import VQD
from qiskit.algorithms.minimum_eigensolvers import VQE
from qiskit.algorithms.optimizers import Optimizer, Minimizer
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import EfficientSU2
from qiskit.providers import BackendV1


class BandCalQ():
    """
    Description: 
    ...
    """
    def __init__(
        self,
        orbital_number: int,
        displacement: float,
        hopping_matrix: np.ndarray,
        backend: BackendV1,
        *,
        interacting_sites: int = 1,
        ansatz: QuantumCircuit | None=None, 
        optimizer: Optimizer | Minimizer | None=None, 
    ) -> None:
        """
        Description:
        Args:
            ...
        """
        
        self.orbital_number = orbital_number
        self.displacement = displacement
        self.hopping_matrix = hopping_matrix
        self.backend = backend
        self.interacting_sites = interacting_sites
        self.ansatz = ansatz
        self.optimizer = optimizer
        
    # Implement getters and setters
    #
    
    def hopping(
        self,
        alpha: int,
        beta: int,
        delta: float,
    ) -> float:
        return self.hopping_matrix[alpha][beta]
    
    def create_hamiltonian(self, k) -> np.ndarray:
        hamiltonian = np.zeros((self.orbital_number, self.orbital_number), dtype=complex)
        for alpha in range(self.orbital_number):
            for beta in range(self.orbital_number):
                for delta in range(1, self.interacting_sites + 1):
                    hamiltonian[alpha][beta] = hamiltonian[alpha][beta] + \
                    self.hopping(alpha, beta, -delta)*cmt.exp(-1j*k*delta*self.displacement) + \
                    self.hopping(alpha, beta, delta)*cmt.exp(1j*k*delta*self.displacement) 
        return hamiltonian

    @classmethod
    def operator_extended_one(
        cls,
        operator: PauliOp, 
        position: int, 
        size: int,
    ) -> PauliOp:
        """Returns PauliOp consisting of I's and operator on provided position of given size"""
        ext_op = I
        for i in range(size):
            if i == 0:
                if i == position:
                    ext_op = operator
            elif i == position:
                ext_op ^= operator
            else:
                ext_op ^= I
        
        return ext_op

    @classmethod
    def operator_extended_two(
        cls,
        operator_1: PauliOp,
        operator_2: PauliOp,
        position_1: int,
        position_2: int,
        size: int,
    ) -> PauliOp:
        """Returns PauliOp consisting of I's, operator_1 on position_1 and operator_2 on position_2 of given size"""
        ext_op = I
        for i in range(size):
            if i == 0:
                if i == position_1:
                    ext_op = operator_1
                elif i == position_2:
                    ext_op = operator_2
            elif i == position_1:
                ext_op ^= operator_1
            elif i == position_2:
                ext_op ^= operator_2
            else:
                ext_op ^= I
        
        return ext_op
    
    # Methods to implement
    def hamiltonian_to_qubit() -> None:
        ...
    
    def compute_band_structure():
        ...

    def plot_band_structure():
        ...