
from __future__ import annotations

import numpy as np

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
        backend: BackendV1,
        *,
        interacting_sites: int = 1,
        ansatz: QuantumCircuit | None=None, 
        optimizer: Optimizer | Minimizer | None=None  
    ) -> None:
        """
        Description:
        Args:
            ...
        """
        #Implement initialization
    
        
    # Implement getters and setters

    # Methods to implement
    def hopping():
        ...
    
    def create_hamiltonian():
        ...

    def hamiltonian_to_qubit():
        ...
    
    def compute_band_structure():
        ...

    def plot_band_structure():
        ...