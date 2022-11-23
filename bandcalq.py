
from __future__ import annotations

import numpy as np
import cmath as cmt

from qiskit.opflow import Z, I, X, Y, PauliOp, PauliSumOp
#from qiskit.algorithms.eigensolvers import VQD
from vqd_fixed import VQD
#
from qiskit.algorithms.eigensolvers import NumPyEigensolver
#
from qiskit.algorithms.minimum_eigensolvers import VQE
from qiskit.algorithms.optimizers import Optimizer, Minimizer, SPSA
from qiskit.algorithms.state_fidelities import ComputeUncompute
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import EfficientSU2, TwoLocal
from qiskit.providers import BackendV1
from qiskit.primitives import BackendEstimator, Sampler



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
    
    def create_hamiltonian(self, momentum) -> np.ndarray:
        hamiltonian = np.zeros((self.orbital_number, self.orbital_number), dtype=complex)
        for alpha in range(self.orbital_number):
            for beta in range(self.orbital_number):
                for delta in range(1, self.interacting_sites + 1):
                    hamiltonian[alpha][beta] = hamiltonian[alpha][beta] + \
                    self.hopping(alpha, beta, -delta)*cmt.exp(-1j*momentum*delta*self.displacement) + \
                    self.hopping(alpha, beta, delta)*cmt.exp(1j*momentum*delta*self.displacement) 
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
    
    def create_hamiltonian_qubit(self, momentum: float) -> None:
        hamiltonian = self.create_hamiltonian(momentum)
        hamiltonian_qubit = 0
        I_op = I
        for i in range(self.orbital_number-1):
            I_op ^= I
        
        for alpha in range(self.orbital_number):
            hamiltonian_qubit += 0.5*hamiltonian[alpha][alpha]*I_op
            hamiltonian_qubit -= 0.5*hamiltonian[alpha][alpha]*self.operator_extended_one(Z, alpha, self.orbital_number)
        
            beta = alpha + 1

            while(beta < self.orbital_number):
                hamiltonian_qubit += 0.5*(hamiltonian[alpha][beta]).real*(self.operator_extended_two(X, X, alpha, beta, self.orbital_number)) + \
                0.5*(hamiltonian[alpha][beta]).real*(self.operator_extended_two(Y, Y, alpha, beta, self.orbital_number))# + \
                0.5*(hamiltonian[alpha][beta]).imag*(self.operator_extended_two(Y, X, alpha, beta, self.orbital_number)) - \
                0.5*(hamiltonian[alpha][beta]).imag*(self.operator_extended_two(X, Y, alpha, beta, self.orbital_number)) 
                beta += 1
        
        self.hamiltonian_qubit = hamiltonian_qubit

        return
    
    def get_betas(self):
        quick_ansatz = TwoLocal(self.orbital_number , rotation_blocks='ry', entanglement_blocks='cz', reps=2)
        quick_optimizer = SPSA(maxiter=50)
        
        initial_point = np.random.random(quick_ansatz.num_parameters)
        
        vqe_algorithm = VQE(ansatz=quick_ansatz, estimator=BackendEstimator(self.backend), optimizer=quick_optimizer, initial_point=initial_point)

        energy_min = vqe_algorithm.compute_minimum_eigenvalue(self.hamiltonian_qubit)
        energy_max = vqe_algorithm.compute_minimum_eigenvalue(-self.hamiltonian_qubit)

        betas = np.zeros([2*self.orbital_number])

        for i in range(2*self.orbital_number):
            betas[i] = 2*(-np.real(energy_max.eigenvalue) - np.real(energy_min.eigenvalue))
        
        return betas

    def compute_band_structure(
        self, 
        momentum_range: float,
        momentum_points_amount: int,
    ) -> None:
        '''Description'''
        dk = 2*momentum_range/(momentum_points_amount - 1)
        self.momentum_array = np.zeros(momentum_points_amount)
        for i in range(momentum_points_amount):
            self.momentum_array[i] = -momentum_range + dk*i
        self.eigenvalues_array = np.zeros((momentum_points_amount, 2*self.orbital_number))

        #default parameters - modify to use setters
        self.ansatz = EfficientSU2(self.orbital_number, su2_gates=['rx', 'rz', 'ry'], entanglement='full', reps=2)
        self.optimizer = SPSA(maxiter=400)

        for i in range(momentum_points_amount):
            self.create_hamiltonian_qubit(self.momentum_array[i])
            beta_parameters = self.get_betas()

            vqd_algorithm = VQD(ansatz=self.ansatz, estimator=BackendEstimator(self.backend), optimizer=self.optimizer,
                                fidelity=ComputeUncompute(sampler=Sampler()), k=2*self.orbital_number, betas=beta_parameters)
            vqd_result = vqd_algorithm.compute_eigenvalues(self.hamiltonian_qubit)
            self.eigenvalues_array[i,:] =   np.real(vqd_result.eigenvalues)

        return
    # Methods to implement 
    def plot_band_structure():
        ...