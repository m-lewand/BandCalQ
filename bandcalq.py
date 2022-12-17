from __future__ import annotations

import numpy as np
import cmath as cmt
import matplotlib.pyplot as plt
from qiskit.opflow import Z, I, X, Y, PauliOp
from vqd_custom import VQD
from qiskit.algorithms import NumPyEigensolver
from qiskit.algorithms.optimizers import Optimizer, Minimizer, SPSA
from qiskit.algorithms.state_fidelities import ComputeUncompute
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import EfficientSU2
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
        lattice_constant: float,
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
        self.lattice_constant = lattice_constant
        self.hopping_matrix = hopping_matrix
        self.backend = backend
        self.interacting_sites = interacting_sites
        self.ansatz = ansatz
        self.optimizer = optimizer
        
    @property
    def ansatz(self):
        return self._ansatz
    
    @property
    def optimizer(self):
        return self._optimizer

    @ansatz.setter
    def ansatz(self, arg):
        if arg == None:
            self.ansatz = EfficientSU2(self.orbital_number, su2_gates=['rx', 'rz', 'ry'], entanglement='full', reps=2)
        else:
            self._ansatz = arg

    @optimizer.setter
    def optimizer(self, arg):
        if arg == None:
            self._optimizer = SPSA(maxiter=400)
        else:
            self._optimizer = arg
            
    
    def hopping(
        self,
        alpha: int,
        beta: int,
        delta: float,
    ) -> float:
        return self.hopping_matrix[alpha][beta]
    
    def create_hamiltonian(self, momentum: float) -> np.ndarray:
        hamiltonian = np.zeros((self.orbital_number, self.orbital_number), dtype=complex)
        for alpha in range(self.orbital_number):
            for beta in range(self.orbital_number):
                for delta in range(1, self.interacting_sites + 1):
                    hamiltonian[alpha][beta] = hamiltonian[alpha][beta] + \
                    self.hopping(alpha, beta, -delta)*cmt.exp(-1j*momentum*delta*self.lattice_constant) + \
                    self.hopping(alpha, beta, delta)*cmt.exp(1j*momentum*delta*self.lattice_constant) 
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

    def compute_band_structure(
        self, 
        momentum_min: float,
        momentum_max: float,
        momentum_points_amount: int,
        *,
        theoretical_points: bool = False,
    ) -> None:
        '''Description'''
        self.momentum_points_amount = momentum_points_amount
        self.eigenvalues_array = np.zeros((2*self.orbital_number, self.momentum_points_amount))
        self.momentum_array = np.linspace(momentum_min/(np.pi/self.lattice_constant), momentum_max/(np.pi/self.lattice_constant), momentum_points_amount)
        if theoretical_points:
            self.eigenvalues_array_theoretical = np.zeros((2*self.orbital_number, self.momentum_points_amount))
            solver = NumPyEigensolver(k=2*self.orbital_number)
        
        for i in range(self.momentum_points_amount):
            self.create_hamiltonian_qubit(self.momentum_array[i])
            if theoretical_points:
                self.eigenvalues_array_theoretical[:,i] = solver.compute_eigenvalues(self.hamiltonian_qubit).eigenvalues

            vqd_algorithm = VQD(ansatz=self.ansatz, estimator=BackendEstimator(self.backend), optimizer=self.optimizer,
                                fidelity=ComputeUncompute(sampler=Sampler()), k=2*self.orbital_number)
            vqd_result = vqd_algorithm.compute_eigenvalues(self.hamiltonian_qubit)
            self.eigenvalues_array[:,i] =   np.real(vqd_result.eigenvalues)
            
        if theoretical_points:
            self.theory_computed = True
        else:
            self.theory_computed = False
        return
 
    def plot_band_structure(
        self,
        *,
        theoretical_points: bool=False,
        save_png: bool=False,
        png_name: str=""
        ):
        
        if theoretical_points and not(self.theory_computed):
            self.eigenvalues_array_theoretical = np.zeros((2*self.orbital_number, self.momentum_points_amount))
            solver = NumPyEigensolver(k=2*self.orbital_number)
            for i in range(self.momentum_points_amount):
                self.create_hamiltonian_qubit(self.momentum_array[i])
                self.eigenvalues_array_theoretical[:,i] = solver.compute_eigenvalues(self.hamiltonian_qubit).eigenvalues
            self.theory_computed = True
        
        plt.figure(1, figsize=(5,5))
        plt.rcParams.update({'font.size' : 12})
        ax = plt.gca()

        for i in range(2*self.orbital_number):
            if theoretical_points:
                color = next(ax._get_lines.prop_cycler)['color']
                plt.plot(self.momentum_array, self.eigenvalues_array[i], 
                marker='o', markersize=5, color=color,mfc='white', linestyle='None')
                plt.plot(self.momentum_array, self.eigenvalues_array_theoretical[i], 
                marker='v', markersize=4, color=color, mfc='black', linestyle='-', alpha=0.9)
            else:
                plt.plot(self.momentum_array, self.eigenvalues_array[i], 
                marker='o', markersize=4, mfc='white',linestyle='--')
        plt.grid()
        plt.xlabel('$k (\pi/a)$')
        plt.ylabel('$E(E_h)$')
        if save_png:
            fig = plt.gcf()
            if png_name == '':
               fig.savefig('band.png', format='png', dpi=1200, facecolor='w')
            else:
               fig.savefig(str(png_name)+'.png', format='png', dpi=1200, facecolor='w')  

        plt.show()
        return
