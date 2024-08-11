"""
Contains classes and functions related to the optimization of laminate stacking sequences based on generic stacks (see README.md for more information).
"""

import math
import time
import warnings
from typing import Optional, Any, Callable, Final

import numpy as np
import scipy.optimize as opt

import auxs.configparser_auxs as cp_aux
import auxs.scipy_auxs as sp_aux
import clt
import file_handling as fh


class GenericPanel:
    """
    Naming convention:
    x is the list of parameters to be optimized. It contains the continuous thickness values of the generic plies.
    n is the number of ply groups.
    """
    constr: clt.LaminateConstraints
    mat: clt.TransverselyIsotropicMaterial
    dims: clt.Dimensions
    x: Optional[np.ndarray] = None

    def __init__(self, constr: clt.LaminateConstraints, mat: clt.TransverselyIsotropicMaterial, dims: clt.Dimensions) -> None:
        self.mat = mat
        self.constr = constr
        self.dims = dims

    def weight(self, x: np.ndarray) -> float:
        return 2 * sum(x) * self.dims.a * self.dims.b * self.mat.rho

    def create_RF_biaxial_shear(self, n: int, lc: clt.MembraneLoadCase) -> Callable[[np.ndarray], float]:
        def RF_biaxial_shear(x: np.ndarray) -> float:
            return clt.Panel(clt.Laminate([clt.Ply(self.mat, t) for t in x], angles=[theta for theta in n * self.constr.ply_group], symmetric=True), self.dims).RF_biaxial_shear(lc)

        return RF_biaxial_shear

    def create_scipy_bounds(self, n: int) -> opt.Bounds:
        lower_bounds: list[float] = ([min(max(self.constr.min_num_of_covering_plies, 1), self.constr.max_num_of_contiguous_plies) * self.constr.ply_thickness] +
                                     [self.constr.ply_thickness] * (n * len(self.constr.ply_group) - 1))
        upper_bounds: list[float] = ([self.constr.max_num_of_contiguous_plies * self.constr.ply_thickness] * (n * len(self.constr.ply_group) - 1) +
                                     [self.constr.max_num_of_contiguous_plies * self.constr.ply_thickness / 2])
        # Note: 'keep_feasible' is currently not set to True because it results in a ValueError ('x0 violates bound constraints').
        # I don't know why this is the case as every element in x0 is the midpoint between its respective lower and upper bound, so the bounds cannot be not violated.
        return opt.Bounds(lb=lower_bounds, ub=upper_bounds)

    def create_scipy_constraints(self, n: int, lc: clt.MembraneLoadCase) -> list[opt.LinearConstraint | opt.NonlinearConstraint]:
        # The combined reserve factor must be greater or equal to 1:
        constr: list[opt.LinearConstraint | opt.NonlinearConstraint] = [opt.NonlinearConstraint(fun=self.create_RF_biaxial_shear(n, lc), lb=1, ub=np.inf)]

        # Every ply angle must account for at least 'minimum ply share' of the stack's thickness (usually the '10% rule'):
        if self.constr.min_ply_share > 0:
            angles: set[float] = set(self.constr.ply_group)

            # Create A1 such that
            # 0 mm = lower bounds lb <= A1.dot(x1) = total thickness - minimum thickness = excess thickness <= upper bounds ub = infinity
            # for a single ply group x1:
            A1: np.ndarray = - np.full(shape=(len(angles), len(self.constr.ply_group)), fill_value=self.constr.min_ply_share)
            for i, theta_i in enumerate(angles):
                for j, theta_j in enumerate(self.constr.ply_group):
                    if theta_i == theta_j:
                        A1[i, j] += 1

            # Create An for n of ply groups:
            An: np.ndarray = np.block([A1] * n)

            # Append the constraint to the list of constraints:
            constr.append(opt.LinearConstraint(A=An, lb=0, ub=np.inf, keep_feasible=True))

        # The laminate must be balanced:
        if self.constr.balanced:
            angles: set[float] = set(np.abs(self.constr.ply_group)) - {0, 90, 180, 360}

            # Create A1 such that
            # -0.001 mm <= A1.dot(x1) = difference in total thickness between +- plies <= 0.001 mm
            # for a single ply group x1:
            A1: np.ndarray = np.zeros(shape=(len(angles), len(self.constr.ply_group)))
            for i, theta_i in enumerate(angles):
                for j, theta_j in enumerate(self.constr.ply_group):
                    if theta_i == abs(theta_j):
                        A1[i, j] += np.sign(theta_j)

            # Create An for n of ply groups:
            An: np.ndarray = np.block([A1] * n)

            # Append the constraint to the list of constraints:
            constr.append(opt.LinearConstraint(A=An, lb=-self.constr.balance_tolerance, ub=self.constr.balance_tolerance, keep_feasible=True))

        # Return the list of constraints:
        return constr

    def estimate_n(self, lc: clt.MembraneLoadCase, m_max: int = 5, n_max: int = 3) -> int:
        # Calculate the homogenized engineering constants of an approximately equivalent quasi-isotropic laminate:
        iso_lam = clt.Laminate(plies=len(self.constr.ply_group) * [clt.Ply(mat=self.mat, t=1)], angles=[theta for theta in self.constr.ply_group], symmetric=True)
        # noinspection PyTypeChecker
        iso_ABD_inv: np.ndarray[float] = np.linalg.inv(iso_lam.ABD)
        E_iso = 1 / max(iso_ABD_inv[0, 0], iso_ABD_inv[1, 1]) / len(self.constr.ply_group)
        nu_iso = - iso_ABD_inv[0, 1] / max(iso_ABD_inv[0, 0], iso_ABD_inv[1, 1])

        # Calculate the buckling factor for biaxial loading using the formula for isotropic materials:
        k_biaxial = np.inf
        alpha = self.dims.a / self.dims.b
        beta = lc.Ny / lc.Nx
        for m in range(1, m_max + 1):
            for n in range(1, n_max + 1):
                if (m ** 2 + beta * n ** 2 * alpha ** 2) == 0:
                    continue
                k_biaxial_new = (m ** 2 + (n * alpha) ** 2) ** 2 / alpha ** 2 / (m ** 2 + beta * (n * alpha) ** 2)
                if 0 < k_biaxial_new < k_biaxial:
                    if m == m_max:
                        print('Warning: m_max in gs_opt.GenericPanel.estimate_n has been reached.\n')
                    if n == n_max:
                        print('Warning: n_max in gs_opt.GenericPanel.estimate_n has been reached.\n')
                    k_biaxial = k_biaxial_new

        # Calculate the buckling factor for shear loading using the formula for isotropic materials:
        k_shear = 4 + 5.34 / alpha ** 2 if alpha < 1 else 5.34 + 4 / alpha ** 2

        # Calculate the thickness of the isotropic panel resulting in a combined reserve factor of 1 using the formula for isotropic materials:
        t_est = (12 * self.dims.b ** 2 * (1 + nu_iso ** 2) * abs(lc.Nx) / np.pi ** 2 / E_iso * (1 / 2 / k_biaxial + np.sqrt(1 / 4 / k_biaxial + 1 / k_shear ** 2))) ** (1 / 3)

        # 0.5 is a correction factor to account for the fact that probably not all generic plies will have the maximum thickness.
        # This value worked well in most of the cases I tried. Still, the estimation of the number of ply groups is subject to future improvements.
        # TODO: Improve the estimation of the number of ply groups.
        n_est = round(t_est / (2 * self.constr.ply_thickness * 0.5 * self.constr.max_num_of_contiguous_plies * len(self.constr.ply_group)))
        print(f'Estimated number of generic ply groups: {n_est}\n')
        return n_est

    def optimize(self, lc: clt.MembraneLoadCase, options: Optional[dict[str, Any]] = None) -> opt.OptimizeResult:
        """
        Optimizes the continuous thickness values `x` of the generic panel using SciPy's trust-constr method.
        :param lc: the load case to be optimized the panel for as an instance of clt.LoadCase
        :param options: solver options, most notably 'gtol', 'xtol', 'barrier_tol' and 'verbose'
        :return: scipy.optimize.OptimizeResult instance containing the result of the optimized thickness values and metadata
        """
        # Set the default solver options:
        if options is None:
            '''
            0 (default): work silently.
            1: display a termination report.
            2: display progress during iterations.
            3: display progress during iterations (more detailed).
            
            Source: SciPy 1.14.0 API reference
            '''
            options = {'verbose': 1}

        result: Optional[opt.OptimizeResult] = None
        n_est: int = self.estimate_n(lc)
        n_add: Final[int] = 0  # TODO: Run test cases to see if trying additional ply groups yields significant weight savings and is worth the time penalty.
        n: int = n_est
        print('Solver output:')
        # Repeat until successful and try at least `n_add` additional ply groups:
        while n <= n_est + n_add:
            # Optimize for minimum weight:
            '''
            Why I chose SciPy's 'trust-constr' method for now:
            
            Method trust-constr is a trust-region algorithm for constrained optimization. It switches between two implementations depending on the problem definition. It is the most versatile 
            constrained minimization algorithm implemented in SciPy and the most appropriate for large-scale problems. For equality constrained problems it is an implementation of Byrd-Omojokun 
            Trust-Region SQP method [...]. When inequality constraints are imposed as well, it switches to the trust-region interior point method [...]. This interior point algorithm, in turn, 
            solves inequality constraints by introducing slack variables and solving a sequence of equality-constrained barrier problems for progressively smaller values of the barrier parameter. 
            The previously described equality constrained SQP method is used to solve the subproblems with increasing levels of accuracy as the iterate gets closer to a solution.
            
            Source: SciPy 1.14.0 API reference
            '''
            bounds: opt.Bounds = self.create_scipy_bounds(n)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                # noinspection PyTypeChecker
                result = opt.minimize(
                    fun=self.weight,
                    x0=np.array(n * len(self.constr.ply_group) * [min(max(self.constr.min_num_of_covering_plies, 1), self.constr.max_num_of_contiguous_plies) * self.constr.ply_thickness]),
                    method='trust-constr',
                    hessp=sp_aux.zero_hessp,
                    bounds=bounds,
                    constraints=self.create_scipy_constraints(n, lc),
                    options=options,
                )

            # Note: 'if result.success' does not work here because result.success will be true once the `xtol` termination condition is satisfied
            # despite constraints such as the minimum reserve factor being violated.
            if result.constr_violation == 0:
                # Store the optimized continuous thickness values:
                if self.x is None or sum(result.x) < sum(self.x):
                    self.x = result.x
            else:
                n_est += 1
            n += 1

        # Return the solver output in case it's of interest:
        return result

    def export(self) -> clt.Panel:
        """
        Exports the generic panel as is.
        :return: clt.Panel instance with optimal (continuous) ply thicknesses
        """
        if self.x is None:
            raise Exception('The generic panel must be optimized before exporting it.\n')
        return clt.Panel(
            clt.Laminate(plies=[clt.Ply(self.mat, t) for t in self.x], angles=[theta for theta in (int(len(self.x) / len(self.constr.ply_group)) * self.constr.ply_group)], symmetric=True), self.dims)

    def discretize(self, lc: Optional[clt.MembraneLoadCase] = None) -> clt.Panel:
        """
        Discretizes the ply thicknesses by first rounding the laminate thickness up to the nearest multiple of the self.constr.ply_thickness and then
        approximating the continuous thickness distribution with discrete plies.
        :param lc: The load cases can be given as an argument to ensure that the discretization does not result in a reserve factor below 1.
        :return: clt.Panel instance with almost optimal and manufacturable (discrete) ply thicknesses
        """

        if self.x is None:
            raise Exception('The generic panel must be optimized before discretizing it.\n')

        t_disc = self.constr.ply_thickness
        odd = math.ceil(2 * sum(self.x) / t_disc) % 2 == 1
        z = -sum(self.x)
        disc_bounds = []  # discrete ply boundaries
        for t_cont in self.x:
            # `round(..., 4)` is required to avoid imprecisions due to floating point arithmetic.
            # `ndigits=4` because a ply thickness of e.g. 0.125 in combination with an odd number of plies will result in z-coordinates with four decimal places.
            disc_bounds.append(round(t_disc * round((z + (0.5 * t_disc if odd else 0)) / t_disc) - (0.5 * t_disc if odd else 0), 4))
            z += t_cont
        disc_bounds.append(0.5 * t_disc if odd else 0)
        cont_angles = self.constr.ply_group * int(len(self.x) / len(self.constr.ply_group))
        disc_angles = []
        for i, theta in enumerate(cont_angles):
            disc_angles.extend(round((disc_bounds[i + 1] - disc_bounds[i]) / t_disc) * [theta])

        if lc is not None:
            # If RF < 1, add plies in the center until the discretized laminate is safe.
            # The increase in RF is mainly due to the larger z-coordinate of the outer plies, not due to the stiffness of the added ply itself.
            # Hence, it doesn't make much of a difference which angle is appended and the last angle of the ply group is chosen to avoid disorientation.
            while clt.Panel(clt.Laminate(plies=[clt.Ply(self.mat, t_disc)] * len(disc_angles), angles=disc_angles, symmetric=True, odd=odd), self.dims).RF_biaxial_shear(lc) < 1:
                disc_angles.append(0)
                odd = not odd

        return clt.Panel(clt.Laminate(plies=[clt.Ply(self.mat, t_disc)] * len(disc_angles), angles=disc_angles, symmetric=True, odd=odd), self.dims)


# For testing:

if __name__ == '__main__':
    # Start timer:
    start = time.time()

    # Read the config file:
    config1 = cp_aux.read_config('data/config.ini')

    # Initialize panel and load case:
    gp1 = fh.generic_panel_from_config(config1)
    lc1 = fh.lc_from_config(config1)
    print(gp1.constr)

    # Optimize the panel:
    gp1.optimize(lc1)

    # Print the results:
    p1_cont = gp1.export()
    print(f'\nWeight of the generic panel: {format(round(p1_cont.weight(), 3), '.3f')} kg')
    print(f'\nCombined reserve factor of the generic panel: {format(round(p1_cont.RF_biaxial_shear(lc1), 3), '.3f')}')
    print(f'\nLaminate of the generic panel:\n{p1_cont.laminate}')

    p1_disc = gp1.discretize(lc1)
    print(f'\nWeight of the discretized panel: {format(round(p1_disc.weight(), 3), '.3f')} kg')
    print(f'\nCombined reserve factor of the discretized panel: {format(round(p1_disc.RF_biaxial_shear(lc1), 3), '.3f')}')
    print(f'\nLaminate of the discretized panel:\n{p1_disc.laminate}')

    # Stop timer:
    print(f'\nElapsed time: {int(round(time.time() - start, 3) * 1000)} ms')
