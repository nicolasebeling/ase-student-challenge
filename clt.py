"""
Contains classes and functions related to classical laminate theory (CLT).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Final, NamedTuple

import numpy as np

# The Reuter matrix converts the tensor shear strain (epsilon) to the engineering shear strain (gamma).
Reuter: Final[np.ndarray] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])


# T is the linear transformation from the material to the problem coordinate system.
def T(theta_deg: float) -> np.ndarray:
    theta_rad = np.deg2rad(theta_deg)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([[c ** 2, s ** 2, 2 * s * c], [s ** 2, c ** 2, -2 * s * c], [-s * c, s * c, c ** 2 - s ** 2]])


class Direction(Enum):
    # In the material coordinate system:
    M1 = 1
    M2 = 2
    M12 = 3
    # In the problem coordinate system:
    Px = 4
    Py = 5
    Pxy = 6


class Dimensions(NamedTuple):
    a: float
    b: float


class MembraneLoadCase(NamedTuple):
    Nx: float = 0
    Ny: float = 0
    Nxy: float = 0


@dataclass
class TransverselyIsotropicMaterial:
    E1: float
    E2: float
    G12: float
    nu12: float
    rho: float

    @cached_property
    def nu21(self) -> float:
        return self.nu12 * self.E2 / self.E1


@dataclass
class Ply:
    mat: TransverselyIsotropicMaterial
    t: float

    @cached_property
    def Q(self) -> np.ndarray:
        Q11 = self.mat.E1 / (1 - self.mat.nu12 * self.mat.nu21)
        Q12 = self.mat.nu12 * self.mat.E2 / (1 - self.mat.nu12 * self.mat.nu21)
        Q22 = self.mat.E2 / (1 - self.mat.nu12 * self.mat.nu21)
        Q66 = self.mat.G12
        return np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])

    def Q_bar(self, theta: float) -> np.ndarray:
        return np.linalg.inv(T(theta)) @ self.Q @ Reuter @ T(theta) @ np.linalg.inv(Reuter)


@dataclass
class Laminate:
    plies: list[Ply]
    angles: list[float]

    # If a laminate is symmetric, only half of the stacking sequence has to be specified.
    symmetric: bool

    # If a symmetric laminate has an odd number of plies, the last ply in plies is assumed to be the mid ply and considered only once in the ABD calculation.
    # Laminate.odd only has an effect on symmetric laminates.
    odd: bool = False

    @cached_property
    def t(self) -> float:
        if self.symmetric:
            return 2 * sum([ply.t for ply in self.plies]) - (self.plies[-1].t if self.odd else 0)
        else:
            return sum([ply.t for ply in self.plies])

    @cached_property
    def ABD(self) -> np.ndarray:
        A = np.zeros(shape=(3, 3))
        B = np.zeros(shape=(3, 3))
        D = np.zeros(shape=(3, 3))
        z: float = -self.t / 2
        if self.symmetric:
            if self.odd:
                for ply, theta in zip(self.plies[:-1], self.angles[:-1]):
                    A = A + ply.Q_bar(theta) * ply.t
                    D = D + ply.Q_bar(theta) * ((z + ply.t) ** 3 - z ** 3)
                    z += ply.t
                A *= 2
                D *= 2 / 3
                A = A + self.plies[-1].Q_bar(self.angles[-1]) * self.plies[-1].t
                D = D + self.plies[-1].Q_bar(self.angles[-1]) * ((z + self.plies[-1].t) ** 3 - z ** 3) / 3
            else:
                for ply, theta in zip(self.plies, self.angles):
                    A = A + ply.Q_bar(theta) * ply.t
                    D = D + ply.Q_bar(theta) * ((z + ply.t) ** 3 - z ** 3)
                    z += ply.t
                A *= 2
                D *= 2 / 3
        else:
            for ply, theta in zip(self.plies, self.angles):
                A = A + ply.Q_bar(theta) * ply.t
                B = B + ply.Q_bar(theta) * ((z + ply.t) ** 2 - z ** 2)
                D = D + ply.Q_bar(theta) * ((z + ply.t) ** 3 - z ** 3)
                z += ply.t
            B /= 2
            D /= 3
        # noinspection PyTypeChecker
        return np.block([[A, B], [B, D]])

    def validate(self, lam: Laminate, constr: LaminateConstraints) -> bool:
        # TODO: Implement laminate validation.
        pass

    def __str__(self) -> str:
        output: str = 'Angle [deg]\t\tThickness [mm]'
        for ply, theta in zip(self.plies, self.angles):
            output += f'\n{round(theta, 3): }\t\t\t\t{format(ply.t, '.3f')}'
        if self.symmetric:
            if self.odd:
                output += '\t(mid plane)'
            else:
                output += '\n---------------------\t(mid plane)'
            for ply, theta in zip(list(reversed(self.plies))[1 if self.odd else 0:], list(reversed(self.angles))[1 if self.odd else 0:]):
                output += f'\n{round(theta, 3): }\t\t\t\t{format(ply.t, '.3f')}'
        return f'{output}'


class LaminateConstraints(NamedTuple):
    ply_thickness: float = 0.128
    balanced: bool = True
    balance_tolerance: float = 0.1
    symmetric: bool = True
    min_ply_share: float = 0.1
    ply_group: list[float] = [45, 0, -45, 90]
    allow_zero_thickness_plies: bool = False
    min_num_of_covering_plies: int = 1
    max_num_of_contiguous_plies: int = 4

    def __str__(self) -> str:
        output = ''
        output += f'The ply thickness is {self.ply_thickness} mm.\n'
        output += f'The laminate has to be balanced with a tolerance of +/-{self.balance_tolerance} mm.\n' if self.balanced else 'The laminate does not have to be balanced.\n'
        output += 'The laminate has to be symmetric.\n' if self.symmetric else 'The laminate does not have to be symmetric.\n'
        output += f'Every ply angle must account for at least {int(self.min_ply_share * 100)}% of the stack.\n'
        output += f'The laminate is built from {self.ply_group} degree plies. ' + 'Zero-thickness plies are allowed.' if self.allow_zero_thickness_plies else 'Zero-thickness plies are not allowed.\n'
        output += f'The minimum number of {self.ply_group[0]}-degree covering plies is {self.min_num_of_covering_plies}.\n'
        output += f'The maximum number of contiguous plies is {self.max_num_of_contiguous_plies}.'
        return f'{output}\n'


@dataclass
class Panel:
    laminate: Laminate
    dimensions: Dimensions

    def sigma_crit(self, sigma_x: float, sigma_y: float, m_max: int = 3, n_max: int = 3) -> float:
        sigma_crit = np.inf
        if sigma_x > 0:
            return sigma_crit
        alpha = self.dimensions.a / self.dimensions.b
        beta = sigma_y / sigma_x
        ABD = self.laminate.ABD
        for m in range(1, m_max + 1):
            for n in range(1, n_max + 1):
                if ((m / alpha) ** 2 + beta * n ** 2) == 0:
                    continue
                sigma_crit_new = (np.pi ** 2 / self.dimensions.b ** 2 / self.laminate.t / ((m / alpha) ** 2 + beta * n ** 2) *
                                  (ABD[3, 3] * (m / alpha) ** 4 + 2 * (ABD[3, 4] + ABD[5, 5]) * (m * n / alpha) ** 2 + ABD[4, 4] * n ** 4))
                if 0 < sigma_crit_new < sigma_crit:
                    if m == m_max:
                        print('Warning: m_max in clt.Panel.sigma_x_crit_biaxial has been reached.')
                    if n == n_max:
                        print('Warning: n_max in clt.Panel.sigma_x_crit_biaxial has been reached.')
                    sigma_crit = sigma_crit_new
        return sigma_crit

    def tau_crit(self) -> float:
        ABD = self.laminate.ABD
        t = self.laminate.t
        b = self.dimensions.b
        delta = np.sqrt(ABD[3, 3] * ABD[4, 4]) / (ABD[3, 4] + 2 * ABD[5, 5])
        if delta >= 1:
            tau_crit = 4 / t / b ** 2 * ((ABD[3, 3] * ABD[4, 4] ** 3) ** 0.25 * (8.12 + 5.05 / delta))
        else:
            tau_crit = 4 / t / b ** 2 * (np.sqrt(ABD[4, 4] * (ABD[3, 4] + 2 * ABD[5, 5])) * (11.7 + 0.532 * delta + 0.938 * delta ** 2))
        return tau_crit

    def RF_biaxial(self, lc: MembraneLoadCase, m_max: int = 3, n_max: int = 3) -> float:
        sigma_x = lc.Nx / self.laminate.t
        sigma_y = lc.Ny / self.laminate.t
        return self.sigma_crit(sigma_x, sigma_y, m_max, n_max) / abs(sigma_x)

    def RF_shear(self, lc: MembraneLoadCase):
        if lc.Nxy == 0:
            return np.inf
        return self.tau_crit() / abs(lc.Nxy / self.laminate.t)

    def RF_biaxial_shear(self, lc: MembraneLoadCase, m_max: int = 5, n_max: int = 3) -> float:
        return 1 / (1 / self.RF_biaxial(lc, m_max, n_max) + 1 / self.RF_shear(lc) ** 2)

    def weight(self):
        area_specific_weight: float = 0
        for ply in self.laminate.plies:
            area_specific_weight += ply.mat.rho * ply.t
        if self.laminate.symmetric:
            area_specific_weight *= 2
            if self.laminate.odd:
                area_specific_weight -= self.laminate.plies[-1].mat.rho * self.laminate.plies[-1].t
        return area_specific_weight * self.dimensions.a * self.dimensions.b


# For testing:

if __name__ == '__main__':
    pass
