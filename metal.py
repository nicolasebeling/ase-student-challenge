"""
Shortened version of metal.py in ase-project-1: Contains classes and functions related to stress analysis of metallic structures.
"""

from dataclasses import dataclass
from functools import cached_property
from typing import NamedTuple

import numpy as np


class Dimensions(NamedTuple):
    a: float
    b: float
    t: float


class IsotropicMaterial(NamedTuple):
    E: float
    G: float
    nu: float


class MembraneLoadCase(NamedTuple):
    Nx: float = 0
    Ny: float = 0
    Nxy: float = 0


@dataclass
class Panel:
    mat: IsotropicMaterial
    dimensions: Dimensions

    @cached_property
    def sigma_e(self):
        return self.mat.E * np.pi ** 2 / 12 / (1 - self.mat.nu ** 2) * (self.dimensions.t / self.dimensions.b) ** 2

    def tau_crit(self):
        alpha = self.dimensions.a / self.dimensions.b
        if alpha < 1:
            return (4 + 5.34 / alpha ** 2) * self.sigma_e
        return (5.34 + 4 / alpha ** 2) * self.sigma_e

    def sigma_crit(self, sigma_x: float, sigma_y: float, m_max: int = 3, n_max: int = 3) -> float:
        if sigma_x > 0:
            return np.inf
        k_biaxial = np.inf
        alpha = self.dimensions.a / self.dimensions.b
        beta = sigma_y / sigma_x
        for m in range(1, m_max + 1):
            for n in range(1, n_max + 1):
                if (m ** 2 + beta * n ** 2 * alpha ** 2) == 0:
                    continue
                k_biaxial_new = (m ** 2 + n ** 2 * alpha ** 2) ** 2 / alpha ** 2 / (m ** 2 + beta * n ** 2 * alpha ** 2)
                if 0 < k_biaxial_new < k_biaxial:
                    if m == m_max:
                        print('Warning: m_max in metal.Panel.sigma_x_crit_biaxial has been reached.')
                    if n == n_max:
                        print('Warning: n_max in metal.Panel.sigma_x_crit_biaxial has been reached.')
                    k_biaxial = k_biaxial_new
        return k_biaxial * self.sigma_e

    def RF_shear(self, lc: MembraneLoadCase):
        return self.tau_crit() / abs(lc.Nxy / self.dimensions.t)

    def RF_biaxial(self, lc: MembraneLoadCase, m_max: int = 3, n_max: int = 3):
        sigma_x = lc.Nx / self.dimensions.t
        sigma_y = lc.Ny / self.dimensions.t
        return self.sigma_crit(sigma_x, sigma_y, m_max, n_max) / abs(sigma_x)

    def RF_biaxial_shear(self, lc: MembraneLoadCase, m_max: int = 5, n_max: int = 3):
        return 1 / (1 / self.RF_biaxial(lc, m_max, n_max) + 1 / self.RF_shear(lc) ** 2)
