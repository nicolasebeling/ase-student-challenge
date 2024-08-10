"""
Contains auxiliary functions meant to be used with SciPy.
"""

from typing import Callable

import numpy as np
import scipy.optimize as opt


def grad_within_bounds(bounds: opt.Bounds, f: Callable[[np.ndarray], float], eps: float = 1e-8) -> Callable[[np.ndarray], np.ndarray]:
    """
    Returns a function computing the gradient of a function `f` as a function of `x` without exceeding certain bounds.
    Useful if `f` is not computable outside these bounds, due to negative arguments being passed to numpy.sqrt for example.
    :param bounds: bounds as an instance of scipy.optimize.Bounds
    :param f: function `f(x)` returning a float
    :param eps: gradient step size (default: 1e-8)
    :return: gradient of `f` as a function of `x`
    """

    def grad(x):
        g = np.zeros(len(x))
        for k in range(len(x)):
            d = np.zeros(len(x))
            d[k] = eps if x[k] + eps <= bounds.ub[k] else -eps
            g[k] = (f(x + d) - f(x)) / d[k]
        return g

    return grad


def zero_hessp(x: np.ndarray, _: np.ndarray) -> np.ndarray:
    """
    Can be handed to e.g. scipy.optimize.minimize when minimizing linear functions to speed up computations and prevent the optimizer from exceeding the bounds when evaluating the hessian.
    :param x: parameters
    :param _: vector `p` (not used)
    :return: a numpy.ndarray the length of `x` filled with zeros
    """
    return np.zeros(len(x))
