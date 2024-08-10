# Laminate Stacking Sequence Optimization

A first attempt at optimizing laminate stacking sequences with respect to stability failure (i.e. buckling of composite panels with simply supported edges under combined biaxial and shear loading).
Based on generic stacks as explained by Ntourmas et al. in [1]. Uses a SciPy solver for optimization. Blending (the more difficult part) is not performed (yet).

_Note: This code has not been sufficiently tested to say with certainty that the calculated solutions are indeed optimal (or almost optimal). In the current version, it is also not very efficient, even with respect to Python's standards, because new `clt.Panel` instances are created for every constraint evaluation, for example. This is because I reused code from other projects and may be improved in the future._

## Usage

`numpy` and `scipy` must be installed.

1. In the current version, you can only optimize one composite panel for a single load case at a time. Enter the following data in `data/config.ini` or keep the default values:
   - Prepreg properties (after curing)
     - `E1` (elastic modulus in fiber direction) \[MPa\]
     - `E2` (transverse elastic modulus) \[MPa\]
     - `G12` (in-plane shear modulus) \[MPa\]
     - `nu12` (major Poisson's ratio) \[-\]
     - `rho` (density) \[t / mm<sup>3</sup>\]
     - `t` (thickness) \[mm\]
   - Panel dimensions \[mm\]
     - `a` (direction of the major compressive load)
     - `b` (direction of the minor compressive or tensile load)
   - Membrane loads \[N / mm\]
     - `Nx` (major compressive normal load)
     - `Ny` (minor compressive or tensile normal load)
     - `Nxy` (in-plane shear load)
   - Laminate constraints
     - `symmetric` (does __not__ have an effect yet; all optimized laminates are symmetric) \[`True` or `False`\]
     - `minimum ply share` (enter 0.1 to enforce the 10% rule, for example) \[-\]
     - `ply group` (repeating sequence of ply angles in the laminates; the order matters) \[°\]
     - `allow zero-thickness plies` (does __not__ have an effect yet; it is assumed to be false) \[`True` or `False`\]
       - Be aware that allowing zero-thickness plies increases the design freedom but may also lead to larger disorientations, resulting in higher interlaminar stresses.
     - `maximum number of contiguous plies` \[-\]
       - Be aware that increasing the maximum number of contiguous plies may lead to higher interlaminar stresses.
     - `minimum number of covering plies` (covering plies will have the first ply angle specified in `ply group`) \[-\]
       - Minimum considered value: 0 if `allow zero-thickness plies` is `True`, else 1
       - Maximum considered value: `maximum number of contiguous plies`
2. Run `gs_opt.py`

## Potential Improvements

### Larger Design Space

Optimize asymmetric laminates and allow zero-thickness generic plies for more design freedom.

### Validation

Validate stacking sequences taking into account the specified constraints as well as additional measures such as the non-dimensional anisotropic coefficients introduced by Nemeth in [2].

### Strength Constraints

Take strength constraints such as 
- maximum equivalent (von Mises) stress for metals and
- Puck failure criteria for composites

into account during the optimization.

### CSV Support

Read multiple materials, load cases, dimensions and constraints from CSV files in order to
- compare different materials for a certain combination of panel dimensions and loads,
- optimize a panel for a set of different load cases,
- optimize multiple panels in a single run,
- do all of the above simultaneously.

### Visualization

Export an image or a GIF that shows how the optimizer converges to the optimal stacking sequence. I'm planning to implement this using the `callback` function of `scipy.optimize.minimize` to collect all intermediate stacking sequences, `matplotlib` to plot the sequences using, for example, a stacked bar chart or `fill_between`, and Pillow (`PIL`) to create the GIF. May be useful for educational purposes.

## References
- [1] G. Ntourmas, F. Glock, F. Daoud, G. Schuhmacher, D. Chronopoulos, and E. Özcan, “Generic stacks and application of composite rules for the detailed sizing of laminated structures,” Composite Structures, vol. 276. Elsevier BV, p. 114487, Nov. 2021. doi: 10.1016/j.compstruct.2021.114487.
- [2] Nemeth, M. P. (1986). Importance of anisotropy on buckling of compression-loaded symmetriccomposite plates. In AIAA Journal (Vol. 24, Issue 11, pp. 1831–1835). American Institute of Aeronautics and Astronautics (AIAA). https://doi.org/10.2514/3.9531