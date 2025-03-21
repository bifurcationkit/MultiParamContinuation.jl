# MultiParamContinuation

| **Documentation** | **Build Status** | **Downloads** |
|:-----------------:|:----------------:|:-------------:|
| [![docs-dev][docs-dev-img]][docs-dev-url] |  [![Build Status](https://github.com/bifurcationkit/MultiParamContinuation.jl/workflows/CI/badge.svg)](https://github.com/bifurcationkit/MultiParamContinuation.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/bifurcationkit/MultiParamContinuation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bifurcationkit/MultiParamContinuation.jl)|  |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://bifurcationkit.github.io/MultiParamContinuation.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://bifurcationkit.github.io/MultiParamContinuation.jl/dev/


Perform covering of an immersed manifold $F(u)=0$ where

$$F:\mathbb R^{n} \to \mathbb R^m. $$

> [!WARNING]
> For now, only 2d immersed manifold are handled, ie n = m+2

## Example

![](https://github.com/rveltz/MultiParamContinuation.jl/blob/main/examples/torus.png?raw=true)

## 📦 Installation

This package requires Julia >= v1.3.0

To install it, please run
 
`] add https://github.com/bifurcationkit/MultiParamContinuation.jl`

To install the bleeding edge version, please run

`] add https://github.com/bifurcationkit/MultiParamContinuation.jl#master`

## 📚 Support and citation
If you use `BifurcationKit.jl` in your work, we ask that you cite the following paper on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib). Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future.


## Related:

- [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)

## Acknowledgment

Romain Veltz wishes to thank Mike Henderson for fruitful discussion.

## References

[1] Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451–76. https://doi.org/10.1142/S0218127402004498.

[2] Henry, Damennick B., and Daniel J. Scheeres. “Fully Numerical Computation of Heteroclinic Connection Families in the Spatial Three-Body Problem.” Communications in Nonlinear Science and Numerical Simulation 130 (March 2024): 107780. https://doi.org/10.1016/j.cnsns.2023.107780.

[3] Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013. https://doi.org/10.1137/1.9781611972573.
