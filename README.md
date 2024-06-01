# MultiParamContinuation

Perform covering of an immersed defined manifold $F(u)=0$ where

$$F:\mathbb R^{n} \to \mathbb R^m. $$

> [!WARNING]
> For now, only 2d immersed manifold are handled, ie n = m+2

## Example

![](https://github.com/rveltz/MultiParamContinuation.jl/blob/main/examples/torus.png?raw=true)

## Support and citation
If you use this package in your work, we ask that you cite the following paper. Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future. It is referenced on HAL-Inria as follows:

```
@misc{veltz:hal-02902346,
  TITLE = {{BifurcationKit.jl}},
  AUTHOR = {Veltz, Romain},
  URL = {https://hal.archives-ouvertes.fr/hal-02902346},
  INSTITUTION = {{Inria Sophia-Antipolis}},
  YEAR = {2020},
  MONTH = Jul,
  KEYWORDS = {pseudo-arclength-continuation ; periodic-orbits ; floquet ; gpu ; bifurcation-diagram ; deflation ; newton-krylov},
  PDF = {https://hal.archives-ouvertes.fr/hal-02902346/file/354c9fb0d148262405609eed2cb7927818706f1f.tar.gz},
  HAL_ID = {hal-02902346},
  HAL_VERSION = {v1},
}
```


## Related:

- [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)

## Acknowledgment

Romain Veltz wishes to thank Mike Henderson for fruitful discussion.

## References

[1] Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451–76. https://doi.org/10.1142/S0218127402004498.

[2] Henry, Damennick B., and Daniel J. Scheeres. “Fully Numerical Computation of Heteroclinic Connection Families in the Spatial Three-Body Problem.” Communications in Nonlinear Science and Numerical Simulation 130 (March 2024): 107780. https://doi.org/10.1016/j.cnsns.2023.107780.

[3] Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013. https://doi.org/10.1137/1.9781611972573.