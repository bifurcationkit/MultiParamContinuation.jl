# MultiParamContinuation

Perform covering of an immersed manifold $F(u)=0$ where

$$F:\mathbb R^{n} \to \mathbb R^m. $$

> For now, only 2d immersed manifold are handled, ie n = m+2

> The algorithm is quadratic in the number $N$ of charts. This will be improved to a $N\log N$ algorithm using a tree.

## Example

![](https://github.com/rveltz/MultiParamContinuation.jl/blob/main/examples/torus.png?raw=true)

## Related:

- [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)

## Acknowledgment

Romain Veltz wishes to thank Mike Henderson for fruitful discussion.

## References

[1] Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451–76. https://doi.org/10.1142/S0218127402004498.

[2] Henry, Damennick B., and Daniel J. Scheeres. “Fully Numerical Computation of Heteroclinic Connection Families in the Spatial Three-Body Problem.” Communications in Nonlinear Science and Numerical Simulation 130 (March 2024): 107780. https://doi.org/10.1016/j.cnsns.2023.107780.

[3] Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013. https://doi.org/10.1137/1.9781611972573.