# MultiParamContinuation.jl

This Julia package aims at performing **multi parameter continuation** of possibly large dimensional equations $F(u, par) = 0$ where $F:\mathbb R^n\times \mathbb R^{par}\to\mathbb R^m$ by taking advantage of iterative methods, dense / sparse formulation and specific hardwares (*e.g.* GPU).

It incorporates a continuation algorithm [^Henderson][^Dankowicz] based on a Newton method to correct a predictor step and a Matrix-Free/Dense/Sparse eigensolver can be used to compute stability and bifurcation points.

## Limitations

As this is early work, the following limitations need to be addressed.

- It is *partially* optimized for speed (allocations, static arrays, etc)
- It is not suitable as is for large scale problems although it is very simple to address this. Note that the interface for jacobian free computation has yet to been pushed.
- It only computes 2d manifolds for now, *i.e.* $n=m+2$
- It allows loose detection of continuous events (no bisection)
- One needs to improve interface for using BVH search tree in large dimensions.

## Installation

To install it, please run

`] add MultiParamContinuation`

To install the bleeding edge version, please run

`] add MultiParamContinuation#master`

## Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on HAL-Inria as follows:

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

## Other softwares

There are many good softwares already available.

- For continuation in small dimension, there is only [COCO](https://sourceforge.net/projects/cocotools/) which is very reliable and for now more general than MultiParamContinuation.jl as it can compute k-d manifolds. Of course, we have the original C++ implementation by M. Henderson [Multifario](https://multifario.sourceforge.io).

- For large scale problems, there is only [Trilinos-LOCA](https://trilinos.github.io/nox_and_loca.html)

In Julia, there is no other alternative to the current package.

## Reproducibility


## References

[^Henderson]:> Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451–76. https://doi.org/10.1142/S0218127402004498.

[^Dankowicz]:> Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013. https://doi.org/10.1137/1.9781611972573.
