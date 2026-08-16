# Henderson covering manifold

```@contents
Pages = ["henderson.md"]
Depth = 3
```

We want to approximte 

$$\mathcal M=\left\{u \in \mathbb{R}^n \mid F(u)=0, \quad F: \mathbb{R}^n \rightarrow \mathbb{R}^{m}\right\}.$$

## Tangent space basis

We get an orthonormal basis of the tangent space at a space $u_i\in\mathcal M$ by looking for a $n\times 2$ matrix $\Phi_i$ such that 

$$\binom{F_u\left(u_i\right)}{\Phi_i^T} \Phi_i=\binom{0}{I}.$$

The Euclidean metric can be replaced by a diagonal metric $D$, see [Weights (metric)](@ref weights-metric).

## Projecting from the tangent space

A point on $\mathcal M$ corresponding to a vector $s$ in the tangent space is solution to:

$$\begin{aligned}
F\left(u_i(s)\right) & =0 \\
\Phi_i^T\left(u_i(s)-\left(u_i+\Phi_i s\right)\right) & =0.
\end{aligned}$$

Hence, $u_i+\Phi_i s$ is projected on $\mathcal M$ at the point $u_i$.

The Euclidean metric can be replaced by a diagonal metric $D$, see [Weights (metric)](@ref weights-metric).
 
## [Weights (metric)](@id weights-metric)

By default, the bordered systems of the previous sections use the Euclidean metric. It can be replaced by a diagonal metric $D = \operatorname{diag}(w)$ through the field `weights = Weight(w)` of the problem (`Weight(TrivialWeight())` by default). The tangent basis is then found by solving

$$\binom{F_u\left(u_i\right)}{\Phi_i^T D} \Phi_i=\binom{0}{I}$$

and the projection on $\mathcal M$ is found by solving $F\left(u_i(s)\right)=0$ with the constraint $\Phi_i^T D\,\left(u_i(s)-\left(u_i+\Phi_i s\right)\right)=0$. The metric also enters:

- the norm $\lVert u\rVert_D = \lVert\sqrt{w}\cdot u\rVert_2$ and the distance between charts $\sum_j w_j (u^1_j-u^2_j)^2$,
- the orthonormalization of the tangent space: $\Phi^T D \Phi = I$,
- the curvature $K = \max_{i,j}\lVert X_{ij}\rVert_D$ and the adaptive radius $R = \mathrm{radius\_factor}\,\sqrt{2\epsilon/K}$.

The convention is the following: the **smaller** $w_j$, the more the coordinate $j$ is damped in the metric. A typical choice is $w_j = 1/\mathrm{scale}_j^2$. Beware: large coefficients on a coordinate amplify its contribution to the curvature $K$ and make the radius of the charts collapse to `Rmin`.

**Example.** For a state made of a spatial profile of $N$ points followed by 2 parameters, the grid is weighted by $1/N$ (as in `MFCreateWeightedNSpace` of multifario):

```@example TUTWEIGHTS
using MultiParamContinuation

F(u, p) = [u[3]]
prob = ManifoldProblem(F, [0., 0, 0], nothing;
                        weights = MultiParamContinuation.Weight([1., 1, 1e-4]),
                        )
```

**To keep in mind**: the weights only damp the curvature in the coordinates they weight. The curvature of the manifold in the parameter directions is not reduced by weights acting on the spatial profile; in that case, one weights the parameters as well or rescales the problem.

## Algorithm

The algorithm [`Henderson`](@ref) described in [^Henderson] covers a manifold $\mathcal M$ with local approximations of the manifold which are called **charts**. It starts with a point $u_0$ on $\mathcal M$, a k-d ball $B_{R_0}(0)$ on the tangent space at $u_0$, a convex polygon $\mathbb P$ on the tangent space. 

A chart $C$ is defined as $C := (u_0,B_{R_0}(0), \mathbb P)$. $R_0$ is the radius of validity of the chart meaning that the distance of the ball (on the tangent) space to the manifold $\mathcal M$ is less than a prescribed value $\epsilon$.

The algorithm selects a set of charts $(C_i)_i$ which covers $\mathcal M$ meaning that the union of the projection of the balls $B_{R_i}(0)$ on the manifold covers (parts of) the manifold.

Charts intersect if the projection of their validity balls on the manifold intersect. Their polygons are then trimmed in order to remove the intersection which can move some vertices inside the validity ball.

The algorithm ends when all the charts have their polygon inside the validity ball.

In the case this is not the case, the algorithm selects such a chart $C = (u_0,B_{R_0}(0), \mathbb P)$ and define a new chart $C_{new}$ by using an exterior vertex $s$ of $\mathbb P$. $s$ is projected on the manifold $\mathcal M$, its tangent space is computed, its radius is $R_0$ and the new polygon is the cube. The new chart $C_{new}$ is then intersected with the previously computed charts.

## Search tree

When the projection is very quickly computed, the limiting factor of the algorithm becomes finding the neighbors of a chart. This can become slow when the atlas is large giving a quadratic $N^2$ complexity in total. This can be alleviated by using a BVH tree to obtain a $N\log N$ complexity as explained in [^Henderson]. This is implemented in `MultiParamContinuation.jl` and can be set up in the struct `Henderson`, see its associated doc.

## References

[^Henderson]:> Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451–76. https://doi.org/10.1142/S0218127402004498.

[^Dankowicz]:> Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013. https://doi.org/10.1137/1.9781611972573.
