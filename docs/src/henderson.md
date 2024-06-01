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

## Projecting from the tangent space

A point on $\mathcal M$ corresponding to a vector $s$ in the tangent space is solution to:

$$\begin{aligned}
F\left(u_i(s)\right) & =0 \\
\Phi_i^T\left(u_i(s)-\left(u_i+\Phi_i s\right)\right) & =0.
\end{aligned}$$

Hence, $u_i+\Phi_i s$ is projected on $\mathcal M$ at the point $u_i$.
 
## Algorithm

The algorithm [`Henderson`](@ref) described in [^Henderson] covers a manifold $\mathcal M$ with local approximations of the manifold which are called **charts**. It starts with a point $u_0$ on $\mathcal M$, a k-d ball $B_{R_0}(0)$ on the tangent space at $u_0$, a convex polygon $\mathbb P$ on the tangent space. 

A chart $C$ is defined as $C := (u_0,B_{R_0}(0), \mathbb P)$. $R_0$ is the radius of validity of the chart meaning that the distance of the ball (on the tangent) space to the manifold $\mathcal M$ is less than a prescribed value $\epsilon$.

The algorithm selects a set of charts $(C_i)_i$ which covers $\mathcal M$ meaning that the union of the projection of the balls $B_{R_i}(0)$ on the manifold covers (parts of) the manifold.

Charts intersect if the projection of their validity balls on the manifold intersect. Their polygons are then trimmed in order to remove the intersection which can move some vertices inside the validity ball.

The algorithm ends when all the charts have their polygon inside the validity ball.

In the case this is not the case, the algorithm selects such a chart $C = (u_0,B_{R_0}(0), \mathbb P)$ and define a new chart $C_{new}$ by using an exterior vertex $s$ of $\mathbb P$. $s$ is projected on the manifold $\mathcal M$, its tangent space is computed, its radius is $R_0$ and the new polygon is the cube. The new chart $C_{new}$ is then intersected with the previously computed charts.




## References

[^Henderson]:> Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451–76. https://doi.org/10.1142/S0218127402004498.

[^Dankowicz]:> Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013. https://doi.org/10.1137/1.9781611972573.
