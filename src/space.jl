"""
Define a hypercube C of radius r.

C = { u ∈ Rⁿ | -r <= u <= r}

You can define one using 

```
cube = Cube(1)
```

You can test whether a point is in the cube by doing

```
cube(rand(3))
```
"""
struct Cube{T <: Real}
    radius::T
end

function (f::Cube)(u, p = nothing)
    R = f.radius
    for ui in u 
        if ~(-R <= ui <= R)
            return false
        end
    end
    return true
end


"""
$TYPEDEF

Structure to define a domain defined by a set of intervals

## Fields

$TYPEDFIELDS

## Usage

For example, to create the domain [1,2] x [3,4], you can does

```
space = ProductSpace([1,3], [2,4])
```

and test whether a point belongs to this space:

```
space(rand(2))
```
"""
struct ProductSpace{T <: AbstractVector}
    "lower parts of all intervals"
    lower::T
    "upper parts of all intervals"
    upper::T
end

function (f::ProductSpace)(u, p = nothing)
    for i in eachindex(u) 
        if ~(f.lower[i] <= u[i] <= f.upper[i])
            return false
        end
    end
    return true
end