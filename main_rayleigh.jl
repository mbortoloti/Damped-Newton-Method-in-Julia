#
#    Damped Newton Method for solving 
#    The Rayleigh Quotienthe on the Sphere
#

using Random;
using LinearAlgebra;

# Dimension setting
    n = 2;


# Symmetric positive definite matrix setting
    rng = MersenneTwister(1234);
    A = randn(rng,n,n);
    A = A' * A;
    A = A + n * Matrix{Float64}(I, n, n);
    
# Rayleigh quotient definition
    f(x) = - x' * A * x;

# Euclidean gradient of f
    ∇f(x) = -2.0 * A * x;

# Orthogonal projection function
    proj(u,x) = (I - x * x') * u;

# Riemannian gradient of f
    gradf(x) = proj(∇f(x),x);



