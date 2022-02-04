#
#    Damped Newton Method for solving 
#    The Rayleigh Quotienthe on the Sphere
#

using Random;
using LinearAlgebra;
using Printf;

include("damped_rayleigh.jl");

# Dimension setting
    n = 2;


# Symmetric positive definite matrix setting
    rng = MersenneTwister(1234);
    A = randn(rng,n,n);
    A = A' * A;
    A = A + n * I;
    
# Rayleigh quotient definition
    f(x) = - x' * A * x;

# Euclidean gradient of f
    ∇f(x) = -2.0 * A * x;

# Orthogonal projection function
    proj(u,x) = (I - x * x') * u;

# Riemannian gradient of f
    gradf(x) = proj(∇f(x),x);

# Retraction setting
    function ret(x,v)
        xpv = x + v;
        return xpv / norm(xpv);
    end



    maxiter = 100;
    ϵ = 1.e-7;
    x0 = rand(n);

    damped_rayleigh(x0,maxiter,ϵ);

