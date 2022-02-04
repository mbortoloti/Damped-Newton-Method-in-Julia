#
#   Damped Newton Method for solving Rayleigh quotient problem
#

function damped_rayleigh(x,maxiter,ϵ)

    iter = 0;
    info_error = 0;

    println("Iteration      || Grad f(x) || ");
    while true
        ng = norm(gradf(x));
        @printf("%5d      %18.14e\n",iter,ng);
        if ng < ϵ
            println("Solution was found!");
            return(x);
        end
        return(x);
    end

end