#
#   Damped Newton Method for solving Rayleigh quotient problem
#

function damped_rayleigh(x,maxiter,ϵ)

    iter = 0;
    info_error = 0;
    dir = "   ";
    println("Iteration      || Grad f(x) ||   Descent direction");
    while true
        ng = norm(gradf(x));
        @printf("%5d      %18.14e     %5s\n",iter,ng,dir);
        if ng < ϵ
            println("Solution was found!");
            return(iter,x,info_error);
        end
        
        iter = iter + 1;
        
        if iter > maxiter
            info_error = 1;
            println("Maximum number of iterations was achieved. Stopping...");
            return(iter,x,info_error);
        end


        # Newton equation 
        K = [hessf(x);x'];
        b = [-gradf(x);0];
        v = K \ b;

        # Gradient of the merit function
        gradphi = hessf(x) * gradf(x);

        # Verification if x belongs to that tangent space
        if abs(v'*x) > ϵ
            dir = "grd";
            v   = -gradphi;
        else
            dir = "new";
        end

        # Linesearch
        stp = armijo(x,v,gradphi);

        # Updating x
        x = ret(x,stp * v);
    end
    return(iter,x,info_error);
end


####################################################################################
#
#   Armijo linesearch
#
###################################################################################
function armijo(x,v,gradphi)

ϕ(x) = 0.5 * norm(gradf(x))^2;

stp = 1.0;
ϕx = ϕ(x);
σ = 1.e-4;
while true
    ϕy = ϕ(ret(x,stp * v));
    stptest = ϕy - ϕx - stp * σ * v' * gradphi;
    if stptest > 0
        stp = 0.5 * stp;
    else
        return(stp);
    end
end

end