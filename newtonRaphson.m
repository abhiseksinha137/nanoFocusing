function [n] = newtonRaphson(em,ed, k0,R)
    epsilon=1e-9;
    i=1;
    n0=1.1;
    while(i<10000)
        fval=f(n0, em,ed, k0,R);
        
        n=n0- fval./df(n0, em,ed, k0,R);
        
        n0=n;

        error=fval;
        
        if (abs(error)<=epsilon)
            break;
        end
        i=i+1;
    end
    if i==10000
        disp(num2str(n))
        throw(MException('Max iteration reached', ...
        'Solution not converged not found'));
    end

end


