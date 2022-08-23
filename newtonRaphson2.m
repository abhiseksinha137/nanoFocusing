function [n] = newtonRaphson2(em,ed, k0,R)
    epsilon=1e-10;
    i=1;
    n0=1.1;
    gval0=g(n0, em,ed, k0,R);
    while(i<10000)
        gval=g(n0, em,ed, k0,R);
        
        n=n0- gval./dg(n0, em,ed, k0,R);
        
        n0=n;

        error=gval/gval0;
        
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

function [ret]=g(n, em,ed, k0,R)
    km=sqrt(n^2-em);
    kd=sqrt(n^2-ed);
    ret= em/km * besseli(1, k0*km*R)/besseli(0, k0*km*R) + ...
        ed/kd * besselk(1,k0*kd*R)/besselk(0, k0*kd*R);

end

function [ret]=dg(n, em,ed, k0,R)
    dn=0.001;
    g2=g(n+dn, em,ed, k0,R);
    g1=g(n, em,ed, k0,R);
    g0=g(n-dn, em,ed, k0,R);
    
    D1=(g2-g1)/dn;
    D2=(g1-g0)/dn;
    ret=mean([D2,D1]);
end

