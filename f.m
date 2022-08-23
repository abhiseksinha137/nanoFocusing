function [ret] = f(n, em,ed, k0,R)
    km=sqrt(n.^2-em);
    kd=sqrt(n.^2-ed);
    
    ret = em/km .* besseli(1, k0*km*R)./besseli(0, k0*km*R) + ...
        ed/kd .* besselk(1, k0*kd*R)./besselk(0, k0*kd*R);

end