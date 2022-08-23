
function [ret] =df(n, em,ed, k0,R)
    km=sqrt(n.^2-em);
    kd=sqrt(n.^2-ed);
    DkmDn=n/em;
    DkdDn=n/ed;
    
    xm=k0*km*R;
    xd=k0*kd*R;

    ret=em*(besseli(1,xm)/besseli(0,xm)*(-1/km^2)*DkmDn + ...
        1/km * (besseli(0,xm)*Dbesseli(1,xm)- Dbesseli(0,xm)*besseli(1,xm))/(besseli(0,xm))^2 * k0*R*DkmDn) + ...
        ed*(besselk(1,xd)/besselk(0,xd)*(-1/kd^2)*DkdDn + ...
        1/kd * (besselk(0,xd)*Dbesselk(1,xd)- Dbesselk(0,xd)*besselk(1,xd))/(besselk(0,xd))^2 * k0*R*DkdDn);

end