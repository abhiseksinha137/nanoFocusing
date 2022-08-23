%% Constants
Reset;
lda0=630e-9;
ldaPrime=lda0/(2*pi); % Vacuum reduced wavelength
k0=1/ldaPrime;
ed=1;                % Dielectric permittivity
em=-18.096+1i*0.48422;      % Metal Permittivity
gamma= 0.57721;      % Euler constant
c=3e8;
oa=0.02;             % Half opening angle
L=2.5e-6;             % Length Scale
N=1000;             % Number of grid points

x=linspace(-L,L,N);
z=linspace(-L,-2e-9, N);

[X,Z]=meshgrid(x,z);

z1=0; z2=-2.5e-6; x1=2e-9; x2=50e-9;
m_const=(x2-x1)/(z2-z1);  c_const=x1-m_const*z1;
R= (m_const*Z +c_const);

plot(z/1e-6,R/1e-9)






nFlat=sqrt(em*ed/(em+ed))
n= 1./(k0.*R)*1/sqrt((-em/(2*ed))*(log(sqrt(-4*em/ed))-gamma));
% plot(z,n)
% hold on;
% plot(z, ones(size(z))*nFlat)
% mesh(X,Z,(real(n))); view(0,90); colorbar()

r=abs(X);
kappam= sqrt(n.^2-em) ;
kappad= sqrt(n.^2-ed) ;

B=besseli(0,k0*kappam.*R)./(besselh(0,k0*kappam.*R));
Ez= THETA(R-r).* besseli(0, k0*kappam.*r) + THETA(r-R).*B.*besselh(0, k0*kappad.*r);

mesh(X,Z,log(abs(Ez))); view(0,90); colorbar()
figure();
mesh(X,Z,real(n)); view(0,90);




n= @(R) 1./(k0.*R)*1/sqrt((-em/(2*ed))*(log(sqrt(-4*em/ed))-gamma));
kappam=@(R) sqrt(n(R).^2-em) ;
kappad=@(R) sqrt(n(R).^2-ed) ;
%% Computing the preexponential factor
for i=1:length(z)
    rad=R(i,1);
    r1=linspace(0,rad,100);
    int1=trapz(r1,abs(besseli(1,k0*kappam(rad)*r1).^2 ).*r1);

    r2=linspace(rad,20*rad,200);
    int2=trapz(r2,abs(besselh(1,k0*kappad(rad)*r2).^2 ).*r2);

    A(i)=real(conj(n(rad)*em)/abs(kappam(rad))^2 .* abs(besselh(0, k0*kappad(rad)*rad )).^2 *int1+ ...
        conj(n(rad)*ed)/abs(kappad(rad))^2* abs(besseli(0,k0*kappam(rad)*rad))^2  *int2    ).^(-1/2);

end
figure()
mesh(X,Z,(abs(Ez.*A))); view(0,90); colorbar()
figure()
plot(z,A/max(A), z,R(:,1)/max(R(:,1))); legend('A', 'R')
