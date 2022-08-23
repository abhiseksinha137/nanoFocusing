Reset;

% x=linspace(-5,5,1000);
% 
% y=THETA(4-x);

% plot(x,y, 'k', 'LineWidth',2)
% ylim([-2,2]);
% grid on;





% % Bessel
% I0=besseli(0,x);
% I1=besseli(1,x);
% plot(x,I0, 'k', x,I1,'r');


%% Constants
lda0=630e-9;
ldaPrime=lda0/(2*pi); % Vacuum reduced wavelength
k0=1/ldaPrime;
ed=1;                % Dielectric permittivity
em=-18.096+1i*0.48422;      % Metal Permittivity
gamma= 0.57721;      % Euler constant
c=3e8;
oa=0.02;             % Half opening angle
L=20e-6;             % Length Scale
N=1000;             % Number of grid points



z1=2e-9; z2=50e-9; x1=0; x2=2.5e-6;
m_const=(z1-z2)/(x1-x2);  c_const=z1-m_const*x1;
R=@(z) m_const*z -c_const;


% R=@(z) Rmin+tan(pi-oa)*z;
zList=linspace(-L,-100e-9,N); 
xList=linspace(-L,L,N);
% % plot
% plot(z,R(z), 'k', 'LineWidth',2)
% ylim([0, max(R(z))])

% The Refractive Index
nFlat=sqrt(em*ed/(em+ed))
% n=@(R, k0) 1./(k0.*R)*1/sqrt((-em/(2*ed))*(log(sqrt(-4*em/ed))-gamma));
n=@(R,k0) 1./(k0*R) * sqrt(-2*ed/em) /(log(sqrt(-4*em/ed))-gamma)
n(10e-6, k0)



% Group velocity
% vg=c/(d(n*w)/dw)
for j=1:length(zList)
    % Finding dw
    w =@(k) c*k;
    Dlda0=0.1e-9;
    Dk=-2*pi/lda0^2*Dlda0;
    Dw=w(Dk);
    %***********
    z=zList(j);
    D=(n(R(z),k0)*w(k0) - n(R(z),k0+Dk)*w(k0+Dk))/Dw;
    vg(j)=c/D;


    vp(j)=c/n(R(z),k0);
end
% figure()
% plot(zList,(vg/c))
% hold on;
% plot(zList,(vg/c))
% plot(zList, n(R(zList), k0))

%% Electric Field
kappam=@(z) sqrt(n(R(z),k0).^2-em) ;
kappad=@(z) sqrt(n(R(z),k0).^2-ed) ;


[X,Z]=meshgrid(xList,zList);
RG=abs(X);
B=besseli(0,k0*kappam(Z).*R(Z))./(besselk(0,k0*kappam(Z).*R(Z)));

% Ez= THETA(R(Z)-RG).* besseli(0, k0*kappam(Z).*RG) + THETA(RG-R(Z)).*B.*besselk(0, k0*kappad(Z).*RG);
Ez= THETA(R(Z)-RG).* besseli(0, k0*kappam(Z).*RG) + THETA(RG-R(Z)).*B.*besselk(0, k0*kappad(Z).*RG);
% Ez(isinf(Ez))=0;

% Ex=THETA(R(Z)-RG)*1i.*n(R(Z),k0)./kappam(Z).*besseli(1,k0*kappam(Z).*RG) + ...
%     THETA(RG-R(Z))*1i.*n(R(Z),k0)./kappad(Z).*B.*besselk(1,k0*kappad(Z).*RG);
Ex=THETA(R(Z)-RG)*1i.*n(R(Z),k0)./kappam(Z).*besseli(1,k0*kappam(Z).*RG) + ...
    THETA(R(Z)-RG)*1i.*n(R(Z),k0)./kappad(Z).*B.*besselk(1,k0*kappad(Z).*RG);
% Ex(isinf(Ex))=0;

surf(X,Z,log(abs(Ex)), 'LineStyle','none'); view(15,90)
set(gca, 'colorscale', 'log')
xlabel('X'); ylabel('Z')
colormap('jet')
colorbar()
title('E_x')

figure()
surf(X,Z,log(abs(Ez)), 'LineStyle','none'); view(15,90)
set(gca, 'colorscale', 'log')
xlabel('X'); ylabel('Z')
colormap('jet')
colorbar()
title('E_z')

RET=abs(Ez);

%% plot kappa
% figure()
% plot(zList, kappam(zList)./kappad(zList))
% figure()
% plot(zList, kappad(zList), '-k',zList, kappam(zList), 'r', zList, n(R(zList),k0), 'b')
% legend('\kappa_d', '\kappa_m', 'n')
% xlabel('z')
% set(gca, 'LineWidth',2)


%% Plot THETA
figure()
mesh(X,Z,real(n(R(Z),k0)));
view(0,90)
xlabel('X')
ylabel('Y');
colormap('jet')
colorbar()
