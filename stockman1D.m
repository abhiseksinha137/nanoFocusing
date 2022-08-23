%% Constants
Reset;
lda0=630e-9;
ldaPrime=lda0/(2*pi); % Vacuum reduced wavelength
k0=1/ldaPrime;
ed=1;                % Dielectric permittivity
em=-18+1i*0.6;      % Metal Permittivity
gamma= 0.57721;      % Euler constant
c=3e8;
L=2.5e-6;             % Length Scale
N=1000;             % Number of grid points

nFlat=sqrt(em*ed/(em+ed))
z=linspace(-L,-2e-9, N);

z1=0; z2=-2.5e-6; x1=2e-9; x2=50e-9;
m_const=(x2-x1)/(z2-z1);  c_const=x1-m_const*z1;
R= (m_const*z +c_const);

% n=1./(k0*R)* sqrt(-2*ed/em)/(log(sqrt(-4*em/ed))-gamma);
n= 1./(k0.*R)*1/sqrt((-em/(2*ed))*(log(sqrt(-4*em/ed))-gamma));

km=sqrt(n.^2-em);
kd=sqrt(n.^2-ed);

dR=diff(R)./diff(z);
delta=abs(dR*sqrt(-em/(2*ed))*(log(sqrt(-4*em/ed))-gamma));

% plot(z/ldaPrime,1./n)
figure()

plot(z(1:end-1),log(delta)); % ylim([0, 0.1]); title('\delta')
hold on;

plot(z, log((1./n))); xlabel('x'); title('1/(n)')


%% Ez
r=R+2e-9;
B=besseli(0,k0*km.*R)./besselh(0,k0*km.*R) ;
Ez=THETA(R-r).*besseli(0, k0*km.*R) + THETA(r-R).*B.*besselh(0,k0*kd.*R);
% Ez=B.*besselh(0,k0*kd.*r);

figure();
subplot(1,3,1)
plot(z/ldaPrime ,abs(Ez), 'k', 'LineWidth',1.5)
xlabel('z'); ylabel('Ez'); title('Ez')

%% %% Ex
Ex= THETA(R-r)*1i.*n./km.*besseli(1,k0*km.*r) +THETA(r-R)*1i.* n./kd .*B.* besselh(1,k0*kd.*r);

subplot(1,3,2)
plot(z/ldaPrime ,abs(Ex), 'k', 'LineWidth',1.5)
xlabel('z'); ylabel('Ex'); title('Ex')


%% Intesity

subplot(1,3,3)
plot(z/ldaPrime ,abs(Ex).^2 + abs(Ez).^2 , 'k', 'LineWidth',1.5)
xlabel('z'); ylabel('intesity'); title('Intensity')
set(gcf, 'Position', [46         404        1256         175])
exportgraphics(gcf,'1D.png', Resolution=300);







