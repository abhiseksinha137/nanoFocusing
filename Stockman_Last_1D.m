%% Constants
Reset;
lda0=630e-9;
ldaPrime=lda0/(2*pi); % Vacuum reduced wavelength
k0=1/ldaPrime;
ed=1;                % Dielectric permittivity
em=-18+1i*0.6;      % Metal Permittivity
gamma= 0.57721;      % Euler constant
c=3e8;
w0=c*k0;
L=2.5e-6;             % Length Scale
N=1000;             % Number of grid points

nFlat=sqrt(em*ed/(em+ed))
z=linspace(-L,-2e-9, N);
x=linspace(-2e-6,2e-6,N);


z1=0; z2=-2.5e-6; x1=2e-9; x2=50e-9;

% RList=abs((x2-2e-9)/z2^2 *z.^2)+2e-9;
m_const=(x2-x1)/(z2-z1);  c_const=x1-m_const*z1;
RList= (m_const*z +c_const);
figure()
plot(z,RList)

for i=1:length(RList)
    n(i)=newtonRaphson(em,ed, k0,RList(i));
    n1(i)=newtonRaphson2(em,ed, k0,RList(i));

    %% Group Velocity
    d=dnBYdk(em,ed, k0,RList(i));
    vg(i)=c/(n(i)+k0*d);
%     vg(i)=c/(n(i)+w0*d);
    
%     displayLoop(i, length(RList))
end
figure()
plot(z, real(1./n), 'k', 'LineWidth',1.5);
xlabel('z');

hold on;
plot(z, real(vg/c), '--k', 'LineWidth',1.5);

Rp=-0.02;
delta=abs(Rp/k0*diff(1./n)./diff(RList));

hold on;
plot(z(1:end-1), delta*10, 'r', 'LineWidth',1.5)

legend('v_p', '10\delta')

% formatPlot(gcf);
% exportgraphics(gcf, 'StockmanData/QuadraticTip.png', 'Resolution',600)

%% Field

for i=1:length(z)
    r=RList(i)+0.2e-9;
    R=RList(i);
    kd=sqrt(n(i)^2-ed);
    km=sqrt(n(i)^2-em);
    B=besseli(0, k0*km*R)/besselk(0, k0*km*R);
    Ex(i)= B*besselk(0, k0*kd*r);
end

figure()
plot(z, real(Ex))

