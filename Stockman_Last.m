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
N=2000;             % Number of grid points

nFlat=sqrt(em*ed/(em+ed))
z=linspace(-L,-2e-9, N);
x=linspace(-2e-6,2e-6,N);


z1=0; z2=-2.5e-6; x1=2e-9; x2=50e-9;
m_const=(x2-x1)/(z2-z1);  c_const=x1-m_const*z1;
RList= (m_const*z +c_const);

for i=1:length(RList)
    n(i)=newtonRaphson(em,ed, k0,RList(i));
    
    clc;
    displayLoop(i, length(RList))
end

plot(z, 1./n, 'k', 'LineWidth',1.5);
xlabel('z');


Rp=-0.02;
delta=abs(Rp/k0*diff(1./n)./diff(RList));

hold on;
plot(z(1:end-1), delta*10, 'r', 'LineWidth',1.5)

legend('v_p', '10\delta')

formatPlot(gcf)
exportgraphics(gcf, 'StockmanData/velocity.png', 'Resolution',600);
%% 2D
[XG, ZG]=meshgrid(x,z);
RG=abs(XG);
z=reshape(z, [length(z),1]);
n=reshape(n, [length(n),1]);
r=reshape(RList, [length(RList),1]);

Z=repmat(z, [1, length(x)]);
N=repmat(n, [1, length(n)]);
R=repmat(r, [1, length(r)]);

km=sqrt(N.^2-em);
kd=sqrt(N.^2-ed);

figure()
surf(XG,ZG,real(N), 'LineStyle','none');view(69,60);
title('R')
xlabel('X'); ylabel('Z')
colormap('jet'); colorbar();
set(gca, 'colorscale', 'log')
formatPlot(gcf);

%% Fields
B=besseli(0,k0*km.*R)./besselk(0,k0*km.*R);
Ez=THETA(R-RG).*besseli(0, k0*km.*RG)  + THETA(RG-R).* B.*besselk(0,k0*kd.*RG);

Ex=THETA(R-RG).* 1i .* N./km .* besseli(1,k0*km.*RG) + ...
    THETA(RG-R) *1i .*N./kd .* B .* besselk(1, k0*kd.*RG);

figure()
% surf(XG,ZG,abs(Ex), 'LineStyle','none');view(90,0);
imagesc(real(Ez));
xlabel('X'); ylabel('Z')
colormap('jet'); colorbar();
set(gca, 'colorscale', 'log')
formatPlot(gcf)


%% Pre-Exponential Factor
RList=reshape(RList,[length(RList),1]);

for i=1:length(RList)
    RCone=RList(i);
    rInt=linspace(0,RCone,100);
    rInt2=linspace(RCone, L,500);
    kappam=km(i,1);
    kappad=kd(i,1);
    
    xm=k0*kappam*RCone;
    xd=k0*kappad*RCone;
    nR=n(i);
    A(i)=real(conj(nR.*em)/abs(kappam)^2 * abs(besselk(0,xd))^2 * trapz(rInt, besseli(1, k0*kappam*rInt).^2 .* rInt)...
        + conj(nR*ed)/abs(kappad)^2 * abs(besseli(0, xd))^2 * trapz(rInt2, besselk(1, k0*kappad*rInt2).^2 .*rInt2))^(-1/2);

end
A=reshape(A, [length(A),1]);
ExFinal=Ex.*A;
figure()
surf(XG,ZG,real(ExFinal), 'LineStyle','none');view(90,0);
xlabel('X'); ylabel('Z')
colormap('jet'); colorbar();
set(gca, 'colorscale', 'log')
formatPlot(gcf)
title('ExFinal')

