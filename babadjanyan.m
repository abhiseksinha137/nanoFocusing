Reset;
savePath=['r0=6/']; mkdir(savePath);
alpha=0.1/2;
e1=1;
e2=-18+1i*0.6;
lda0=0.633e-6;
c=3e8;
w=2*pi*c/lda0;
eta=5.2+1i*0.12;

r0=6e-6;

theta=linspace(0,pi/2,1000);
r=linspace(6e-6,1e-9,1000);

[R,TH]=meshgrid(r,theta);

%% Field Calculations
E_theta=1i*c/(w*e1) ./R.^(3/2) * (eta+1i/2).* besselk(1, eta*TH).*...
    exp(1i*eta*log(R/r0));

E_theta(TH<alpha)=nan;


E_r=1i*c/(w*e1)*eta./R.^(3/2) .* besselk(0, eta*TH) .* exp(1i*eta*log(R/r0));
E_r(TH<alpha)=nan;

X=R.*cos(TH);
Y=R.*sin(TH);

%% Plottings
figETheta_R_TH=figure();
subPlotFun(figETheta_R_TH,  R,TH,E_theta, '\theta', 'R', 'E_\theta')
exportgraphics(figETheta_R_TH,[savePath,'Etheta_R-TH.png'], Resolution=300);
%%
figETheta_X_Y=figure();
subPlotFun(figETheta_X_Y,  X,Y,E_theta, 'X', 'Y', 'E_\theta')
exportgraphics(figETheta_X_Y,[savePath,'Etheta_X-Y.png'], Resolution=300);
%%
figEr_R_TH=figure();
subPlotFun(figEr_R_TH,  R,TH,E_r, '\theta', 'R', 'E_r')
exportgraphics(figEr_R_TH,[savePath,'Er_R-TH.png'], Resolution=300);
%%
figEr_X_Y=figure();
subPlotFun(figEr_X_Y,  X,Y,E_r, 'X', 'Y', 'E_r')
exportgraphics(figEr_X_Y,[savePath,'Er_X-Y.png'], Resolution=300);
%%

%% Write Max
Er_max=max(max(abs(E_r)))
Etheta_max=max(max(abs(E_theta)))
lines=['E_r_Max, ', 'E_theta_Max \n',  num2str(Er_max, '%2.2e'), ',', num2str(Etheta_max, '%2.2e')];
file=fopen([savePath,'Max.txt'],'w');
fprintf(file, lines);
fclose(file);
% Write



% 
% 
% figure();
% subplot(1,2,1)
% mesh(X,Y,real(E_theta)); view(0,90); xlabel('X'); ylabel('Y')
% colormap(cmap); colorbar();
% title('E_\theta')
% 
% 
% 
% 
% 
% subplot(1,2,2)
% mesh(X,Y,real(E_r)); view(0,90); xlabel('X'); ylabel('Y')
% colormap(cmap); colorbar();
% title('E_r')
% 
% set(gcf, 'Position', [0         135        1357         531])
% exportgraphics(gcf, 'babadjanyan.png', Resolution=600)

