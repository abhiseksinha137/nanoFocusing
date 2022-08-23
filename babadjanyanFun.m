function[E_TH_Max, E_R_Max]=babadjanyanFun(r0)
alpha=0.1/2;
e1=1;
e2=-18+1i*0.6;
lda0=0.633e-6;
c=3e8;
w=2*pi*c/lda0;
eta=5.2+1i*0.12;


theta=linspace(0,pi/2,1000);
r=linspace(6e-6,1e-9,1000);

[R,TH]=meshgrid(r,theta);

%% Field Calculations
E_theta=1i*c/(w*e1) ./R.^(3/2) * (eta+1i/2).* besselk(1, eta*TH).*...
    exp(1i*eta*log(R/r0));

E_theta(TH<alpha)=nan;


E_r=1i*c/(w*e1)*eta./R.^(3/2) .* besselk(0, eta*TH) .* exp(1i*eta*log(R/r0));
E_r(TH<alpha)=nan;

E_TH_Max=max(max(abs(E_theta)));
E_R_Max=max(max(abs(E_r)));




