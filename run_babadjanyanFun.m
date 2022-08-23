% Run babadjanyanFun
Reset;
r0List=linspace(10e-9,6e-6,40);
for i=1:length(r0List)
    r0=r0List(i);
    [E_TH_Max(i), E_R_Max(i)]=babadjanyanFun(r0);
    i
end
%% Plot
figure()
subplot(2,2,1)
plot(r0List,E_TH_Max, 'k', LineWidth=1.5)
xlabel('r_0'); ylabel('max(E_\theta)')
subplot(2,2,2)
plot(r0List,E_R_Max, 'k', LineWidth=1.5)
xlabel('r_0'); ylabel('max(E_r)')

subplot(2,2,3)
loglog(r0List,E_TH_Max, 'k', LineWidth=1.5)
xlabel('r_0'); ylabel('max(E_\theta)')
subplot(2,2,4)
loglog(r0List,E_R_Max, 'k', LineWidth=1.5)
xlabel('r_0'); ylabel('max(E_r)')
% Save
exportgraphics(gcf, 'E max variation with r0.png', 'Resolution',300)

