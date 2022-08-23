function[]=subPlotFun(fig,X,Y,Z, xl, yl, tit)

cmap='jet';

figure(fig);

subplot(1,3,1)

surf(X,Y,abs(Z), 'FaceAlpha',1, 'EdgeColor','none'); 
view(0,90);
set(gca, 'colorscale', 'log')
xlabel(xl); ylabel(yl);
colormap(cmap);colorbar();
title(['|', tit,'|'])

subplot(1,3,2)

surf(X,Y,real(Z), 'FaceAlpha',1, 'EdgeColor','none'); 
view(0,90);
set(gca, 'colorscale', 'log')
xlabel(xl); ylabel(yl);
colormap(cmap);colorbar();
title(['Re[', tit,']'])
subplot(1,3,3)

surf(X,Y,angle(Z), 'FaceAlpha',1, 'EdgeColor','none'); 
view(0,90);
xlabel(xl); ylabel(yl);
colormap(cmap);colorbar();
title([ tit,' Phase'])


set(fig,'position',[85 332 1220 327])
