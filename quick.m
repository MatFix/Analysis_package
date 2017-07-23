x = linspace(-Ny*c/2,Ny*c/2,Ny);
y = Mmat(32,:);

Xqp = linspace(-Ny*c/2,Ny*c/2,10*Ny);
Yqp = spline(x,y,Xqp);

str = int2str(Ny*c*1e9);
plot(Xqp*1e9,Yqp,'k','linewidth',1.3)
save(['vars_' str], 'Xqp', 'Yqp')
