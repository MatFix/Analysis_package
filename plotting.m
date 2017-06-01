function plotting( results )
%% Plotting function
%  receives a struct _results_ containing all the values needed
%

R = results{1};
time = results{2};
N = results{3};
V = results{4};
p = results{5};
folder = results{6};

%% Trajectory 1

x = R(:,1)'*1e9;
y = R(:,2)'*1e9;

figure(1)
Map = [0.8500,0.3250,0.0980;0,0.4470,0.7410];

colormap(Map)
col = p';
col(col == 1) = 0;

z = zeros(size(x));
surface([x;x],[y;y],[z;z],[col;col],...
'facecol','no',...
'edgecol','flat',...
'linew',2);
axis equal
grid on
xlabel('x [nm]')
ylabel('y [nm]')
hold on
a = scatter(x(1),y(1),'r','filled');
b = scatter(x(end),y(end),'b','filled');
ll = plot(0,0,'color',[0,0.4470,0.7410],'linewidth',2);
rr = plot(0,0,'color',[0.8500,0.3250,0.0980],'linewidth',2);

title('Vortex core trajectory')
legend([a b ll rr],'Startpoint','Endpoint','\itp\rm = +1','\itp\rm = -1')

phi = 0:0.001:2*pi;
maxis = 300 * cos(phi);
Maxis = 500 * sin(phi);
plot(maxis,Maxis,'--','color',[0,0.4470,0.7410],'linewidth',2)

saveas(gcf, [folder '\trajectory1'], 'fig')

%% Trajectory 2

figure(2)

plot(x(ceil(N/3):end),y(ceil(N/3):end),'linewidth',1.5)
axis equal
grid on
xlabel('x [nm]')
ylabel('y [nm]')
hold on
a = scatter(x(end),y(end),'b','filled');
legend(a,'Endpoint')

saveas(gcf, [folder '\trajectory2'], 'fig')

%% x(t), y(t)

figure(3)

[hAx,hLine1,hLine2] = plotyy(time,x,time,y);

hLine1.LineWidth = 1.1;
hLine2.LineWidth = hLine1.LineWidth;

ylim(hAx(1), [min([x y]), max([x y])*1.1])
hAx(2).YLim = hAx(1).YLim;

box off

title('Time resolved vortex core movement')
xlabel('time [s]')

ylabel(hAx(1),'x [nm]') % left y-axis
ylabel(hAx(2),'y [nm]') % right y-axis

saveas(gcf, [folder '\xtyt'], 'fig')

%% Core velocity

figure(4)

plot(time,V,'linewidth',1.1)
title('Vortex core absolute velocity')
xlabel('time [s]')
ylabel('velocity [m \cdot s^{-1}]')
box off

saveas(gcf, [folder '\velocity'], 'fig')
end

