Bstep = 2.5e-3;
Bmax = 150e-3;
%%

dailyFolder = 'D:\Program Files\mumax\Simulazioni\SIMULAZIONI 19-5\';
simulationFolder = 'nanodots_hysteresis_9_dots_30nm_thick\';

folder = [dailyFolder simulationFolder];            % folder containing table

fid = fopen([folder 'table.txt']);
fgets(fid);                                     % skip first line
matrix = cell2mat(textscan(fid, '%f %f %f %f %f %f %f%*[^\n]'));
fclose(fid);

mx = matrix(:,2);
B = matrix(:,5);
clear matrix

semiCourse = Bmax/Bstep;

saturation = 1:semiCourse;
downGoing = semiCourse+1:(3*semiCourse + 1);
upGoing = max(downGoing):length(B);

%% figure

semiCourse = Bmax/Bstep;

saturation = 1:semiCourse;
downGoing = semiCourse+1:(3*semiCourse + 1);
upGoing = max(downGoing):length(B);

figure

a = plot(B(saturation),mx(saturation),'k-o','LineWidth',1.5,...
'MarkerFaceColor','k','MarkerSize',4);
hold on
b = plot(B(downGoing),mx(downGoing),'b-o','LineWidth',1.5,...
'MarkerFaceColor','b','MarkerSize',4);
c = plot(B(upGoing),mx(upGoing),'r-o','LineWidth',1.5,...
'MarkerFaceColor','r','MarkerSize',4);

title('Hysteresis curve in the X direction for 1 dot (d = 200 nm, t = 30nm)')
xlabel('B  [T]')
ylabel('M_{\itx\rm}/M_{\its\rm}')

xlim([-Bmax Bmax])
ylim([-1 1])

legend([a b c],'Saturation','Down','Up')

saveas(gcf, [folder '\hysteresis'], 'fig')