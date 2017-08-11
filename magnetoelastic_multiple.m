
clear variables
close all

%% Simulation parameters

Nx = 160;
Ny = 160;
c = 5e-9;                          % cell size

k = -1;                            % handedness

init = 1;
% files renaming
rename = true;

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\NUOVE\elastic\static\normal\';
simulationFolder = 'magnetoelastic_static_rings_800nm\';

folder = [dailyFolder simulationFolder];            % folder containing files
PythonScript = 'batchRenamer.py';                   % Python rename script

% calls the MATLAB wrapper to a pyhton function (requires python 3.5 or
% above)

if rename
    try
        renameFiles(folder,PythonScript);
    catch
        disp('Error: files already renamed. Continuing execution...')
    end
end

%% Extract Ku1 from table

fid = fopen([folder 'table.txt']);
fgets(fid);                                     % skip first line
temp = cell2mat(textscan(fid, '%f %f %f %f %f%*[^\n]'));   % read time column, ignore rest
Ku1 = temp(:,end);
N = length(Ku1);
fclose(fid);

% cut K vector to the number of analyzed points N
if length(Ku1)> N
    Ku1 = Ku1(1:N);
end

theta = [];
phi = [];


for kk = init:N
    
    % files reading
    fid = fopen([folder 'A_' num2str(kk) '.ovf']);
    
    % read the M columns
    M = cell2mat(textscan(fid,'%f %f %*[^\n]','CommentStyle','#'));
    fclose(fid);
    
    Mx = M(:,1);
    My = M(:,2);
    
    % put Mx,My into matrix form
    MmatX = reshape(Mx,[Nx,Ny]);
    MmatX(MmatX == 0) = NaN;
    %%
    MmatX = rot90(MmatX);
    
    MmatY = reshape(My,[Nx,Ny]);
    MmatY(MmatY == 0) = NaN;
    %%
    MmatY = rot90(MmatY);
     
    %%
    clear Mx My M
    
    for ii = 1:Nx
        for jj = 1:Ny
            if (abs(MmatX(ii,jj)) > 1e-10) && (abs(MmatY(ii,jj)) > 1e-10)
                mx = MmatX(ii,jj);
                my = MmatY(ii,jj);
                dx =  ii - (Nx + 1)/2;
                dy =  jj - (Ny + 1)/2;
                theta = [theta; atan2d(dy,dx)];
                phi = [phi; atan2d(my,mx)];
            end
        end
    end
end 
%%
thetaMatrix = reshape(theta,[],N);
phi = reshape(phi,[], N);

thetaMatrix = thetaMatrix + 180;
[thetaVector,I] = sort(thetaMatrix(:,1));
phi_natural = thetaVector + k*90;

for ind = init:N
    p_temp = phi(:,ind);
    p_temp2 = p_temp(I) + k*90;
    p_temp2(thetaVector >= 180) = 360 + p_temp2(thetaVector >= 180);
    phi(:,ind) = p_temp2;    
end

delta_phi = phi - phi_natural*ones(1,N);
%%

baRatio = zeros(N,1);
XX = thetaVector(1:end/2);

Ku1 = Ku1*1e-3;         % J to kJ

figure
str = [];
for index = init:3:N
plot(thetaVector,phi(:,index),'Marker','.','LineWidth',1.3,'Markerfacecolor','auto','MarkerSize',15);
hold on
str = strvcat(str,sprintf('%0.1f kJ \\cdot m^{-3}\n', Ku1(index)));
end

for index = 1:N
    YY = phi(1:end/2,index);
    fittedCurve = fit(XX,YY,'atan2d(-c^2 * cosd(x),sind(x)) + p*90','StartPoint', [1, 1]);
    baRatio(index) = fittedCurve.c;
end

xlim([min(thetaVector) max(thetaVector) + 0.1])
xlabel('\theta  [°]')
ylabel('\phi  [°]')
box off
lgd = legend(str,'location','southeast');
legend('boxoff')
hold off
title('Spin angle')
saveas(gcf, [folder '\phi'], 'fig')

% save to file the fitted b/a ratio

stringValues = sprintf('%f     %f\r\n',[baRatio,Ku1]');

fileID = fopen([folder '\ba_ratio.txt'],'w');
    fprintf(fileID,'B/A ratio -- Ku1 (kJ/m^3) \r\n');
    fprintf(fileID,stringValues);
fclose(fileID);


%%
figure

% This is to get the "true" unstressed state
% useful in thick dots 
phi_natural = phi(:,1);
%

for index = init:3:N
plot(thetaVector,phi(:,index) - phi_natural,'Marker','.','LineWidth',1.3,'Markerfacecolor','auto','MarkerSize',15);
hold on
end

g = 0.5;
plot(thetaVector, zeros(length(thetaVector)),':','linewidth',1.5,'color',g*[1 1 1])
xlim([min(thetaVector) max(thetaVector)])
xlabel('\theta  [°]')
ylabel('\Delta\phi  [°]')
box off
lgd = legend(str,'location','southeast');
legend('boxoff')
hold off
title('Deviation in spin angle from the unstressed state') 
saveas(gcf, [folder '\deltaphi'], 'fig')
   