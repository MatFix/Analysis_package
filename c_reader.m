
clear variables
close all

%% Simulation parameters

Nx = 128;
Ny = 128;
c = 3.90625e-9;                          % cell size

k = -1;                            % handedness

init = 1;
% files renaming
rename = true;

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\NUOVE\elastic\static\normal\';
simulationFolder = 'magnetoelastic_static_rings_500nm_ME0e-6\';

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

%% Loop

N = Nx/2;

theta = [];
phi = [];

r = zeros(N,1);

for kk = 1:N
    
    thetaTemp = [];
    phiTemp = [];

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
    MmatX = rot90(MmatX,3);
    
    MmatY = reshape(My,[Nx,Ny]);
    MmatY(MmatY == 0) = NaN;
    %%
    MmatY = rot90(MmatY,3);
     
    %%
    clear Mx My M
    
    for ii = 1:Nx
        for jj = 1:Ny
            if (abs(MmatX(ii,jj)) > 1e-10) && (abs(MmatY(ii,jj)) > 1e-10)
                mx = MmatX(ii,jj);
                my = MmatY(ii,jj);
                dx =  ii - (Nx + 1)/2;
                dy =  jj - (Ny + 1)/2;
                phiTemp = [phiTemp; atan2d(my,mx)];
                thetaTemp = [thetaTemp; atan2d(dy,dx)];
            end
        end
    end
    
    r(kk) = c*sqrt(dx^2 + dy^2);
    
    if kk == 1
        NdataPoints = length(thetaTemp);
    end
    
    NcurrentPoints = length(thetaTemp);
    
    % zero padding to match all lengths
%     if length(thetaTemp) < NdataPoints
        phiTemp = [phiTemp;NaN(NdataPoints-NcurrentPoints,1)];
        thetaTemp = [thetaTemp;NaN(NdataPoints-NcurrentPoints,1)];
%     end
    theta = [theta;thetaTemp];
    phi = [phi;phiTemp];
end 
%%
thetaMatrix = reshape(theta,[],N);
phi = reshape(phi,[], N);
%%
% 
phiSorted = zeros(size(phi));
thetaSorted = zeros(size(phi));
for jj = 1:N
[phiSorted(:,jj),I] = sort(phi(:,jj));
thetaTemp = thetaMatrix(:,jj);
thetaSorted(:,jj) = thetaTemp(I);
end
% values = find(phiSorted<0);
% phiSorted(values) = phiSorted(values)+360;

% p_temp2 = phiSorted + k*90;
% p_temp2(thetaSorted >= 180) = 360 + p_temp2(thetaSorted >= 180);
% phi = p_temp2;
% 
% theta = thetaSorted;

 %% find the c fitting coefficient

baRatio = nan(N,1);
for index = 1:N
    value = (find(isnan(thetaSorted(:,index)),1) - 1)/2;
    
    if isempty(value)
        value = N/2;
    end
    
    if value < 10
        continue
    end
    
    XX = phiSorted(1:value,index);
    YY = thetaSorted(1:value,index);
    
    fittedCurve = fit(XX,YY,'atan2d(c^2 * sind(x),cosd(x))','StartPoint', [1]);
    baRatio(index) = fittedCurve.c;
end

[r,I] = sort(r);
baRatio = baRatio(I);


%%
figure
plot(r*1e9,baRatio,'linewidth',1.1)

xlabel('{\itr}  [nm]')
ylabel('b/a')
box off
title('b/a ratio as a function of distance from centre')
saveas(gcf, [folder '\ba_ratio_radius'], 'fig')


% save to file the fitted b/a ratio
stringValues = sprintf('%f     %f\r\n',[baRatio,r*1e9]');

fileID = fopen([folder '\ba_ratio.txt'],'w');
    fprintf(fileID,'B/A ratio -- r (m) \r\n');
    fprintf(fileID,stringValues);
fclose(fileID);
% 
% 
% %%
% figure
% 
% % This is to get the "true" unstressed state
% % useful in thick dots 
% phi_natural = phi(:,1);
% %
% 
% for index = init:3:N
% plot(thetaVector,phi(:,index) - phi_natural,'Marker','.','LineWidth',1.3,'Markerfacecolor','auto','MarkerSize',15);
% hold on
% end
% 
% g = 0.5;
% plot(thetaVector, zeros(length(thetaVector)),':','linewidth',1.5,'color',g*[1 1 1])
% xlim([min(thetaVector) max(thetaVector)])
% xlabel('\theta  [°]')
% ylabel('\Delta\phi  [°]')
% box off
% lgd = legend(str,'location','southeast');
% legend('boxoff')
% hold off
% title('Deviation in spin angle from the unstressed state') 
% saveas(gcf, [folder '\deltaphi'], 'fig')
   