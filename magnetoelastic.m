
clear variables
close all

%% Simulation parameters

Nx = 64;
Ny = 64;
c = 5e-9;                          % cell size

N = 2;                           % number of data files

% files renaming
rename = false;

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\SIMULAZIONI 25-5\';
simulationFolder = 'nanodots_magnetoelastic_coupling_1\';

folder = [dailyFolder simulationFolder];            % folder containing files
PythonScript = 'batchRenamer.py';                   % Python rename script

% calls the MATLAB wrapper to a pyhton function (requires python 3.5 or
% above)

if rename
    renameFiles(folder,PythonScript);
end

%%

% preallocation
phi = zeros(Nx*Ny,1);
phi2 = zeros(Nx*Ny,1);
theta = zeros(Nx*Ny,1);

for kk = 1:N
    fid = fopen([folder 'A_' num2str(kk) '.ovf']);
    
    % read the Mx, My columns
    M = cell2mat(textscan(fid,'%f %f %*[^\n]','CommentStyle','#'));
    
    fclose(fid);
    
    mx = M(:,1);
    my = M(:,2);
    
    clear M
    
    for ii = 1:length(mx)
        if (abs(mx(ii)) > 1e-10) && (abs(my(ii)) > 1e-10)
            
            % current cell index in the (x,y) matrix
            
            cellx = mod(ii - 1,Nx) + 1;
            celly = floor((ii - 1)/Nx) + 1;
            
            % calculate distance, correct for half cell
            
            lx = cellx - (Nx/2);
            lx = lx - sign(lx)*1/2;
            ly = celly - (Ny/2);
            ly = ly - sign(ly)*1/2;
            
            % angles
            
            phi(ii,kk) = atan2(my(ii),mx(ii));
            
            theta(ii,kk) = atan2(ly,lx);
            
        else
            phi(ii,kk) = NaN;
        end
    end
end

theta(not(any(phi,2)),:) = [];
phi(not(any(phi,2)),:)   = [];


% zz = theta(:,1) <  0;
% 
% theta(zz,:) = theta(zz,:) + 2*pi;

theta = rad2deg(theta);
phi = rad2deg(phi) - theta;

[theta_sorted,I] = sort(theta);

phi_sorted = phi(I);

%% figures

figure
scatter(theta_sorted(:,1),phi_sorted(:,1))
figure
scatter(theta_sorted(:,2),phi_sorted(:,2))