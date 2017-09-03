%% Spin waves analyzer

clear variables
close all

%% Simulation parameters

Nx = 128;
Ny = 128;
c = 2.5e-9;                          % cell size

% N = 201;                           % number of data files
% columns = 1;                       % number of columns of data in the files

xSlicel = 1;                       % data slicing parameters
ySliceu = 1;                       % put l,u = 0 and r,d = 1 for no slicing

xSlicer = 0;
ySliced = 0;

maxVelocity = 550;                  % max velocity of vortex core [m/s]
                                    % use to limit range of vortex search

nFit = 2;                           % number of fitting points (nFit*2 + 1)

% files renaming
rename = true;

% saving of 3D stills
stills = false;
skip = 2;                           % save one still every _skip_ files
quality = 'hq';                     % high quality ('hq') or low quality
                                    %('lq') stills
rf = 3;                             % refinement in hq stills (2 -- 5)

% T simulation yes/no
T = 'yes';

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\NUOVE\gradient+gaussianspot\';
simulationFolder = 'nanodot_320nm_thermal_150K_Gradient_cell_2.5nm\';

folder = [dailyFolder simulationFolder];            % folder containing files
PythonScript = 'batchRenamer.py';                   % Python rename script

% calls the MATLAB wrapper to a pyhton function (requires python 3.5 or
% above)

if rename
    try
        renameFiles(folder,PythonScript);
    catch
        fprintf('\nError: files already renamed. Continuing execution...\n')
    end
end

%% Check number of header lines

fid = fopen([folder 'A_1.ovf']);
tline = fgets(fid);
header = 0;
while strcmp(tline(1),'#')          % header lines in files have leading #
    tline = fgets(fid);
    header = header + 1;
end
fclose(fid);

%% Extract time from table

fid = fopen([folder 'table.txt']);
fgets(fid);                                     % skip first line
time = cell2mat(textscan(fid, '%f %*[^\n]'));   % read time column, ignore rest
N = length(time);
fclose(fid);

% cut time vector to the number of analyzed points N
if length(time)> N
    %     time(1) = [];
    time = time(1:N);
end

% calculate max interval of vortex search
if maxVelocity ~=0
    dt = mean(gradient(time));
    dMax = ceil(maxVelocity*dt/c);              % max no. of cells moved
end

%%

% preallocation
Mmatx = zeros(Nx,Ny);
endGradient = zeros(N,1);
beginGradient = zeros(N,1);

% main loop
for ii = 1:N
    % files reading
    fid = fopen([folder 'A_' num2str(ii) '.ovf']);
    
    % read the M columns
    M = cell2mat(textscan(fid,'%f %f %f','CommentStyle','#'));
    fclose(fid);
    
    % select the z coordinate
%     Mz = M(:,end);
    
    % select the x coordinate
    Mx = M(:,1);
    
    % select the y coordinate
    My = M(:,2);
    
    % select the z coordinate
    Mz = M(:,3);
    
    % put Mx into matrix form
    Mmatx = reshape(Mx,[Nx,Ny]);
    
    %%
    Mmatx = rot90(Mmatx,3);
    
    % put My into matrix form
    Mmaty = reshape(My,[Nx,Ny]);
    
    %%
    Mmaty = rot90(Mmaty,3);
    
    % put Mz into matrix form
    Mmatz = reshape(Mz,[Nx,Ny]);
    
    %%
    Mmatz = rot90(Mmatz,3);
    
    
    
    %%
    clear M Mx
%     figure
    
%     mesh(Mmatx)
%     colorbar
%     xlabel('R')
%     ylabel('C')
%     hold on


%     if ii == 1
%         constB = Mmatx(:,2);
%         constE = Mmatx(:,end-1);
%     end
%     
%     beginGradient(ii) = sum(Mmatx(:,2) - constB);
%     endGradient(ii) = sum(Mmatx(:,end-1) - constE);

    
    if ii == 1
        constBx = Mmatx(:,2);
        constEx = Mmatx(:,end-1);
        
        constBy = Mmaty(:,2);
        constEy = Mmaty(:,end-1);
        
        constBz = Mmatz(:,2);
        constEz = Mmatz(:,end-1);
    end
    
    beginGradient(ii) = sum(sqrt((Mmatx(:,2) - constBx).^2 + (Mmaty(:,2) - constBy).^2 + (Mmatz(:,2) - constBz).^2));
    endGradient(ii) = sum(sqrt((Mmatx(:,end-1) - constBx).^2 + (Mmaty(:,end-1) - constBy).^2 + (Mmatz(:,end-1) - constBz).^2));
    
    
%     clear Mmatx
    
end

%%
%     figure 
%     mesh(Mmatx)
%     colorbar
%     xlabel('R')
%     ylabel('C')
%     hold on

bG = beginGradient - mean(beginGradient);
eG = endGradient - mean(endGradient);
figure
plot(time,bG)
hold on
plot(time,eG)
legend('Beginning of gradient','End of gradient')
legend('Boxoff')

%% Fourier analysis

bG = bG - mean(bG);
eG = eG - mean(eG);

dt = mean(gradient(time));
L = 2^(nextpow2(length(bG)));
f = 1/dt*(0:(L/2 - 1))/L;

% calculate fft

Q = fft(bG,L);
P2 = (abs(Q)).^2/L;
BG = P2(1:L/2);

Q = fft(eG,L);
P2 = (abs(Q)).^2/L;
EG = P2(1:L/2);

% convert to GHz
F = f*1e-9;

figure
plot(F,BG,'linewidth',1.1)
hold on
plot(F,EG,'linewidth',1.1)
xlabel('Frequency [GHz]')
ylabel('Absolute value [arb. units]')
legend('X','Y')
xlim([min(F) max(F)])
