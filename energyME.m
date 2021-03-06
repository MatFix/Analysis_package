%% Vortex-antivortex magnetization state peak analyzer
%
%  peakFinder reads from the A_ii.ovf files generated by mumax3 a matrix
%  corresponding to the mz values; performs a five-point 2D Gaussian least
%  squares fit to locate the peak, and plots its position in time.
%
%  To help the program locate a peak when there's "turbulence" in mz, the
%  xSlicel, xSlicer, ySliceu, ySliced parameters can be used to exclude
%  regions containing unwanted peaks.
%
%  For accessory functions, see also: GAUSSIANPEAK RENAMEFILES
%

clear variables
close all

%% Simulation parameters

% Nx = 128;
% Ny = 128;
c = 1e-9;                          % cell size
cz = 5e-9;
nz = 4;

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
stills = 'n';

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\NUOVE\elastic\static\normal\';
simulationFolder = 'magnetoelastic_energy_320nm_prova_2\tot\';

folder = [dailyFolder simulationFolder];            % folder  python 3.5 or
% above)

if rename
PythonScript = 'batchRenamer.py';                   % Python rename script
    try
% calls the MATLAB wrapper to a pyhton function (requires
        renameFiles(folder,PythonScript);
    catch
        disp('Error: files already renamed. Continuing execution...')
    end
end

%% Extract data from table

fid = fopen([folder 'table.txt']);
fgets(fid);                                     % skip first line
matrix = cell2mat(textscan(fid, '%f %f %f %f %f %*[^\n]'));   % read time column, ignore rest
Ku = matrix(:,5);
Nfiles = length(Ku);
fclose(fid);

%%
%read the size of matrix
% files reading
fid = fopen([folder 'A_1.ovf']);

% read the M columns
Evec = cell2mat(textscan(fid,'%f%*[^\n]','CommentStyle','#'));
fclose(fid);

Evec = Evec(1:(length(Evec)/nz));

N = sqrt(length(Evec));
clear Evec

x = linspace(-c/2*(N-1),c/2*(N-1),N);
y = linspace(-c/2*(N-1),c/2*(N-1),N);
[X,Y] = meshgrid(x,y);

if stills == 'y'
    [Xq,Yq] = meshgrid(linspace(min(x),max(x),length(x)*2)*1e9);
end

r = linspace(0,c*(N-1)/2,N/2); % The vector of values of r to integrate

% Pre-allocate for speed
cumEnergy = zeros(size(r)); % This will store the volumes for each r

% main loop
for ii = 1:Nfiles
    % files reading
    fid = fopen([folder 'A_' num2str(ii) '.ovf']);
    
    % read the M columns
    Evec = cell2mat(textscan(fid,'%f%*[^\n]','CommentStyle','#'));
    fclose(fid);
    
    Evec = Evec(1:(length(Evec)/nz));
    
    N = sqrt(length(Evec));
    
    % put Mz into matrix form
    Ed = reshape(Evec,[N,N]);
    Ed = rot90(Ed,3);
    %%
%     EdTensor(:,:,ii) = Ed;
    %%
    clear Evec
    %%
    
    % For each element in rVolume, create a filter function that is zero outside
    % of the radius of interest. Then multiply the original data by the filter before
    % integration. This way, the data can be easily integrated over the cylindrical
    % region of interest.
    for jj=1:length(r)
        % Create a filter for the original data
        filter = sqrt(X.^2+Y.^2) <= r(jj); % The filter is 1 inside radius, 0 outside
        % imagesc(x,y,filter);
        cumEnergy(jj,ii)=4*cz*c^2*trapz(trapz(Ed.*filter));
    end
    
    energyCirc(:,ii) = gradient(cumEnergy(:,ii),r);
    
    if stills == 'y'
%         surf(Xq(1,:),Yq(:,1),griddata(X*1e9,Y*1e9,Ed,Xq,Yq,'cubic'),'EdgeColor','none');
%         hold on
        imagesc(Xq(1,:),Yq(:,1),griddata(X*1e9,Y*1e9,Ed,Xq,Yq,'cubic'))
        box off
%         axis equal
        colormap jet
        colorbar
        caxis([0 10]*1e4);
        xlim([min(x) max(x)]*1e9)
        ylim([min(y) max(y)]*1e9)
        grid off
        xlabel('x [nm]')
        ylabel('y [nm]')
        title('Total energy density (interpolated) [J \cdot m^{-3}]')
        print([folder 'A_' num2str(ii)],'-r300','-dpng')
        close(figure)
    end
end

%%

figure
for ii = 1:1:Nfiles
plot(r,cumEnergy(:,ii),'linewidth',1.4);
% Plot the cumulative integrated area vs radius
hold on
end
box off
xlabel('r [nm]')
ylabel('Cumulative energy [J]')
hold off
saveas(gcf, [folder '\Ed_tot_r'], 'fig')

figure
for ii = 1:1:Nfiles
plot(r*1e9,energyCirc(:,ii),'linewidth',1.4);
% Plot the cumulative integrated area vs radius
hold on
end
box off
xlabel('r [nm]')
ylabel('Energy density per unit radius [J \cdot m^{-1}]')
hold off
saveas(gcf, [folder '\Ed_radius'], 'fig')

figure
for ii = 1:1:Nfiles
plot(r*1e9,energyCirc(:,ii)./(2*pi*r'),'linewidth',1.4);
% Plot the cumulative integrated area vs radius
hold on
end
box off
xlabel('r [nm]')
ylabel('E [J \cdot m^{-2}]')
title('Demagnetizing energy per unit area (averaged around the circumference)')
hold off
saveas(gcf, [folder '\Ed_area'], 'fig')
%%
% 
% CC = c*ones(N,1);
% Etot = squeeze(4*cz*c^2*trapz(trapz(EdTensor)));
% 
% figure
% plot(Ku,Etot,'linewidth',1.2)
% box off
% xlabel('K_{u1} [J \cdot m^{-3}]')
% ylabel('Total E_d [J]')
% saveas(gcf, [folder '\Ed_tot'], 'fig')