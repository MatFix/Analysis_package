clear variables
close all

global xx yy Mx My rCoeff

%% Simulation parameters

Nx = 40;
Ny = 40;
c = 5e-9;                          % cell size
t = 30e-9;

% rhoTot - average resistivity
% rhoT - transveral resistivity
% rhoL - longitudinal resistivity
% deltaRho/rhoTot - AMR ratio
%
% data from from Bogart and Atkinson (2009) and McGuire and Potter (1975)
% valid for Py with Ni80Fe20

t_nm = t*1e9;
rhoTot = 1.2935*t_nm^(-0.4716) * 1e-6;  % Ohm.m
deltaRho = 0.03*rhoTot;
rhoT = (3*rhoTot - deltaRho)/3;
rCoeff = deltaRho/rhoT;

% files renaming
rename = false;

xx = linspace(-1,1,Nx);
yy = linspace(-1,1,Ny);

%% Files folder and files rename

dailyFolder = 'D:\Program Files\mumax\Simulazioni\AMR\';
simulationFolder = 'nanodots_hysteresis_1_dot_longer\';

folder = [dailyFolder simulationFolder];            % folder containing files
PythonScript = 'batchRenamer.py';                   % Python rename script

% calls the MATLAB wrapper to a pyhton function (requires python 3.5 or
% above) which renames files according to the scheme A_ii.ovf, ii = 1, ii++

if rename
    renameFiles(folder,PythonScript);
end

%%

fid = fopen([folder 'table.txt']);
fgets(fid);                                     % skip first line
matrix = cell2mat(textscan(fid, '%f %f %f %f %f %f %f%*[^\n]'));
fclose(fid);
B = matrix(:,5);
clear matrix

N = length(B);

%%

% preallocation

R = zeros(size(B));
R_simple = zeros(size(B));

%% Create a PDE Model with a single dependent variable

numberOfPDE = 1;
pdem = createpde(numberOfPDE);

a = 0;
f = 0;
deltaV = 2;

load('variables')

g = decsg(gd,sf,ns);
geometryFromEdges(pdem,g);

applyBoundaryCondition(pdem,'Edge',5, 'u', +deltaV/2);
applyBoundaryCondition(pdem,'Edge',6, 'u', +deltaV/2);
applyBoundaryCondition(pdem,'Edge',7, 'u', -deltaV/2);
applyBoundaryCondition(pdem,'Edge',8, 'u', -deltaV/2);

generateMesh(pdem,'Hmax',0.03);

[X,Y] = meshgrid(xx);
querypoints = [X(:),Y(:)]';

% line of integration

lx = 21;

%%
n = 0;

for kk = 1:N
    fid = fopen([folder 'A_' num2str(kk) '.ovf']);
    
    % read the Mx, My columns
    Mread = cell2mat(textscan(fid,'%f %f %*[^\n]','CommentStyle','#'));
    
    fclose(fid);
    
    mx = Mread(:,1);
    my = Mread(:,2);
    
    % Mx and My are arranged so that their indices M(m,n) go like y,x:
    % mx = [1 2 3 4 5 6 7 8]'
    %
    % Mx = [5 6 7 8;
    %       1 2 3 4]
    % This ensures that the PDE results and the matrices are arranged in
    % the same way
    
    Mx = flip((reshape(mx, [Ny,Nx]))');
    My = flip((reshape(my, [Ny,Nx]))');
 
    clear mx my Mread
    
    u = pdenonlin(pdem,'cfunction_altered_prova(x,y,ux,uy)',a,f);
    
    result = createPDEResults(pdem,u);
    uintrp = interpolateSolution(result,querypoints);
    uintrp = reshape(uintrp,size(X));
    
    [Fx, Fy] = gradient(uintrp);
    
    Fx(isnan(Fx)) = 0;
    Fy(isnan(Fy)) = 0;
    
    modulus_squared = Fx(:,lx).^2 + Fy(:,lx).^2;
    
    if mean(modulus_squared) == 0
        modulus_squared = 1;
    end
    
    sigma_line = 1./(rhoT + deltaRho./modulus_squared.*(Mx(:,lx).*Fx(:,lx) + My(:,lx).*Fy(:,lx)).^2);
    
    sigma_line(isnan(sigma_line)) = 0;
    
    S = - t*trapz(Fx(:,lx).*sigma_line);
    
    R(kk) = abs(deltaV/S);
    
    % Print state to console
    fprintf_r('Step completed: %i',kk);
end

fprintf_r('reset');
fprintf('\n');

%% Plotting

% FEM mesh and structure

figure
pdemesh(pdem);
axis equal

title('FEM mesh on the structure')

% for hysteresis plotting

G = gradient(B);
goingUp   = G>=0;
goingDown = G<0;

%
figure
a = plot(B(goingUp),R(goingUp),'b-o','LineWidth',1.5,...
    'MarkerFaceColor','b','MarkerSize',4);
hold on
b = plot(B(goingDown),R(goingDown),'r-o','LineWidth',1.5,...
    'MarkerFaceColor','r','MarkerSize',4);

legend([a b], 'B ascending', 'B descending')

xlabel('B [T]')
ylabel('R [\Omega]')
title('Anisotropic magnetoresistance with current parallel to applied field')
%%
% potential distribution at B = max
figure
pdeplot(pdem,'xydata',u,'zdata',u)
colormap parula
grid on
string = sprintf('Potential distribution at B = %i mT', max(B)*1e3);
title(string)