clear;
close all;
clc;
clearvars;

addpath('..\..\Tools');
addpath('..\SharedTools');

topographyFilename = ...
    "..\..\..\..\Data\Ceres_Dawn_FC_HAMO_DTM_DLR_Global_60ppd_Oct2016.tif";
gravityFilename = "..\..\..\..\Data\Ceres18C_KonoplivEtAl2018.txt";

% % Run if topography has to be reread from the original
% getTopo(topographyFilename);

topography = matfile("topography.mat").topography;

G = 6.6743e-11;
g0 = 0.284; % assuming homogeneous density
radius = 470e3; % (park, 2016)
gravitationalParameter = 62.62736e9; % (konopliv, 2018)
mass = gravitationalParameter / G;
momentOfInertia = 0.375 * mass * radius ^ 2; % (mao and mckinnon, 2018)

SHbounds = [5 14];
[gravity, gravitySHCoefficients] = readGravity(gravityFilename, ...
    SHbounds, gravitationalParameter, radius);

[n,gravityDegreeVariance] = degreeVariance(gravitySHCoefficients);

baseCompensationDepth = 70e3;
baseReferenceDensity = 1570;

compensationDepths = 40e3:5e3:100e3;
referenceDensities = [1000:100:2000, 1570];

figure(1)
yline(1)
xticks(1:(SHbounds(2) - SHbounds(1)));
xticklabels(string(SHbounds(1):SHbounds(2)));
xlabel("Spherical Harmonics Degree [-]")
ylabel("Normalized Degree Variance [-]")
hold on

for compensationDepth = compensationDepths
    prattCrust = pratt(topography, baseReferenceDensity, compensationDepth);
    
    %% Create model (from GSH package)
    Model = struct();
    
    Model.number_of_layers = 1;
    
    % Additional variables
    Model.GM = gravitationalParameter;
    Model.Re_analyse = radius;
    Model.Re = radius;
    Model.geoid = 'none';
    Model.nmax = 18;
    Model.correct_depth = 0;
    
    
    % % Topo layer
    Model.l1.bound = topography;
    Model.l1.dens  = prattCrust;
    
    % % Mantle
    Model.l2.bound = repmat(-compensationDepth, 180, 360);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km
    
    ModelSHCoefficients = segment_2layer_model(topography, ...
                        Model.l2.bound, -compensationDepth, prattCrust, 0, ...
                        compensationDepth / 5, Model);
    
    [n,modelDegreeVariance] = degreeVariance(ModelSHCoefficients);
    modelDegreeVariance(1:2) = [0, 0];
    
    normalizedDegreeVariance = modelDegreeVariance ./ gravityDegreeVariance;
    if compensationDepth == compensationDepths(1)
        color = "green";
    elseif compensationDepth == baseCompensationDepth
        color = "red";
    else
        color = "blue";
    end
    plot(normalizedDegreeVariance((SHbounds(1)+1):(SHbounds(2)+1)), "Color", color)
end
hold off

figure(2)
yline(1)
xticks(1:(SHbounds(2) - SHbounds(1)));
xticklabels(string(SHbounds(1):SHbounds(2)));
xlabel("Spherical Harmonics Degree [-]")
ylabel("Normalized Degree Variance [-]")
hold on

for referenceDensity = referenceDensities
    prattCrust = pratt(topography, referenceDensity, baseCompensationDepth);
    
    %% Create model (from GSH package)
    Model = struct();
    
    Model.number_of_layers = 1;
    
    % Additional variables
    Model.GM = gravitationalParameter;
    Model.Re_analyse = radius;
    Model.Re = radius;
    Model.geoid = 'none';
    Model.nmax = 18;
    Model.correct_depth = 0;
    
    
    % % Topo layer
    Model.l1.bound = topography;
    Model.l1.dens  = prattCrust;
    
    % % Mantle
    Model.l2.bound = repmat(-baseCompensationDepth, 180, 360);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km
    
    ModelSHCoefficients = segment_2layer_model(topography, ...
                        Model.l2.bound, -baseCompensationDepth, prattCrust, 0, ...
                        baseCompensationDepth / 5, Model);
    
    [n,modelDegreeVariance] = degreeVariance(ModelSHCoefficients);
    modelDegreeVariance(1:2) = [0, 0];
    
    normalizedDegreeVariance = modelDegreeVariance ./ gravityDegreeVariance;
    if referenceDensity == referenceDensities(1)
        color = "green";
    elseif referenceDensity == baseReferenceDensity
        color = "red";
    else
        color = "blue";
    end
    plot(normalizedDegreeVariance((SHbounds(1)+1):(SHbounds(2)+1)), "Color", color)
end
hold off