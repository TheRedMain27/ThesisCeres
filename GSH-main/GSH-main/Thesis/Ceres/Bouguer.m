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

G = 6.6743e-11;
g0 = 0.284; % assuming homogeneous density
radius = 470e3; % (park, 2016)
gravitationalParameter = 62.62736e9; % (konopliv, 2018)
mass = gravitationalParameter / G;
momentOfInertia = 0.375 * mass * radius ^ 2; % (mao and mckinnon, 2018)

directoryName = "BouguerTrialOptimization1";
if ~exist(directoryName, "dir")
    mkdir(directoryName)
end

SHbounds = [3 18];
save(directoryName + "/SHbounds.mat","SHbounds")

topography = matfile("topography.mat").topography;
gravity = readGravity(gravityFilename, SHbounds, gravitationalParameter,...
    radius);

crustThickness = 470e3;

runOptimization = false;
printResults = false;
showModel = true;

if runOptimization
    crustDensities = 910:1:1000;
    sizeCrustDensities = size(crustDensities);
    save(directoryName + "/crustDensities.mat","crustDensities")

    errors = zeros(1,sizeCrustDensities(2));
    weights = repmat(transpose(sind(0.5:179.5)), 1, 360);

    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km

    Model = struct();
        
    Model.number_of_layers = 1;
    
    % Additional variables
    Model.GM = gravitationalParameter;
    Model.Re_analyse = radius;
    Model.Re = radius;
    Model.geoid = 'none';
    Model.nmax = 18;
    Model.correct_depth = 0;

    Model.l1.bound = topography;
    Model.l2.bound = repmat(-crustThickness, 180, 360);

    for densityIndex = 1:sizeCrustDensities(2)
        crustDensity = crustDensities(densityIndex);
        
        Model.l1.dens  = crustDensity;
        
        SH_coefficients = model_SH_analysis(Model);
        data = model_SH_synthesis(lonLim, latLim, height, SHbounds, ...
            SH_coefficients, Model);
        
        bouguerAnomaly = gravity - flip(data.vec.R, 1);

        errors(densityIndex) = ...
            rmse(bouguerAnomaly(125:145,235:265), zeros(21,31), "all");
    end
    save(directoryName + "/errors.mat", "errors");
end

if printResults
    errors = matfile(directoryName + "/errors.mat").errors;
    crustDensities = ...
        matfile(directoryName + "/crustDensities.mat").crustDensities;
    
    minError = min(errors);
    index = find(errors == minError);
    crustDensity = crustDensities(index);

    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Crust density @ minimum RMSE: " + string(crustDensity) + " [kg/m3]")
    disp(" ")
end

if showModel
    crustDensity = 930;
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km

    Model = struct();
        
    Model.number_of_layers = 1;
    
    % Additional variables
    Model.GM = gravitationalParameter;
    Model.Re_analyse = radius;
    Model.Re = radius;
    Model.geoid = 'none';
    Model.nmax = 18;
    Model.correct_depth = 0;

    Model.l1.bound = topography;
    Model.l1.dens  = crustDensity;

    Model.l2.bound = repmat(-crustThickness, 180, 360);

    SH_coefficients = model_SH_analysis(Model);
    data = model_SH_synthesis(lonLim, latLim, height, SHbounds, ...
        SH_coefficients, Model);

    bouguerAnomaly = gravity - flip(data.vec.R, 1);

    latitudeTicks = 0:30:180;
    latitudeTickLabels = string(flip(-90:30:90));
    longitudeTicks = 0:60:360;
    longitudeTickLabels = string(-180:60:180);
    
    figurePosition = get(groot, 'DefaultFigurePosition');
    figurePosition(1) = figurePosition(1) - (2 * figurePosition(4) - ...
        figurePosition(3)) / 2;
    figurePosition(3) = 2 * figurePosition(4);
    
    leftx = 235;
    rightx = 265;
    uppery = 125;
    lowery = 145;

    figure('Position', figurePosition)
    colormap("turbo")
    imagesc(topography / 1e3)
    axis image
    bar = colorbar;
    bar.Label.String = "Topography [km]";
    xline(leftx)
    xline(rightx)
    yline(uppery)
    yline(lowery)
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);

    figure('Position', figurePosition)
    colormap("turbo")
    imagesc(gravity * 1e5)
    axis image
    bar = colorbar;
    bar.Label.String = "Gravitational Acceleration [mGal]";
    xline(leftx)
    xline(rightx)
    yline(uppery)
    yline(lowery)
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
    
    figure('Position', figurePosition)
    colormap("turbo")
    imagesc(flip(data.vec.R, 1) * 1e5)
    axis image
    bar = colorbar;
    bar.Label.String = "Bouguer Gravity [mGal]";
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
    
    figure('Position', figurePosition)
    colormap("turbo")
    imagesc(bouguerAnomaly * 1e5)
    axis image
    bar = colorbar;
    bar.Label.String = "Bouguer Anomaly [mGal]";
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
end