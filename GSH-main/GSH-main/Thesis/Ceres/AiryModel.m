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
% gravity = readGravity(gravityFilename, 3);

G = 6.6743e-11;
g0 = 0.284; % assuming homogeneous density
radius = 470e3; % (park, 2016)
gravitationalParameter = 62.62736e9; % (konopliv, 2018)
mass = gravitationalParameter / G;
momentOfInertia = 0.375 * mass * radius ^ 2; % (mao and mckinnon, 2018)

mantleDensity = 2367; % (King, 2018)
crustDensity = 910; % from optimization
referenceDepth = 11.7e3; % from optimization

crustGravity = ...
    calculateLowerCrustGravity(radius, referenceDepth, crustDensity, mass);
airyCrust = airyEqualPressures(topography, crustDensity, ...
    mantleDensity, referenceDepth, g0, crustGravity);

%% Create model (from GSH package)
Model = struct();

Model.number_of_layers = 2;

% Additional variables
Model.GM = gravitationalParameter;
Model.Re_analyse = radius;
Model.Re = radius;
Model.geoid = 'none';
Model.nmax = 18;
Model.correct_depth = 0;


% % Topo layer
Model.l1.bound = topography;
Model.l1.dens  = crustDensity;

% % Mantle
Model.l2.bound = topography - airyCrust;
Model.l2.dens = mantleDensity;

% % Bottom
Model.l3.bound = repmat(-radius, 180, 360);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];
lonLim =    [0.5 359.5 1];
height =    0; % LAMO is 375km
SHbounds =  [3 18];

SH_coefficients = model_SH_analysis(Model);
data = model_SH_synthesis(lonLim, latLim, height, SHbounds, ...
    SH_coefficients, Model);

latitudeTicks = 0:30:180;
latitudeTickLabels = string(flip(-90:30:90));
longitudeTicks = 0:60:360;
longitudeTickLabels = string(-180:60:180);

figurePosition = get(groot, 'DefaultFigurePosition');
figurePosition(1) = figurePosition(1) - (2 * figurePosition(4) - ...
    figurePosition(3)) / 2;
figurePosition(3) = 2 * figurePosition(4);

figure('Position', figurePosition)
colormap("turbo")
imagesc(airyCrust / 1e3)
axis image
bar = colorbar;
bar.Label.String = "Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/AiryCrustThickness")
saveas(gcf, "Images/PNG/AiryCrustThickness.png")

figure('Position', figurePosition)
colormap("turbo")
imagesc(flip(data.vec.R, 1) * 1e5)
axis image
bar = colorbar;
bar.Label.String = "Model Gravitational Acceleration [mGal]";
title("Airy")
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/AiryGravity")
saveas(gcf, "Images/PNG/AiryGravity.png")