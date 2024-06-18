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
SHbounds =  matfile("PrattModel/SHbounds.mat").SHbounds;

[gravity, gravitySHCoefficients] = readGravity(gravityFilename, ...
    SHbounds, gravitationalParameter, radius);

compensationDepth = ...
    matfile("PrattModel/compensationDepth.mat").compensationDepth;
referenceDensity = ...
    matfile("PrattModel/referenceDensity.mat").referenceDensity;

prattCrust = matfile("PrattModel/prattCrust.mat").prattCrust;
gravityResidual = matfile("PrattModel/gravityResidual").gravityResidual;

mantleDensity = 2520;

dr = -gravityResidual .* compensationDepth ^ 2 * 180 ^ 2 ./ ...
    (G * (mantleDensity - prattCrust) * pi ^ 2 * (radius - compensationDepth) ^ 4);
crustThickness = dr + compensationDepth;

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
imagesc(crustThickness / 1e3)
axis image
bar = colorbar;
bar.Label.String = "Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);