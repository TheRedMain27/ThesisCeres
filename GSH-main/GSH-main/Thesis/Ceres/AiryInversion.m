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

referenceDepth = ...
    matfile("AiryModel/referenceDepth.mat").referenceDepth;
crustDensity = ...
    matfile("AiryModel/crustDensity.mat").crustDensity;
mantleDensity = ...
    matfile("AiryModel/mantleDensity.mat").mantleDensity;

airyCrust = matfile("AiryModel/airyCrust.mat").airyCrust;
gravityResidual = matfile("AiryModel/gravityResidual").gravityResidual;

drho = (90 / (pi * G)) * gravityResidual ./ ((airyCrust - topography) / 2);
crustDensity = crustDensity + drho;

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
imagesc(crustDensity)
axis image
bar = colorbar;
bar.Label.String = "Crust Density [kg/m3]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);