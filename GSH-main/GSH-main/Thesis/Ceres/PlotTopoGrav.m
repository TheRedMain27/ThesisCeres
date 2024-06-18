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

radius = 470e3; % (park, 2016)
gravitationalParameter = 62.62736e9; % (konopliv, 2018)
SHbounds = [1 10];

topography = matfile("topography.mat").topography;
gravity = readGravity(gravityFilename, SHbounds, gravitationalParameter,...
    radius);

G = 6.6743e-11;

bouguer = bouguerCorrection(topography, 910, G);

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
imagesc(topography / 1e3)
axis image
bar = colorbar;
bar.Label.String = "Topography [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/Topography")
saveas(gcf, "Images/PNG/Topography.png")

figure('Position', figurePosition)
colormap("turbo")
imagesc(gravity * 1e5)
axis image
bar = colorbar;
bar.Label.String = "Gravitational Acceleration [mGal]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/Gravity")
saveas(gcf, "Images/PNG/Gravity.png")
% 
% figure('Position', figurePosition)
% colormap("turbo")
% imagesc((gravity - bouguer) * 1e5)
% axis image
% bar = colorbar;
% bar.Label.String = "Bouguer Anomaly [mGal]";
% xticks(longitudeTicks);
% xticklabels(longitudeTickLabels);
% yticks(latitudeTicks);
% yticklabels(latitudeTickLabels);
% savefig("Images/BouguerAnomaly")
% saveas(gcf, "Images/PNG/BouguerAnomaly.png")