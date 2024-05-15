clear;
close all;
clc;

addpath('..\Tools');

topoFileName = "..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
mohoFileName = "..\..\..\Data\MohoCrustalThicknessSzwillusEtAl2019";  % (Szwillus et al., 2019)

bulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
shearModulus = 81.6e9; % olivine, (Mao er al., 2015)
elasticThickness = 34e3; % PPI slides

mantleDensity = 3320; % olivine
crustDensity = 2810;
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;
compensationDepth = 20e3;

topography = readTopo(topoFileName);

% AiryCrust = airyEqualMasses(topography, crustDensity, ...
%     mantleDensity, compensationDepth);
% crustGravity = calculateLowerCrustGravity(radius, compensationDepth, ...
%     crustDensity, mass);
% AiryCrust = airyEqualPressures(topography, crustDensity, ...
%     mantleDensity, compensationDepth, g0, crustGravity);

% FlexureCrust = flexureInfinitePlate(AiryCrust, bulkModulus, ...
%     shearModulus, elasticThickness, mantleDensity, crustDensity, g0, ...
%     radius);
% FlexureCrust = flexureThinShell(AiryCrust, bulkModulus,...
%     shearModulus, elasticThickness, mantleDensity, crustDensity, g0, ...
%     radius);

mohoDepth = readMoho(mohoFileName);

geoid = getGeoid();

figure(1)
colormap('hot');
imagesc(topography);
colorbar;
hold on

figure(2)
colormap('hot');
imagesc(mohoDepth);
colorbar;

figure(3)
colormap('hot');
imagesc(geoid);
colorbar;
hold off