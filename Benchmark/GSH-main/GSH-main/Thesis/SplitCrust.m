clear;
close all;
clc;

addpath('..\Tools');

topoFileName = "..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
mohoFileName = "..\..\..\Data\MohoDepthSzwillusEtAl2019";  % (Szwillus et al., 2019)

% run if you need to regenerate the mask in the right resolution
% maskName = "..\..\..\Data\MSK2014_landtypes.1min.geod.bin";
% getMask(maskName);
mask = matfile("mask360x180.mat").mask;

topography = readTopo(topoFileName);
mohoDepth = readMoho(mohoFileName);
geoid = getGeoid();

measuredThickness = mohoDepth * 1e3 - geoid + topography;

bulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
shearModulus = 81.6e9; % olivine, (Mao er al., 2015)
elasticThickness = 34e3; % PPI slides

mantleDensity = 3320; % olivine
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;

continentalCrustDensity = 2.65e3; % (Tewari, Prasad and Kumar, 2018)
oceanicCrustDensity  = 2.9e3; % (Tewari, Prasad and Kumar, 2018)
continentalReferenceDepth = mean(measuredThickness(mask), "all");
oceanicReferenceDepth = mean(measuredThickness(~mask), "all");
meanContinentalTopo = mean(topography(mask), "all");
meanOceanicTopo = mean(topography(~mask), "all");

continentalAiryCrust = airyEqualMasses(topography - meanContinentalTopo,...
    continentalCrustDensity, mantleDensity, ...
    continentalReferenceDepth) .* mask;
oceanicAiryCrust = airyEqualMasses(topography - meanOceanicTopo, ...
    oceanicCrustDensity, mantleDensity, ...
    oceanicReferenceDepth) .* ~mask;
airyCrust = continentalAiryCrust + oceanicAiryCrust;

[error, RMS, minError, maxError] = characterizeError(measuredThickness, ...
    airyCrust);

disp("RMS: " + string(RMS / 1e3) + " [km]")
disp("Minimum Error: " + string(minError / 1e3) + " [km]")
disp("Maximum Error: " + string(maxError / 1e3) + " [km]")

figure(1)
colormap('hot');
imagesc(error / 1e3);
colorbar;