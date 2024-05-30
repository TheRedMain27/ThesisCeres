clear;
close all;
clc;
clearvars;

addpath('..\Tools');

% topoFileName = "..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
% mohoFileName = "..\..\..\Data\MohoDepthSzwillusEtAl2019";  % (Szwillus et al., 2019)
% 
% topography = readTopo(topoFileName);
% mohoDepth = readMoho(mohoFileName);
% geoid = getGeoid();
% 
% measuredThickness = mohoDepth * 1e3 - geoid + topography;
% 
% g0 = 9.80665;
% radius = 6371e3;
% mass = 3.986004418e14 / 6.6743e-11;

tic
gravityField = getEarthGravity();
toc

colormap("hot")
imagesc(gravityField)
colorbar()