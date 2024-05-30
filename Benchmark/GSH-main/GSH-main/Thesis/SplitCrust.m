clear;
close all;
clc;
clearvars;

addpath('..\Tools');

topoFileName = "..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
mohoFileName = "..\..\..\Data\MohoDepthSzwillusEtAl2019";  % (Szwillus et al., 2019)
tectonicAgesFileName = "..\..\..\Data\ThermoTectonicAgeGoutorbeEtAl2011"; % (Goutorbe et al., 2011)

[fullCrustMap, mesozoicCenozoicMap, paleozoicMap, proterozoicMap, ...
    archeanMap] = readTectonicAges(tectonicAgesFileName);

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
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;

% % optimal values from gridsearch single crust
% baseReferenceDepth = 24e3;
% baseMantleDensity = 3150;
% baseCrustDensity = 2620;
% baseElasticThickness = 16.6e3;
% 
% depthAlterations = 8e3:8e3;
% sizeDepthAlterations = size(depthAlterations);
% mantleDensityAlterations = -400:-400;
% sizeMantleDensityAlterations = size(mantleDensityAlterations);
% crustDensityAlterations = -660:-660;
% sizeCrustDensityAlterations = size(crustDensityAlterations);
% thicknessAlterations = 62.4e3:62.4e3;
% sizeThicknessAlterations = size(thicknessAlterations);
% 
% depthIndexes = 1:sizeDepthAlterations(2);
% mantleDensityIndexes = 1:sizeMantleDensityAlterations(2);
% crustDensityIndexes = 1:sizeCrustDensityAlterations(2);
% thicknessIndexes = 1:sizeThicknessAlterations(2);
% 
% minError = inf;
% minThickness = 0;
% 
% minReferenceDepth = 0;
% minMantleDensity = 0;
% minCrustDensity = 0;
% minElasticThickness = 0;
% 
% for depthIndex = depthIndexes
%     referenceDepth = baseReferenceDepth + depthAlterations(depthIndex);
%     for mantleIndex = mantleDensityIndexes
%         mantleDensity = baseMantleDensity + ...
%             mantleDensityAlterations(mantleIndex);
%         for crustIndex = crustDensityIndexes
%             crustDensity = baseCrustDensity + ...
%                 crustDensityAlterations(crustIndex);
%             for thicknessIndex = thicknessIndexes
%                 elasticThickness = baseElasticThickness + ...
%                     thicknessAlterations(thicknessIndex);
% 
%                 crustGravity = calculateLowerCrustGravity(radius, ...
%                     referenceDepth, crustDensity, mass);
%                 airyCrust = airyEqualPressures(topography, crustDensity,...
%                     mantleDensity, referenceDepth, g0, crustGravity);
%                 flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
%                     shearModulus, elasticThickness, mantleDensity, ...
%                     crustDensity, g0, radius);
% 
%                 error = selectiveRMSE(flexureCrust, measuredThickness, ...
%                     mask);
% 
%                 if error < minError
%                     minError = error;
% 
%                     minReferenceDepth = referenceDepth;
%                     minMantleDensity = mantleDensity;
%                     minCrustDensity = crustDensity;
%                     minElasticThickness = elasticThickness;
%                 end
%             end
%         end
%     end
% end
% 
% disp("Minimum error: " +  string(minError / 1e3) + " [km]")
% disp("Reference depth @ minimum error: " + string(minReferenceDepth / 1e3) + " [km]")
% disp("Mantle density @ minimum error: " + string(minMantleDensity) + " [kg/m3]")
% disp("Crust density @ minimum error: " + string(minCrustDensity) + " [kg/m3]")
% disp("Elastic thickness @ minimum error: " + string(minElasticThickness / 1e3) + " [km]")

%% Found results:
% Reference Depth = 32 km
% Mantle Density = 2750 kg/m3
% Crust Density = 1960 kg/m3
% Elastic Thickness = 79 km
% Minimum Error = 4.527 km

referenceDepth = 32e3;
mantleDensity = 2750;
crustDensity = 1960;
elasticThickness = 79e3;

crustGravity = calculateLowerCrustGravity(radius, referenceDepth, ...
    crustDensity, mass);
airyCrust = airyEqualPressures(topography, crustDensity,mantleDensity, ...
    referenceDepth, g0, crustGravity);
flexureCrust = flexureThinShell(airyCrust, bulkModulus, shearModulus, ...
    elasticThickness, mantleDensity, crustDensity, g0, radius);

disp("Meso- and cenozoic crust error: " + ...
    string(selectiveRMSE(flexureCrust, ...
    measuredThickness, mesozoicCenozoicMap) / 1e3) + " [km]")
disp("Paleozoic crust error: " + ...
    string(selectiveRMSE(flexureCrust, measuredThickness, ...
    paleozoicMap) / 1e3) + " [km]")
disp("Proterozoic crust error: " + ...
    string(selectiveRMSE(flexureCrust, measuredThickness, ...
    proterozoicMap) / 1e3) + " [km]")
disp("Archean crust error: " + ...
    string(selectiveRMSE(flexureCrust, measuredThickness, ...
    archeanMap) / 1e3) + " [km]")
disp("Combined crust error: " + string(selectiveRMSE(flexureCrust, ...
    measuredThickness, fullCrustMap) / 1e3) + " [km]")

[error, percentageError, wrongRMS, minError, maxError] = characterizeError(...
    measuredThickness, flexureCrust);

flexureCrust(~mask) = mean(flexureCrust, "all");
error(~mask) = 0;
percentageError(~mask) = mean(percentageError, "all");

modelDepth = flexureCrust - topography + geoid;
depthError = abs(modelDepth - mohoDepth * 1e3);
mohoUncertainty = readMohoUncertainty(mohoFileName) * 1e3;
errorUncertaintyRatio = depthError ./ mohoUncertainty;

latitudeTicks = 0:30:180;
latitudeTickLabels = string(flip(-90:30:90));
longitudeTicks = 0:60:360;
longitudeTickLabels = string(-180:60:180);

figure(1)
colormap('hot');
imagesc(measuredThickness / 1e3);
bar = colorbar;
bar.Label.String = "Measured Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SplitCrust/MeasuredCrustThickness")
saveas(gcf, "Images/SplitCrust/PNG/MeasuredCrustThickness.png")

figure(2)
colormap('hot');
imagesc(flexureCrust / 1e3, "AlphaData", mask);
set(gca, 'color', [0 0 1])
bar = colorbar;
bar.Label.String = "Model Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SplitCrust/ModelCrustThickness")
saveas(gcf, "Images/SplitCrust/PNG/ModelCrustThickness.png")

figure(3)
colormap('hot');
imagesc(error / 1e3, "AlphaData", mask);
set(gca, 'color', [0 0 1])
bar = colorbar;
bar.Label.String = "Model Error [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SplitCrust/ModelError")
saveas(gcf, "Images/SplitCrust/PNG/ModelError.png")

percentageError(percentageError > 100) = 100;

figure(4)
colormap('hot');
imagesc(percentageError, "AlphaData", mask);
set(gca, 'color', [0 0 1])
bar = colorbar;
bar.Label.String = "Model Error [%]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SplitCrust/PercentageModelError")
saveas(gcf, "Images/SplitCrust/PNG/PercentageModelError.png")

errorUncertaintyRatio(errorUncertaintyRatio > 2) = 2;

figure(5)
colormap('hot');
imagesc(errorUncertaintyRatio, "AlphaData", mask);
set(gca, 'color', [0 0 1])
bar = colorbar;
bar.Label.String = "Depth Error to Uncertainty Ratio";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SplitCrust/ErrorUncertaintyRatio")
saveas(gcf, "Images/SplitCrust/PNG/ErrorUncertaintyRatio.png")