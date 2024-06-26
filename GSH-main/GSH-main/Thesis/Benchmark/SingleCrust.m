clear;
close all;
clc;
clearvars;

addpath('..\..\Tools');
addpath('..\SharedTools');

topoFileName = "..\..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
mohoFileName = "..\..\..\..\Data\MohoDepthSzwillusEtAl2019";  % (Szwillus et al., 2019)

topography = readTopo(topoFileName);
mohoDepth = readMoho(mohoFileName);
geoid = getGeoid();

measuredThickness = mohoDepth * 1e3 - geoid + topography;

bulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
shearModulus = 81.6e9; % olivine, (Mao er al., 2015)
elasticThickness = 34e3; % (Watts and Moore, 2017)

mantleDensity = 3320; % olivine mineral database
crustDensity = 2835; % (Christensen and Mooney, 1995)
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;
referenceDepth = 24e3;

% airyCrust = airyEqualMasses(topography, crustDensity, ...
%     mantleDensity, referenceDepth);
crustGravity = calculateLowerCrustGravity(radius, referenceDepth, ...
    crustDensity, mass);
airyCrust = airyEqualPressures(topography, crustDensity, ...
    mantleDensity, referenceDepth, g0, crustGravity);

% flexureCrust = flexureInfinitePlate(airyCrust, bulkModulus, ...
%     shearModulus, elasticThickness, mantleDensity, crustDensity, g0, ...
%     radius);
flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
    shearModulus, elasticThickness, mantleDensity, crustDensity, g0, ...
    radius);

[error, percentageError, RMS, minError, maxError] = characterizeError(...
    measuredThickness, flexureCrust);

disp("RMS: " + string(RMS / 1e3) + " [km]")
disp("Minimum error: " + string(minError / 1e3) + " [km]")
disp("Maximum error: " + string(maxError / 1e3) + " [km]")
disp("Negative thicknesses: " + string(min(flexureCrust, [], "all") < 0))
if min(flexureCrust, [], "all") < 0
    disp("Minimum crust thickness: " + ...
        string(min(flexureCrust, [], "all") / 1e3) + " [km]")
end

modelDepth = flexureCrust - topography + geoid;
depthError = abs(modelDepth - mohoDepth * 1e3);
mohoUncertainty = readMohoUncertainty(mohoFileName) * 1e3;
errorUncertaintyRatio = depthError ./ mohoUncertainty;

latitudeTicks = 0:30:180;
latitudeTickLabels = string(flip(-90:30:90));
longitudeTicks = 0:60:360;
longitudeTickLabels = string(-180:60:180);

figurePosition = get(groot, 'DefaultFigurePosition');
figurePosition(1) = figurePosition(1) - (2 * figurePosition(4) - ...
    figurePosition(3)) / 2;
figurePosition(3) = 2 * figurePosition(4);

figure('Position', figurePosition)
colormap('hot');
imagesc(measuredThickness / 1e3);
axis image
bar = colorbar;
bar.Label.String = "Measured Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SingleCrust/MeasuredCrustThickness")
saveas(gcf, "Images/SingleCrust/PNG/MeasuredCrustThickness.png")

figure('Position', figurePosition)
colormap('hot');
imagesc(flexureCrust / 1e3);
axis image
bar = colorbar;
bar.Label.String = "Model Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SingleCrust/ModelCrustThickness")
saveas(gcf, "Images/SingleCrust/PNG/ModelCrustThickness.png")

figure('Position', figurePosition)
colormap('hot');
imagesc(error / 1e3);
axis image
bar = colorbar;
bar.Label.String = "Model Error [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SingleCrust/ModelError")
saveas(gcf, "Images/SingleCrust/PNG/ModelError.png")

percentageError(percentageError > 100) = 100;

figure('Position', figurePosition)
colormap('hot');
imagesc(percentageError);
axis image
bar = colorbar;
bar.Label.String = "Model Error [%]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SingleCrust/PercentageModelError")
saveas(gcf, "Images/SingleCrust/PNG/PercentageModelError.png")

errorUncertaintyRatio(errorUncertaintyRatio > 2) = 2;

figure('Position', figurePosition)
colormap('hot');
imagesc(errorUncertaintyRatio);
axis image
bar = colorbar;
bar.Label.String = "Depth Error to Uncertainty Ratio";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/SingleCrust/ErrorUncertaintyRatio")
saveas(gcf, "Images/SingleCrust/PNG/ErrorUncertaintyRatio.png")