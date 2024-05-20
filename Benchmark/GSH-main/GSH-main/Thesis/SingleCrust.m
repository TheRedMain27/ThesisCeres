clear;
close all;
clc;

addpath('..\Tools');

topoFileName = "..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
mohoFileName = "..\..\..\Data\MohoDepthSzwillusEtAl2019";  % (Szwillus et al., 2019)

topography = readTopo(topoFileName);
mohoDepth = readMoho(mohoFileName);
geoid = getGeoid();

measuredThickness = mohoDepth * 1e3 - geoid + topography;

bulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
shearModulus = 81.6e9; % olivine, (Mao er al., 2015)
elasticThickness = 17e3; % PPI slides

mantleDensity = 3150;
crustDensity = 2620;
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
savefig("Images/MeasuredCrustThickness")
saveas(gcf, "Images/PNG/MeasuredCrustThickness.png")

figure(2)
colormap('hot');
imagesc(flexureCrust / 1e3);
bar = colorbar;
bar.Label.String = "Model Crust Thickness [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/ModelCrustThickness")
saveas(gcf, "Images/PNG/ModelCrustThickness.png")

figure(3)
colormap('hot');
imagesc(error / 1e3);
bar = colorbar;
bar.Label.String = "Model Error [km]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/ModelError")
saveas(gcf, "Images/PNG/ModelError.png")

percentageError(percentageError > 100) = 100;

figure(4)
colormap('hot');
imagesc(percentageError);
bar = colorbar;
bar.Label.String = "Model Error [%]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("Images/PercentageModelError")
saveas(gcf, "Images/PNG/PercentageModelError.png")