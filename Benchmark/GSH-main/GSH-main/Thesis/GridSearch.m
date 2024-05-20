clear;
close all;
clc;
clearvars;

addpath('..\Tools');

topoFileName = "..\..\..\Data\Earth2014.BED2014.5min.geod.bin";  % (Hirt and Rexer, 2015)
mohoFileName = "..\..\..\Data\MohoDepthSzwillusEtAl2019";  % (Szwillus et al., 2019)

topography = readTopo(topoFileName);
mohoDepth = readMoho(mohoFileName);
geoid = getGeoid();

measuredThickness = mohoDepth * 1e3 - geoid + topography;

bulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
shearModulus = 81.6e9; % olivine, (Mao er al., 2015)
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;
referenceDepth = 24e3;

% optimal values from SingleVariableOptimization
baseMantleDensity = 3300;
baseCrustDensity = 2780;
baseElasticThickness = 17e3;

% depthAlterations = -100:50:100;
% sizeDepthAlterations = size(depthAlterations);
densityAlterations = -150:5:-150;
sizeDensityAlterations = size(densityAlterations);
thicknessAlterations = -400:50:-400;
sizeThicknessAlterations = size(thicknessAlterations);

% depthIndexes = 1:sizeDepthAlterations(2);
densityIndexes = 1:sizeDensityAlterations(2);
thicknessIndexes = 1:sizeThicknessAlterations(2);

minError = inf;
minThickness = 0;

% minReferenceDepth = 0;
minMantleDensity = 0;
minCrustDensity = 0;
minElasticThickness = 0;

% for depthIndex = depthIndexes
%     referenceDepth = baseReferenceDepth + depthAlterations(depthIndex);
for mantleIndex = densityIndexes
    mantleDensity = baseMantleDensity + ...
        densityAlterations(mantleIndex);
    for crustIndex = densityIndexes
        crustDensity = baseCrustDensity + ...
            densityAlterations(crustIndex);
        for thicknessIndex = thicknessIndexes
            elasticThickness = baseElasticThickness + ...
                thicknessAlterations(thicknessIndex);

            crustGravity = calculateLowerCrustGravity(radius, ...
                referenceDepth, crustDensity, mass);
            airyCrust = airyEqualPressures(topography, crustDensity,...
                mantleDensity, referenceDepth, g0, crustGravity);
            flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
                shearModulus, elasticThickness, mantleDensity, ...
                crustDensity, g0, radius);

            error = rmse(flexureCrust, measuredThickness, "all");

            if error < minError
                minError = error;
                minThickness = min(flexureCrust, [], "all");

                % minReferenceDepth = referenceDepth;
                minMantleDensity = mantleDensity;
                minCrustDensity = crustDensity;
                minElasticThickness = elasticThickness;
            end
        end
    end
end
% end

disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Minimum crust thickness @ minimum error: " + string(minThickness / 1e3) + " [km]")
% disp("Reference depth @ minimum error: " + string(minReferenceDepth / 1e3) + " [km]")
disp("Mantle density @ minimum error: " + string(minMantleDensity) + " [kg/m3]")
disp("Crust density @ minimum error: " + string(minCrustDensity) + " [kg/m3]")
disp("Elastic thickness @ minimum error: " + string(minElasticThickness / 1e3) + " [km]")

%% Found results:
% Reference Depth = 24 km
% Mantle Density = 3150 kg/m3
% Crust Density = 2620 kg/m3
% Elastic Thickness = 16.6 km
% Minimum Error = 6.1337 km