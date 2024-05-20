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
elasticThickness = 34e3; % PPI slides

baseMantleDensity = 3320; % olivine
baseCrustDensity = 2760;
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;
baseReferenceDepth = mean(measuredThickness, "all");

mantleDensity = baseMantleDensity;
crustDensity = baseCrustDensity;
depthAlterations = -2e3:100:2e3;
sizeDepthAlterations = size(depthAlterations);
errors = zeros(1,sizeDepthAlterations(2));
minimumThicknesses = zeros(1,sizeDepthAlterations(2));
indexes = 1:sizeDepthAlterations(2);
for index = indexes
    referenceDepth = baseReferenceDepth + depthAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        referenceDepth, crustDensity, mass);
    airyCrust = airyEqualPressures(topography, crustDensity, ...
        mantleDensity, referenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
        shearModulus, elasticThickness, mantleDensity, crustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

figure(1)
plot((baseReferenceDepth + depthAlterations) / 1e3, errors)
savefig("Images/compVerr")
figure(2)
plot((baseReferenceDepth + depthAlterations) / 1e3, minimumThicknesses)
savefig("Images/compVminthick")

baseMantleDensity = baseMantleDensity + 50;
baseCrustDensity = baseCrustDensity - 50;
referenceDepth = baseReferenceDepth;
densityAlterations = -200:20:200;
sizeDensityAlterations = size(densityAlterations);
errors = zeros(sizeDensityAlterations(2),sizeDensityAlterations(2));
minimumThicknesses = zeros(sizeDensityAlterations(2),sizeDensityAlterations(2));
indexes = 1:sizeDensityAlterations(2);

for mantleIndex = indexes
    mantleDensity = baseMantleDensity + densityAlterations(mantleIndex);
    for crustIndex = indexes
        crustDensity = baseCrustDensity + densityAlterations(crustIndex);
        crustGravity = calculateLowerCrustGravity(radius, ...
            referenceDepth, crustDensity, mass);
        airyCrust = airyEqualPressures(topography, crustDensity, ...
            mantleDensity, referenceDepth, g0, crustGravity);
        flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
            shearModulus, elasticThickness, mantleDensity, crustDensity, ...
            g0, radius);
        errors(mantleIndex, crustIndex) = rmse(flexureCrust, ...
            measuredThickness, "all");
        minimumThicknesses(mantleIndex, crustIndex) = min(flexureCrust, [], "all");
    end
end

minError = min(errors, [], "all");
[miny, minx] = find(errors == minError);
minMantleDensity = baseMantleDensity + densityAlterations(miny);
minCrustDensity = baseCrustDensity + densityAlterations(minx);
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Mantle density @ minimum error: " + string(minMantleDensity) + " [kg/m3]")
disp("Crust density @ minimum error: " + string(minCrustDensity) + " [kg/m3]")
errors(errors > 8e3) = 8e3;

figure(1)
colormap("hot")
imagesc(flip(errors, 1) / 1e3)
xticks(1:2:sizeDensityAlterations(2))
yticks(1:2:sizeDensityAlterations(2))
xticklabels(string(densityAlterations(1:2:sizeDensityAlterations(2)) + baseCrustDensity))
yticklabels(string(flip(densityAlterations(1:2:sizeDensityAlterations(2)) + baseMantleDensity)))
bar = colorbar;
bar.Label.String = "RMSE [km]";
savefig("Images/densityVerr")

figure(2)
colormap("hot")
imagesc(flip(minimumThicknesses, 1) / 1e3)
xticks(1:2:sizeDensityAlterations(2))
yticks(1:2:sizeDensityAlterations(2))
xticklabels(string(densityAlterations(1:2:sizeDensityAlterations(2)) + baseCrustDensity))
yticklabels(string(flip(densityAlterations(1:2:sizeDensityAlterations(2)) + baseMantleDensity)))
bar = colorbar;
bar.Label.String = "Minimum Crust Thickness [km]";
savefig("Images/densityVminthick")