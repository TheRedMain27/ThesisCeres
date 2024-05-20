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

baseBulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
baseShearModulus = 81.6e9; % olivine, (Mao er al., 2015)
baseElasticThickness = 34e3; % PPI slides

baseMantleDensity = 3320; % olivine
baseCrustDensity = 2760;
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;
baseReferenceDepth = mean(measuredThickness, "all");

%%%%%%%%%%%%%%%%%%%%%%%%%%% Reference Depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%
depthAlterations = -2e3:100:2e3;
sizeDepthAlterations = size(depthAlterations);
errors = zeros(1,sizeDepthAlterations(2));
minimumThicknesses = zeros(1,sizeDepthAlterations(2));
indexes = 1:sizeDepthAlterations(2);

for index = indexes
    referenceDepth = baseReferenceDepth + depthAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        referenceDepth, baseCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
        baseMantleDensity, referenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, baseBulkModulus,...
        baseShearModulus, baseElasticThickness, baseMantleDensity, baseCrustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

minError = min(errors, [], "all");
minReferenceDepth = (baseReferenceDepth + depthAlterations(errors == minError)) / 1e3;
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Reference Depth @ minimum error: " + string(minReferenceDepth) + " [km]")

figure(1)
yyaxis left
plot((baseReferenceDepth + depthAlterations) / 1e3, errors / 1e3)
xlabel("Reference Depth [km]")
ylabel("RMSE [km]")
yyaxis right
plot((baseReferenceDepth + depthAlterations) / 1e3, minimumThicknesses / 1e3)
ylabel("Minimum Crustal Thickness [km]")
savefig("Images/ReferenceDepthOptimization")
saveas(gcf, "Images/PNG/ReferenceDepthOptimization.png")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
densityAlterations = -200:20:200;
sizeDensityAlterations = size(densityAlterations);

%% mantle density
errors = zeros(1,sizeDensityAlterations(2));
minimumThicknesses = zeros(1,sizeDensityAlterations(2));
indexes = 1:sizeDensityAlterations(2);

for index = indexes
    mantleDensity = baseMantleDensity + densityAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        baseReferenceDepth, baseCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
        mantleDensity, baseReferenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, baseBulkModulus,...
        baseShearModulus, baseElasticThickness, mantleDensity, baseCrustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

minError = min(errors, [], "all");
minMantleDensity = (baseMantleDensity + densityAlterations(errors == minError));
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Mantle Density @ minimum error: " + string(minMantleDensity) + " [kg/m3]")

figure(2)
yyaxis left
plot((baseMantleDensity + densityAlterations), errors / 1e3)
xlabel("Mantle Density [kg/m^{3}]")
ylabel("RMSE [km]")
yyaxis right
plot((baseMantleDensity + densityAlterations), minimumThicknesses / 1e3)
ylabel("Minimum Crustal Thickness [km]")
savefig("Images/MantleDensityOptimization")
saveas(gcf, "Images/PNG/MantleDensityOptimization.png")

%% crust density
errors = zeros(1,sizeDensityAlterations(2));
minimumThicknesses = zeros(1,sizeDensityAlterations(2));
indexes = 1:sizeDensityAlterations(2);

for index = indexes
    crustDensity = baseCrustDensity + densityAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        baseReferenceDepth, crustDensity, mass);
    airyCrust = airyEqualPressures(topography, crustDensity, ...
        baseMantleDensity, baseReferenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, baseBulkModulus,...
        baseShearModulus, baseElasticThickness, baseMantleDensity, crustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

minError = min(errors, [], "all");
minCrustDensity = (baseCrustDensity + densityAlterations(errors == minError));
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Crust Density @ minimum error: " + string(minCrustDensity) + " [kg/m3]")

figure(3)
yyaxis left
plot((baseCrustDensity + densityAlterations), errors / 1e3)
xlabel("Crust Density [kg/m^{3}]")
ylabel("RMSE [km]")
yyaxis right
plot((baseCrustDensity + densityAlterations), minimumThicknesses / 1e3)
ylabel("Minimum Crustal Thickness [km]")
savefig("Images/CrustDensityOptimization")
saveas(gcf, "Images/PNG/CrustDensityOptimization.png")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Elastic Thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thicknessAlterations = -24e3:1e3:6e3;
sizeThicknessAlterations = size(thicknessAlterations);
errors = zeros(1,sizeThicknessAlterations(2));
minimumThicknesses = zeros(1,sizeThicknessAlterations(2));
indexes = 1:sizeThicknessAlterations(2);

for index = indexes
    elasticThickness = baseElasticThickness + thicknessAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        baseReferenceDepth, baseCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
        baseMantleDensity, baseReferenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, baseBulkModulus,...
        baseShearModulus, elasticThickness, baseMantleDensity, baseCrustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

minError = min(errors, [], "all");
minElasticThickness = (baseElasticThickness + thicknessAlterations(errors == minError)) / 1e3;
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Elastic Thickness @ minimum error: " + string(minElasticThickness) + " [km]")

figure(4)
yyaxis left
plot((baseElasticThickness + thicknessAlterations) / 1e3, errors / 1e3)
xlabel("Elastic Thickness [km]")
ylabel("RMSE [km]")
yyaxis right
plot((baseElasticThickness + thicknessAlterations) / 1e3, minimumThicknesses / 1e3)
ylabel("Minimum Crustal Thickness [km]")
savefig("Images/ElasticThicknessOptimization")
saveas(gcf, "Images/PNG/ElasticThicknessOptimization.png")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Moduli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modulusAlterations = -50e9:50e8:50e9;
sizeModulusAlterations = size(modulusAlterations);

%% bulk modulus
errors = zeros(1,sizeModulusAlterations(2));
minimumThicknesses = zeros(1,sizeModulusAlterations(2));
indexes = 1:sizeModulusAlterations(2);

for index = indexes
    bulkModulus = baseBulkModulus + modulusAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        baseReferenceDepth, baseCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
        baseMantleDensity, baseReferenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
        baseShearModulus, baseElasticThickness, baseMantleDensity, baseCrustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

minError = min(errors, [], "all");
minBulkModulus = (baseBulkModulus + modulusAlterations(errors == minError)) / 1e9;
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Bulk Modulus @ minimum error: " + string(minBulkModulus) + " [GPa]")

figure(5)
yyaxis left
plot((baseBulkModulus + modulusAlterations) / 1e9, errors / 1e3)
xlabel("Bulk Modulus [GPa]")
ylabel("RMSE [km]")
yyaxis right
plot((baseBulkModulus + modulusAlterations) / 1e9, minimumThicknesses / 1e3)
ylabel("Minimum Crustal Thickness [km]")
savefig("Images/BulkModulusOptimization")
saveas(gcf, "Images/PNG/BulkModulusOptimization.png")

%% crust density
errors = zeros(1,sizeModulusAlterations(2));
minimumThicknesses = zeros(1,sizeModulusAlterations(2));
indexes = 1:sizeModulusAlterations(2);

for index = indexes
    shearModulus = baseShearModulus + modulusAlterations(index);
    crustGravity = calculateLowerCrustGravity(radius, ...
        baseReferenceDepth, baseCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
        baseMantleDensity, baseReferenceDepth, g0, crustGravity);
    flexureCrust = flexureThinShell(airyCrust, baseBulkModulus,...
        shearModulus, baseElasticThickness, baseMantleDensity, baseCrustDensity, ...
        g0, radius);
    errors(index) = rmse(flexureCrust, measuredThickness, "all");
    minimumThicknesses(index) = min(flexureCrust, [], "all");
end

minError = min(errors, [], "all");
minShearModulus = (baseShearModulus + modulusAlterations(errors == minError)) / 1e9;
disp("Minimum error: " +  string(minError / 1e3) + " [km]")
disp("Shear Modulus @ minimum error: " + string(minShearModulus) + " [GPa]")

figure(6)
yyaxis left
plot((baseShearModulus + modulusAlterations) / 1e9, errors / 1e3)
xlabel("Shear Modulus [GPa]")
ylabel("RMSE [km]")
yyaxis right
plot((baseShearModulus + modulusAlterations) / 1e9, minimumThicknesses / 1e3)
ylabel("Minimum Crustal Thickness [km]")
savefig("Images/ShearModulusOptimization")
saveas(gcf, "Images/PNG/ShearModulusOptimization.png")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%