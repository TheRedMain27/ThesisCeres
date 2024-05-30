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

topography = readTopo(topoFileName);
mohoDepth = readMoho(mohoFileName);
geoid = getGeoid();

measuredThickness = mohoDepth * 1e3 - geoid + topography;

bulkModulus = 128.8e9; % olivine, (Mao et al., 2015)
shearModulus = 81.6e9; % olivine, (Mao er al., 2015)
g0 = 9.80665;
radius = 6371e3;
mass = 3.986004418e14 / 6.6743e-11;

runOptimization = false;
runVisualization = true;

if runOptimization
    % optimal values from split crust
    referenceDepth = 32e3;
    mantleDensity = 2700;
    crustDensity = 1900;
    elasticThickness = 79e3;
    
    depthAlterations = -3e3:1e3:3e3;
    mantleDensityAlterations = -1.5e3:100:2.5e3;
    crustDensityAlterations = -1.2e3:100:1.5e3;
    thicknessAlterations = -60e3:1e3:40e3;
    
    optimalValues = zeros([4, 4]);
    
    tic
    %% Meso- and Cenozoic
    [optimalReferenceDepth, optimalMantleDensity, optimalCrustDensity, ...
        optimalElasticThickness] = successiveOptimization(...
        referenceDepth, mantleDensity, crustDensity, elasticThickness, ...
        bulkModulus, shearModulus, radius, mass, g0, topography, ...
        mesozoicCenozoicMap, measuredThickness, depthAlterations, ...
        mantleDensityAlterations, crustDensityAlterations, ...
        thicknessAlterations);
    optimalValues(:,1) = [optimalReferenceDepth; optimalMantleDensity; ...
        optimalCrustDensity; optimalElasticThickness];
    
    crustGravity = calculateLowerCrustGravity(radius, ...
        optimalReferenceDepth, optimalCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, optimalCrustDensity,...
        optimalMantleDensity, optimalReferenceDepth, g0, crustGravity);
    mesozoicCenozoicFlexureCrust = flexureThinShell(airyCrust, bulkModulus,...
        shearModulus, optimalElasticThickness, optimalMantleDensity, ...
        optimalCrustDensity, g0, radius);
    
    %% Paleozoic
    [optimalReferenceDepth, optimalMantleDensity, optimalCrustDensity, ...
        optimalElasticThickness] = successiveOptimization(...
        referenceDepth, mantleDensity, crustDensity, elasticThickness, ...
        bulkModulus, shearModulus, radius, mass, g0, topography, ...
        paleozoicMap, measuredThickness, depthAlterations, ...
        mantleDensityAlterations, crustDensityAlterations, ...
        thicknessAlterations);
    optimalValues(:,2) = [optimalReferenceDepth; optimalMantleDensity; ...
        optimalCrustDensity; optimalElasticThickness];
    
    crustGravity = calculateLowerCrustGravity(radius, ...
        optimalReferenceDepth, optimalCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, optimalCrustDensity,...
        optimalMantleDensity, optimalReferenceDepth, g0, crustGravity);
    paleozoicFlexureCrust = flexureThinShell(airyCrust, bulkModulus,...
        shearModulus, optimalElasticThickness, optimalMantleDensity, ...
        optimalCrustDensity, g0, radius);
    
    %% Proterozoic
    [optimalReferenceDepth, optimalMantleDensity, optimalCrustDensity, ...
        optimalElasticThickness] = successiveOptimization(...
        referenceDepth, mantleDensity, crustDensity, elasticThickness, ...
        bulkModulus, shearModulus, radius, mass, g0, topography, ...
        proterozoicMap, measuredThickness, depthAlterations, ...
        mantleDensityAlterations, crustDensityAlterations, ...
        thicknessAlterations);
    optimalValues(:,3) = [optimalReferenceDepth; optimalMantleDensity; ...
        optimalCrustDensity; optimalElasticThickness];
    
    crustGravity = calculateLowerCrustGravity(radius, ...
        optimalReferenceDepth, optimalCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, optimalCrustDensity,...
        optimalMantleDensity, optimalReferenceDepth, g0, crustGravity);
    proterozoicFlexureCrust = flexureThinShell(airyCrust, bulkModulus,...
        shearModulus, optimalElasticThickness, optimalMantleDensity, ...
        optimalCrustDensity, g0, radius);
    
    %% Archean
    [optimalReferenceDepth, optimalMantleDensity, optimalCrustDensity, ...
        optimalElasticThickness] = successiveOptimization(...
        referenceDepth, mantleDensity, crustDensity, elasticThickness, ...
        bulkModulus, shearModulus, radius, mass, g0, topography, ...
        archeanMap, measuredThickness, depthAlterations, ...
        mantleDensityAlterations, crustDensityAlterations, ...
        thicknessAlterations);
    optimalValues(:,4) = [optimalReferenceDepth; optimalMantleDensity; ...
        optimalCrustDensity; optimalElasticThickness];
    
    crustGravity = calculateLowerCrustGravity(radius, ...
        optimalReferenceDepth, optimalCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, optimalCrustDensity,...
        optimalMantleDensity, optimalReferenceDepth, g0, crustGravity);
    archeanFlexureCrust = flexureThinShell(airyCrust, bulkModulus,...
        shearModulus, optimalElasticThickness, optimalMantleDensity, ...
        optimalCrustDensity, g0, radius);
    
    multiCrust = mesozoicCenozoicFlexureCrust .* mesozoicCenozoicMap + ...
        paleozoicFlexureCrust .* paleozoicMap + ...
        proterozoicFlexureCrust .* proterozoicMap + ...
        archeanFlexureCrust .* archeanMap;
    disp(" ")
    toc
    disp(" ")
    
    save("optimalValues.mat", "optimalValues")
    save("multiCrust.mat", "multiCrust")
    
    disp("Meso- and cenozoic crust error: " + ...
        string(selectiveRMSE(mesozoicCenozoicFlexureCrust, ...
        measuredThickness, mesozoicCenozoicMap) / 1e3) + " [km]")
    disp("Paleozoic crust error: " + ...
        string(selectiveRMSE(paleozoicFlexureCrust, measuredThickness, ...
        paleozoicMap) / 1e3) + " [km]")
    disp("Proterozoic crust error: " + ...
        string(selectiveRMSE(proterozoicFlexureCrust, measuredThickness, ...
        proterozoicMap) / 1e3) + " [km]")
    disp("Archean crust error: " + ...
        string(selectiveRMSE(archeanFlexureCrust, measuredThickness, ...
        archeanMap) / 1e3) + " [km]")
    disp("Combined crust error: " + string(selectiveRMSE(multiCrust, ...
        measuredThickness, fullCrustMap) / 1e3) + " [km]")
end

if runVisualization
    load("optimalValues.mat");
    maps(:,:,1) = mesozoicCenozoicMap;
    maps(:,:,2) = paleozoicMap;
    maps(:,:,3) = proterozoicMap;
    maps(:,:,4) = archeanMap;
    multiCrust = zeros([180, 360]);
    
    for index = 1:4
        optimalReferenceDepth = optimalValues(1,index);
        optimalMantleDensity = optimalValues(2,index);
        optimalCrustDensity = optimalValues(3,index);
        optimalElasticThickness = optimalValues(4,index);
        
        crustGravity = calculateLowerCrustGravity(radius, ...
            optimalReferenceDepth, optimalCrustDensity, mass);
        airyCrust = airyEqualPressures(topography, optimalCrustDensity,...
            optimalMantleDensity, optimalReferenceDepth, g0, crustGravity);
        flexureCrust = flexureThinShell(airyCrust, bulkModulus,...
            shearModulus, optimalElasticThickness, optimalMantleDensity, ...
            optimalCrustDensity, g0, radius);
        
        multiCrust = multiCrust + flexureCrust .* maps(:,:,index);
    end

    disp("Meso- and cenozoic crust error: " + ...
        string(selectiveRMSE(multiCrust, measuredThickness, ...
        mesozoicCenozoicMap) / 1e3) + " [km]")
    disp("Paleozoic crust error: " + ...
        string(selectiveRMSE(multiCrust, measuredThickness, ...
        paleozoicMap) / 1e3) + " [km]")
    disp("Proterozoic crust error: " + ...
        string(selectiveRMSE(multiCrust, measuredThickness, ...
        proterozoicMap) / 1e3) + " [km]")
    disp("Archean crust error: " + ...
        string(selectiveRMSE(multiCrust, measuredThickness, ...
        archeanMap) / 1e3) + " [km]")
    disp("Combined crust error: " + string(selectiveRMSE(multiCrust, ...
        measuredThickness, fullCrustMap) / 1e3) + " [km]")

    [error, percentageError, wrongRMS, minError, maxError] = ...
        characterizeError(measuredThickness, multiCrust);
    
    multiCrust(~fullCrustMap) = mean(multiCrust, "all");
    error(~fullCrustMap) = 0;
    percentageError(~fullCrustMap) = mean(percentageError, "all");
    
    modelDepth = multiCrust - topography + geoid;
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
    savefig("Images/MultiCrust/MeasuredCrustThickness")
    saveas(gcf, "Images/MultiCrust/PNG/MeasuredCrustThickness.png")
    
    figure(2)
    colormap('hot');
    imagesc(multiCrust / 1e3, "AlphaData", fullCrustMap);
    set(gca, 'color', [0 0 1])
    bar = colorbar;
    bar.Label.String = "Model Crust Thickness [km]";
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
    savefig("Images/MultiCrust/ModelCrustThickness")
    saveas(gcf, "Images/MultiCrust/PNG/ModelCrustThickness.png")

    figure(3)
    colormap('hot');
    imagesc(error / 1e3, "AlphaData", fullCrustMap);
    set(gca, 'color', [0 0 1])
    bar = colorbar;
    bar.Label.String = "Model Error [km]";
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
    savefig("Images/MultiCrust/ModelError")
    saveas(gcf, "Images/MultiCrust/PNG/ModelError.png")
    
    percentageError(percentageError > 100) = 100;
    
    figure(4)
    colormap('hot');
    imagesc(percentageError, "AlphaData", fullCrustMap);
    set(gca, 'color', [0 0 1])
    bar = colorbar;
    bar.Label.String = "Model Error [%]";
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
    savefig("Images/MultiCrust/PercentageModelError")
    saveas(gcf, "Images/MultiCrust/PNG/PercentageModelError.png")
    
    errorUncertaintyRatio(errorUncertaintyRatio > 2) = 2;
    
    figure(5)
    colormap('hot');
    imagesc(errorUncertaintyRatio, "AlphaData", fullCrustMap);
    set(gca, 'color', [0 0 1])
    bar = colorbar;
    bar.Label.String = "Depth Error to Uncertainty Ratio";
    xticks(longitudeTicks);
    xticklabels(longitudeTickLabels);
    yticks(latitudeTicks);
    yticklabels(latitudeTickLabels);
    savefig("Images/MultiCrust/ErrorUncertaintyRatio")
    saveas(gcf, "Images/MultiCrust/PNG/ErrorUncertaintyRatio.png")
end