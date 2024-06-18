clear;
close all;
clc;
clearvars;

addpath('..\..\Tools');
addpath('..\SharedTools');

topographyFilename = ...
    "..\..\..\..\Data\Ceres_Dawn_FC_HAMO_DTM_DLR_Global_60ppd_Oct2016.tif";
gravityFilename = "..\..\..\..\Data\Ceres18C_KonoplivEtAl2018.txt";

% % Run if topography has to be reread from the original
% getTopo(topographyFilename);

topography = matfile("topography.mat").topography;

G = 6.6743e-11;
g0 = 0.284; % assuming homogeneous density
radius = 470e3; % (park, 2016)
gravitationalParameter = 62.62736e9; % (konopliv, 2018)
mass = gravitationalParameter / G;
momentOfInertia = 0.375 * mass * radius ^ 2; % (mao and mckinnon, 2018)

directoryName = "AiryTrialOptimization";
if ~exist(directoryName, "dir")
    mkdir(directoryName)
end

SHbounds = [1 18];
save(directoryName + "/SHbounds.mat","SHbounds")
gravity = readGravity(gravityFilename, SHbounds, gravitationalParameter,...
    radius);

mantleDensity = 2367; % (King, 2018)
save(directoryName + "/mantleDensity.mat","mantleDensity")

runOptimization = true;
runVisualization = true;
runFinalVisualization = false;

if runOptimization
    crustDensities = [910, 1000:100:2300];
    sizeCrustDensities = size(crustDensities);
    referenceDepths = 10e3:2.5e3:50e3;
    sizeReferenceDepths = size(referenceDepths);
    
    save(directoryName + "/crustDensities.mat","crustDensities")
    save(directoryName + "/referenceDepths.mat","referenceDepths")
    
    errors = zeros(sizeReferenceDepths(2), sizeCrustDensities(2));
    onlyPositiveThicknesses = ...
        zeros(sizeReferenceDepths(2), sizeCrustDensities(2));
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km
    
    weights = repmat(transpose(sind(0.5:179.5)), 1, 360);
    
    for densityIndex = 1:sizeCrustDensities(2)
        crustDensity = crustDensities(densityIndex);
        for depthIndex = 1:sizeReferenceDepths(2)
            referenceDepth = referenceDepths(depthIndex);
    
            crustGravity = ...
                calculateLowerCrustGravity(radius, referenceDepth, crustDensity, mass);
            airyCrust = airyEqualPressures(topography, crustDensity, ...
                mantleDensity, referenceDepth, g0, crustGravity);
            
            if min(airyCrust, [], "all") > 0
                onlyPositiveThicknesses(depthIndex, densityIndex) = 1;
                
                %% Create model (from GSH package)
                Model = struct();
        
                Model.number_of_layers = 2;
        
                % Additional variables
                Model.GM = gravitationalParameter;
                Model.Re_analyse = radius;
                Model.Re = radius;
                Model.geoid = 'none';
                Model.nmax = 18;
                Model.correct_depth = 0;
        
        
                % % Topo layer
                Model.l1.bound = topography;
                Model.l1.dens  = crustDensity;
        
                % % Mantle
                Model.l2.bound = topography - airyCrust;
                Model.l2.dens = mantleDensity;
        
                % % Bottom
                Model.l3.bound = repmat(-radius, 180, 360);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                SH_coefficients = model_SH_analysis(Model);
                data = model_SH_synthesis(lonLim, latLim, height, ...
                    SHbounds, SH_coefficients, Model);
        
                errors(depthIndex, densityIndex) = ...
                    weightedRMSE(flip(data.vec.R, 1), gravity, weights);
            else
                errors(depthIndex, densityIndex) = NaN;
            end
        end
    end
    save(directoryName + "/errors.mat", "errors");
    save(directoryName + "\onlyPositiveThicknesses.mat", ...
        "onlyPositiveThicknesses")
end

if runVisualization
    onlyPositiveThicknesses = matfile(directoryName + ...
        "/onlyPositiveThicknesses.mat").onlyPositiveThicknesses;
    onlyPositiveThicknesses = logical(onlyPositiveThicknesses);
    errors = matfile(directoryName + "/errors.mat").errors;
    crustDensities = ...
        matfile(directoryName + "/crustDensities.mat").crustDensities;
    referenceDepths = ...
        matfile(directoryName + "/referenceDepths.mat").referenceDepths;
    
    minError = min(errors(onlyPositiveThicknesses), [], "all");
    % minError = min(errors, [], "all");
    [miny, minx] = find(errors == minError);
    referenceDepth = referenceDepths(miny);
    crustDensity = crustDensities(minx);

    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Reference depth @ minimum RMSE: " + ...
        string(referenceDepth / 1e3) + " [km]")
    disp("Crust density @ minimum RMSE: " + string(crustDensity) + " [kg/m3]")
    disp(" ")

    sizeCrustDensities = size(crustDensities);
    densityTickStep = floorDiv(sizeCrustDensities(2), 10);
    if sizeCrustDensities(2) < 10
        densityTickStep = 1;
    end
    sizeReferenceDepths = size(referenceDepths);
    depthTickStep = floorDiv(sizeReferenceDepths(2), 10);
    if sizeReferenceDepths(2) < 10
        depthTickStep = 1;
    end
    
    figure(1)
    colormap("turbo")
    imagesc(errors * 1e5, "AlphaData", onlyPositiveThicknesses)
    % imagesc(errors * 1e5)
    bar = colorbar;
    bar.Label.String = "RMSE [mGal]";
    set(gca, "ColorScale", "log")
    xlabel("Crust Density [kg/m^3]")
    ylabel("Reference Depth [km]")
    yticks(1:depthTickStep:sizeReferenceDepths(2));
    xticks(1:densityTickStep:sizeCrustDensities(2));
    yticklabels(string(referenceDepths(1:depthTickStep:sizeReferenceDepths(2)) / 1e3));
    xticklabels(string(crustDensities(1:densityTickStep:sizeCrustDensities(2))));
    savefig(directoryName + "/errors")
    saveas(gcf, directoryName + "/errors.png")
    
    figure(2)
    colormap("turbo")
    imagesc((errors - minError) / minError, ...
        "AlphaData", onlyPositiveThicknesses)
    % imagesc((errors - minError) / minError)
    bar = colorbar;
    bar.Label.String = "\Delta RMSE / min(RMSE)";
    set(gca, "ColorScale", "log")
    xlabel("Crust Density [kg/m^3]")
    ylabel("Reference Depth [km]")
    yticks(1:depthTickStep:sizeReferenceDepths(2));
    xticks(1:densityTickStep:sizeCrustDensities(2));
    yticklabels(string(referenceDepths(1:depthTickStep:sizeReferenceDepths(2)) / 1e3));
    xticklabels(string(crustDensities(1:densityTickStep:sizeCrustDensities(2))));
    savefig(directoryName + "/errorsNormalized")
    saveas(gcf, directoryName + "/errorsNormalized.png")
end

if runFinalVisualization
    onlyPositiveThicknesses = matfile("AiryErrorsHighLevel/" + ...
        "onlyPositiveThicknesses.mat").onlyPositiveThicknesses;
    onlyPositiveThicknesses = logical(onlyPositiveThicknesses);
    errors = matfile("AiryErrorsHighLevel/errors.mat").errors;
    crustDensities = ...
        matfile("AiryErrorsHighLevel/crustDensities.mat").crustDensities;
    referenceDepths = ...
        matfile("AiryErrorsHighLevel/referenceDepths.mat").referenceDepths;
    
    minError = min(errors(onlyPositiveThicknesses), [], "all");
    [miny, minx] = find(errors == minError);
    referenceDepth = referenceDepths(miny);
    crustDensity = crustDensities(minx);

    disp("HIGH LEVEL")
    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Reference depth @ minimum RMSE: " + ...
        string(referenceDepth / 1e3) + " [km]")
    disp("Crust density @ minimum RMSE: " + string(crustDensity) + " [kg/m3]")
    disp(" ")
    
    figure(1)
    colormap("turbo")
    imagesc(errors * 1e5, "AlphaData", onlyPositiveThicknesses)
    bar = colorbar;
    bar.Label.String = "RMSE [mGal]";
    % set(gca, "ColorScale", "log")
    xlabel("Crust Density [kg/m^3]")
    ylabel("Reference Depth [km]")
    yticks(1:5:50);
    xticks(1:5:25);
    yticklabels(string(referenceDepths(1:5:50) / 1e3));
    xticklabels(string(crustDensities(1:5:25)));
    savefig("Images/AiryErrorsHighLevel")
    saveas(gcf, "Images/PNG/AiryErrorsHighLevel.png")
    
    figure(2)
    colormap("turbo")
    imagesc((errors - minError) / minError, ...
        "AlphaData", onlyPositiveThicknesses)
    bar = colorbar;
    bar.Label.String = "\Delta RMSE / min(RMSE)";
    set(gca, "ColorScale", "log")
    xlabel("Crust Density [kg/m^3]")
    ylabel("Reference Depth [km]")
    yticks(1:5:50);
    xticks(1:5:25);
    yticklabels(string(referenceDepths(1:5:50) / 1e3));
    xticklabels(string(crustDensities(1:5:25)));
    savefig("Images/AiryErrorsNormalizedHighLevel")
    saveas(gcf, "Images/PNG/AiryErrorsNormalizedHighLevel.png")
    
    onlyPositiveThicknesses = matfile("AiryErrorsLowLevel/" + ...
        "onlyPositiveThicknesses.mat").onlyPositiveThicknesses;
    onlyPositiveThicknesses = logical(onlyPositiveThicknesses);
    errors = matfile("AiryErrorsLowLevel/errors.mat").errors;
    crustDensities = ...
        matfile("AiryErrorsLowLevel/crustDensities.mat").crustDensities;
    referenceDepths = ...
        matfile("AiryErrorsLowLevel/referenceDepths.mat").referenceDepths;
    
    minError = min(errors(onlyPositiveThicknesses), [], "all");
    [miny, minx] = find(errors == minError);
    referenceDepth = referenceDepths(miny);
    crustDensity = crustDensities(minx);

    disp("LOW LEVEL")
    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Reference depth @ minimum RMSE: " + ...
        string(referenceDepth / 1e3) + " [km]")
    disp("Crust density @ minimum RMSE: " + string(crustDensity) + " [kg/m3]")
    
    figure(3)
    colormap("turbo")
    imagesc(errors * 1e5, "AlphaData", onlyPositiveThicknesses)
    bar = colorbar;
    bar.Label.String = "RMSE [mGal]";
    % set(gca, "ColorScale", "log")
    xlabel("Crust Density [kg/m^3]")
    ylabel("Reference Depth [km]")
    yticks(1:10:100);
    xticks(1:4:20);
    yticklabels(string(referenceDepths(1:10:100) / 1e3));
    xticklabels(string(crustDensities(1:4:20)));
    savefig("Images/AiryErrorsLowLevel")
    saveas(gcf, "Images/PNG/AiryErrorsLowLevel.png")
    
    figure(4)
    colormap("turbo")
    imagesc((errors - minError) / minError, ...
        "AlphaData", onlyPositiveThicknesses)
    bar = colorbar;
    bar.Label.String = "\Delta RMSE / min(RMSE)";
    set(gca, "ColorScale", "log")
    xlabel("Crust Density [kg/m^3]")
    ylabel("Reference Depth [km]")
    yticks(1:10:100);
    xticks(1:4:20);
    yticklabels(string(referenceDepths(1:10:100) / 1e3));
    xticklabels(string(crustDensities(1:4:20)));
    savefig("Images/AiryErrorsNormalizedLowLevel")
    saveas(gcf, "Images/PNG/AiryErrorsNormalizedLowLevel.png")
end

%% Optimal Values (from low level):
% @ mantle density 2367 [kg/m3]
% @ no negative crust thicknesses
% reference depth = 11.7 [km]
% crust density = 910 [kg/m3]
% rmse = 13.2187 [mGal]