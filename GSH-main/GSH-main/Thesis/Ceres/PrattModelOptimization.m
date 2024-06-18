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

directoryName = "PrattTrialOptimization";
if ~exist(directoryName, "dir")
    mkdir(directoryName)
end

SHbounds = [1 18];
save(directoryName + "/SHbounds.mat","SHbounds")
gravity = readGravity(gravityFilename, SHbounds, gravitationalParameter,...
    radius);

runOptimization = true;
runVisualization = true;
runFinalVisualization = false;

if runOptimization
    referenceDensities = 2000:100:3000;
    sizeReferenceDensities = size(referenceDensities);
    compensationDepths = 4e3:4e3:40e3;
    sizeCompensationDepths = size(compensationDepths);

    save(directoryName + "/referenceDensities.mat","referenceDensities")
    save(directoryName + "/compensationDepths.mat","compensationDepths")
    
    errors = zeros(sizeCompensationDepths(2), sizeReferenceDensities(2));
    heavierThanIce = ...
        zeros(sizeCompensationDepths(2), sizeReferenceDensities(2));
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km
    
    weights = repmat(transpose(sind(0.5:179.5)), 1, 360);
    
    for densityIndex = 1:sizeReferenceDensities(2)
        referenceDensity = referenceDensities(densityIndex);
        for depthIndex = 1:sizeCompensationDepths(2)
            compensationDepth = compensationDepths(depthIndex);

            prattCrust = pratt(topography, referenceDensity, ...
                compensationDepth);
            
            if min(prattCrust, [], "all") > 910
                heavierThanIce(depthIndex, densityIndex) = 1;
                
                %% Create model (from GSH package)
                Model = struct();
        
                Model.number_of_layers = 1;
        
                % Additional variables
                Model.GM = gravitationalParameter;
                Model.Re_analyse = radius;
                Model.Re = radius;
                Model.geoid = 'none';
                Model.nmax = 18;
                Model.correct_depth = 0;
        
        
                % % Topo layer
                Model.l1.bound = topography;
                Model.l1.dens  = prattCrust;
                
                % % Mantle
                Model.l2.bound = repmat(-compensationDepth, 180, 360);
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
    save(directoryName + "\heavierThanIce.mat", ...
        "heavierThanIce")
end

if runVisualization
    heavierThanIce = matfile(directoryName + ...
        "/heavierThanIce.mat").heavierThanIce;
    heavierThanIce = logical(heavierThanIce);
    errors = matfile(directoryName + "/errors.mat").errors;
    referenceDensities = matfile(directoryName + ...
        "/referenceDensities.mat").referenceDensities;
    compensationDepths = matfile(directoryName + ...
        "/compensationDepths.mat").compensationDepths;
    
    minError = min(errors(heavierThanIce), [], "all");
    % minError = min(errors, [], "all");
    [miny, minx] = find(errors == minError);
    compensationDepth = compensationDepths(miny);
    referenceDensity = referenceDensities(minx);

    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Compensation depth @ minimum RMSE: " + ...
        string(compensationDepth / 1e3) + " [km]")
    disp("Reference density @ minimum RMSE: " + string(referenceDensity) + " [kg/m3]")
    disp(" ")

    sizeReferenceDensities = size(referenceDensities);
    densityTickStep = floorDiv(sizeReferenceDensities(2), 10);
    if sizeReferenceDensities(2) < 10
        densityTickStep = 1;
    end
    sizeCompensationDepths = size(compensationDepths);
    depthTickStep = floorDiv(sizeCompensationDepths(2), 10);
    if sizeCompensationDepths(2) < 10
        depthTickStep = 1;
    end
    
    figure(1)
    colormap("turbo")
    imagesc(errors * 1e5, "AlphaData", heavierThanIce)
    bar = colorbar;
    bar.Label.String = "RMSE [mGal]";
    set(gca, "ColorScale", "log")
    xlabel("Reference Density [kg/m^3]")
    ylabel("Compensation Depth [km]")
    yticks(1:depthTickStep:sizeCompensationDepths(2));
    xticks(1:densityTickStep:sizeReferenceDensities(2));
    yticklabels(string(compensationDepths(...
        1:depthTickStep:sizeCompensationDepths(2)) / 1e3));
    xticklabels(string(referenceDensities(...
        1:densityTickStep:sizeReferenceDensities(2))));
    savefig(directoryName + "/errors")
    saveas(gcf, directoryName + "/errors.png")
    
    figure(2)
    colormap("turbo")
    imagesc((errors - minError) / minError, ...
        "AlphaData", heavierThanIce)
    bar = colorbar;
    bar.Label.String = "\Delta RMSE / min(RMSE)";
    set(gca, "ColorScale", "log")
    xlabel("Reference Density [kg/m^3]")
    ylabel("Compensation Depth [km]")
    yticks(1:depthTickStep:sizeCompensationDepths(2));
    xticks(1:densityTickStep:sizeReferenceDensities(2));
    yticklabels(string(compensationDepths(...
        1:depthTickStep:sizeCompensationDepths(2)) / 1e3));
    xticklabels(string(referenceDensities(...
        1:densityTickStep:sizeReferenceDensities(2))));
    savefig(directoryName + "/errorsNormalized")
    saveas(gcf, directoryName + "/errorsNormalized.png")
end

if runFinalVisualization
    heavierThanIce = matfile("PrattErrorsHighLevel/" + ...
        "heavierThanIce.mat").heavierThanIce;
    heavierThanIce = logical(heavierThanIce);
    errors = matfile("PrattErrorsHighLevel/errors.mat").errors;
    referenceDensities = ...
        matfile("PrattErrorsHighLevel/referenceDensities.mat").referenceDensities;
    compensationDepths = ...
        matfile("PrattErrorsHighLevel/compensationDepths.mat").compensationDepths;
    
    minError = min(errors(heavierThanIce), [], "all");
    [miny, minx] = find(errors == minError);
    compensationDepth = compensationDepths(miny);
    referenceDensity = referenceDensities(minx);

    disp("HIGH LEVEL")
    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Compensation depth @ minimum RMSE: " + ...
        string(compensationDepth / 1e3) + " [km]")
    disp("Reference density @ minimum RMSE: " + string(referenceDensity) + " [kg/m3]")
    disp(" ")
    
    figure(1)
    colormap("turbo")
    imagesc(errors * 1e5, "AlphaData", heavierThanIce)
    bar = colorbar;
    bar.Label.String = "RMSE [mGal]";
    % set(gca, "ColorScale", "log")
    xlabel("Reference Density [kg/m^3]")
    ylabel("Compensation Depth [km]")
    yticklabels(string(compensationDepths / 1e3));
    xticklabels(string(referenceDensities));
    savefig("Images/PrattErrorsHighLevel")
    saveas(gcf, "Images/PNG/PrattErrorsHighLevel.png")
    
    figure(2)
    colormap("turbo")
    imagesc((errors - minError) / minError, ...
        "AlphaData", heavierThanIce)
    bar = colorbar;
    bar.Label.String = "\Delta RMSE / min(RMSE)";
    set(gca, "ColorScale", "log")
    xlabel("Reference Density [kg/m^3]")
    ylabel("Compensation Depth [km]")
    yticklabels(string(compensationDepths / 1e3));
    xticklabels(string(referenceDensities));
    savefig("Images/PrattErrorsNormalizedHighLevel")
    saveas(gcf, "Images/PNG/PrattErrorsNormalizedHighLevel.png")
    
    heavierThanIce = matfile("PrattErrorsLowLevel/" + ...
        "heavierThanIce.mat").heavierThanIce;
    heavierThanIce = logical(heavierThanIce);
    errors = matfile("PrattErrorsLowLevel/errors.mat").errors;
    referenceDensities = ...
        matfile("PrattErrorsLowLevel/referenceDensities.mat").referenceDensities;
    compensationDepths = ...
        matfile("PrattErrorsLowLevel/compensationDepths.mat").compensationDepths;
    
    minError = min(errors(heavierThanIce), [], "all");
    [miny, minx] = find(errors == minError);
    compensationDepth = compensationDepths(miny);
    referenceDensity = referenceDensities(minx);

    disp("LOW LEVEL")
    disp("Minimum RMSE: " + string(minError * 1e5) + " [mGal]")
    disp("Compensation depth @ minimum RMSE: " + ...
        string(compensationDepth / 1e3) + " [km]")
    disp("Reference density @ minimum RMSE: " + string(referenceDensity) + " [kg/m3]")
    
    figure(3)
    colormap("turbo")
    imagesc(errors * 1e5, "AlphaData", heavierThanIce)
    bar = colorbar;
    bar.Label.String = "RMSE [mGal]";
    % set(gca, "ColorScale", "log")
    xlabel("Reference Density [kg/m^3]")
    ylabel("Compensation Depth [km]")
    yticks(1:20:201);
    xticks(1:2:21);
    yticklabels(string(compensationDepths(1:20:201) / 1e3));
    xticklabels(string(referenceDensities(1:2:21)));
    savefig("Images/PrattErrorsLowLevel")
    saveas(gcf, "Images/PNG/PrattErrorsLowLevel.png")
    
    figure(4)
    colormap("turbo")
    imagesc((errors - minError) / minError, ...
        "AlphaData", heavierThanIce)
    bar = colorbar;
    bar.Label.String = "\Delta RMSE / min(RMSE)";
    set(gca, "ColorScale", "log")
    xlabel("Reference Density [kg/m^3]")
    ylabel("Compensation Depth [km]")
    yticks(1:20:201);
    xticks(1:2:21);
    yticklabels(string(compensationDepths(1:20:201) / 1e3));
    xticklabels(string(referenceDensities(1:2:21)));
    savefig("Images/PrattErrorsNormalizedLowLevel")
    saveas(gcf, "Images/PNG/PrattErrorsNormalizedLowLevel.png")
end

%% Optimal Values (from low level):
% @ no density lower than 910 [kg/m3]
% reference density = [kg/m3]
% compensation depth = [km]
% RMSE =  [mGal]