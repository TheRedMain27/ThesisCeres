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

directoryName = "AiryDegreeVarianceTry3";
if ~exist(directoryName, "dir")
    mkdir(directoryName)
end

SHbounds = [5 14];
save(directoryName + "/SHbounds.mat","SHbounds")
[gravity, gravitySHCoefficients] = readGravity(gravityFilename, SHbounds, gravitationalParameter,...
    radius);

[n,gravityDegreeVariance] = degreeVariance(gravitySHCoefficients);

runOptimization = true;
printResult = true;

if runOptimization
    referenceDepths = 30e3:1e3:50e3;
    sizeReferenceDepths = size(referenceDepths);
    
    save(directoryName + "/referenceDepths.mat","referenceDepths")
    
    errors = zeros(sizeReferenceDepths(2), 1);
    onlyPositiveThicknesses = ...
        zeros(sizeReferenceDepths(2), 1);
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km

    Model = struct();
        
    Model.number_of_layers = 2;
    
    % Additional variables
    Model.GM = gravitationalParameter;
    Model.Re_analyse = radius;
    Model.Re = radius;
    Model.geoid = 'none';
    Model.nmax = 18;
    Model.correct_depth = 0;

    Model.l1.bound = topography;
    
    for depthIndex = 1:sizeReferenceDepths(2)
        referenceDepth = referenceDepths(depthIndex);
        [mantleDensity, crustDensity] = ...
            determine2LayerDensities(referenceDepth, radius, mass, ...
            momentOfInertia);

        crustGravity = ...
            calculateLowerCrustGravity(radius, referenceDepth, ...
            crustDensity, mass);
        airyCrust = airyEqualPressures(topography, crustDensity, ...
            mantleDensity, referenceDepth, g0, crustGravity);
        
        if min(airyCrust, [], "all") > 0
            onlyPositiveThicknesses(depthIndex, 1) = 1;

            Model.l1.dens  = crustDensity;
            Model.l2.bound = topography - airyCrust;
            Model.l2.dens = mantleDensity;
            Model.l3.bound = repmat(-radius, 180, 360);
    
            ModelSHCoefficients = segment_2layer_model(topography, ...
                Model.l2.bound, -radius, crustDensity, ...
                mantleDensity, 3e3, Model);
            [n,modelDegreeVariance] = ...
                degreeVariance(ModelSHCoefficients);
            modelDegreeVariance(1:2) = [0, 0];
    
            errors(depthIndex, 1) = ...
                rmse(modelDegreeVariance(SHbounds(1):SHbounds(2)), ...
                gravityDegreeVariance(SHbounds(1):SHbounds(2)), "all");
        else
            errors(depthIndex, 1) = NaN;
        end
    end
    save(directoryName + "/errors.mat", "errors");
    save(directoryName + "\onlyPositiveThicknesses.mat", ...
        "onlyPositiveThicknesses")
end

if printResult
    onlyPositiveThicknesses = matfile(directoryName + ...
        "/onlyPositiveThicknesses.mat").onlyPositiveThicknesses;
    onlyPositiveThicknesses = logical(onlyPositiveThicknesses);
    errors = matfile(directoryName + "/errors.mat").errors;
    referenceDepths = ...
        matfile(directoryName + "/referenceDepths.mat").referenceDepths;
    
    minError = min(errors(onlyPositiveThicknesses), [], "all");
    % minError = min(errors, [], "all");
    index = find(errors == minError);
    referenceDepth = referenceDepths(index);
    [mantleDensity, crustDensity] = ...
            determine2LayerDensities(referenceDepth, radius, mass, ...
            momentOfInertia);

    disp("Minimum RMSE: " + string(minError * 1e10) + " [mGal^2]")
    disp("Reference depth @ minimum RMSE: " + ...
        string(referenceDepth / 1e3) + " [km]")
    disp("Crust density @ minimum RMSE: " + string(crustDensity) + " [kg/m3]")
    disp("Mantle density @ minimum RMSE: " + string(mantleDensity) + " [kg/m3]")
    disp(" ")

    figure(1)
    plot(referenceDepths / 1e3, errors(:,1) * 1e10)
    xlabel("Reference Depth [km]")
    ylabel("RMSE [mGal^2]")
    savefig(directoryName + "/errors")
    saveas(gcf, directoryName + "/errors.png")
end

%% Optimal Values (from low level):
% @ mantle density 2367 [kg/m3]
% @ no negative crust thicknesses
% reference depth =  [km]
% crust density =  [kg/m3]
% rmse =  [mGal]