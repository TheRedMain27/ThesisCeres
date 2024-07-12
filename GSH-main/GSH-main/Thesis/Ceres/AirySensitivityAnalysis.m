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

SHbounds = [5 14];
[gravity, gravitySHCoefficients] = readGravity(gravityFilename, ...
    SHbounds, gravitationalParameter, radius);

[n,gravityDegreeVariance] = degreeVariance(gravitySHCoefficients);
gravityDegreeVariance(1:2) = [0, 0];

baseReferenceDepth = 39e3;
baseCrustDensity = 1310;
baseMantleDensity = 2410;

referenceDepths = 10e3:5e3:70e3;
crustDensities = 900:100:2200;
mantleDensities = 1500:100:4000;

% figure(1)
% yline(1)
% xticks(1:(SHbounds(2) - SHbounds(1)));
% xticklabels(string(SHbounds(1):SHbounds(2)));
% xlabel("Spherical Harmonics Degree [-]")
% ylabel("Normalized Degree Variance [-]")
% title("Reference Depth")
% hold on
% 
% for referenceDepth = referenceDepths
%     crustGravity = calculateLowerCrustGravity(radius, referenceDepth, ...
%         baseCrustDensity, mass);
%     airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
%         baseMantleDensity, referenceDepth, g0, crustGravity);
% 
%     %% Create model (from GSH package)
%     Model = struct();
% 
%     Model.number_of_layers = 2;
% 
%     % Additional variables
%     Model.GM = gravitationalParameter;
%     Model.Re_analyse = radius;
%     Model.Re = radius;
%     Model.geoid = 'none';
%     Model.nmax = 18;
%     Model.correct_depth = 0;
% 
% 
%     % % Topo layer
%     Model.l1.bound = topography;
%     Model.l1.dens  = baseCrustDensity;
% 
%     % % Mantle
%     Model.l2.bound = topography - airyCrust;
%     Model.l2.dens = baseMantleDensity;
% 
%     % % Bottom
%     Model.l3.bound = repmat(-radius, 180, 360);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     latLim =    [-89.5 89.5 1];
%     lonLim =    [0.5 359.5 1];
%     height =    0; % LAMO is 375km
% 
%     ModelSHCoefficients = segment_2layer_model(topography, ...
%                     Model.l2.bound, -radius, baseCrustDensity, ...
%                     baseMantleDensity, 3e3, Model);
% 
%     [n,modelDegreeVariance] = degreeVariance(ModelSHCoefficients);
%     modelDegreeVariance(1:2) = [0, 0];
% 
%     normalizedDegreeVariance = modelDegreeVariance ./ gravityDegreeVariance;
%     if referenceDepth == referenceDepths(1)
%         color = "green";
%     elseif referenceDepth == 40e3
%         color = "red";
%     else
%         color = "blue";
%     end
%     plot(normalizedDegreeVariance((SHbounds(1)+1):(SHbounds(2)+1)), "Color", color)
% end
% hold off
% savefig("Images/ReferenceDepthSensitivity")
% saveas(gcf, "Images/PNG/ReferenceDepthSensitivity.png")

figure(2)
yline(1)
xticks(1:(SHbounds(2) - SHbounds(1)));
xticklabels(string(SHbounds(1):SHbounds(2)));
xlabel("Spherical Harmonics Degree [-]")
ylabel("Normalized Degree Variance [-]")
title("Crust Density")
hold on

for crustDensity = crustDensities
    crustGravity = calculateLowerCrustGravity(radius, baseReferenceDepth, ...
        crustDensity, mass);
    airyCrust = airyEqualPressures(topography, crustDensity, ...
        baseMantleDensity, baseReferenceDepth, g0, crustGravity);

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
    Model.l2.dens = baseMantleDensity;

    % % Bottom
    Model.l3.bound = repmat(-radius, 180, 360);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km

    ModelSHCoefficients = segment_2layer_model(topography, ...
                    Model.l2.bound, -radius, crustDensity, ...
                    baseMantleDensity, 3e3, Model);

    [n,modelDegreeVariance] = degreeVariance(ModelSHCoefficients);
    modelDegreeVariance(1:2) = [0, 0];

    normalizedDegreeVariance = modelDegreeVariance ./ gravityDegreeVariance;
    if crustDensity == crustDensities(1)
        color = "green";
    elseif crustDensity == 1300
        color = "red";
    else
        color = "blue";
    end
    plot(normalizedDegreeVariance((SHbounds(1)+1):(SHbounds(2)+1)), "Color", color)
end
hold off
savefig("Images/CrustDensitySensitivity")
saveas(gcf, "Images/PNG/CrustDensitySensitivity.png")

figure(3)
yline(1)
xticks(1:(SHbounds(2) - SHbounds(1)));
xticklabels(string(SHbounds(1):SHbounds(2)));
xlabel("Spherical Harmonics Degree [-]")
ylabel("Normalized Degree Variance [-]")
title("Mantle Density")
hold on

for mantleDensity = mantleDensities
    crustGravity = calculateLowerCrustGravity(radius, baseReferenceDepth, ...
        baseCrustDensity, mass);
    airyCrust = airyEqualPressures(topography, baseCrustDensity, ...
        mantleDensity, baseReferenceDepth, g0, crustGravity);

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
    Model.l1.dens  = baseCrustDensity;

    % % Mantle
    Model.l2.bound = topography - airyCrust;
    Model.l2.dens = mantleDensity;

    % % Bottom
    Model.l3.bound = repmat(-radius, 180, 360);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    height =    0; % LAMO is 375km

    ModelSHCoefficients = segment_2layer_model(topography, ...
                    Model.l2.bound, -radius, baseCrustDensity, ...
                    mantleDensity, 3e3, Model);

    [n,modelDegreeVariance] = degreeVariance(ModelSHCoefficients);
    modelDegreeVariance(1:2) = [0, 0];

    normalizedDegreeVariance = modelDegreeVariance ./ gravityDegreeVariance;
    if mantleDensity == mantleDensities(1)
        color = "green";
    elseif mantleDensity == 2400
        color = "red";
    else
        color = "blue";
    end
    plot(normalizedDegreeVariance((SHbounds(1)+1):(SHbounds(2)+1)), "Color", color)
end
hold off
savefig("Images/MantleDensitySensitivity")
saveas(gcf, "Images/PNG/MantleDensitySensitivity.png")