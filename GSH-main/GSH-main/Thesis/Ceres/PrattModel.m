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
SHbounds =  [5 14];
save("PrattModel/SHbounds.mat","SHbounds")

[gravity, gravitySHCoefficients] = readGravity(gravityFilename, ...
    SHbounds, gravitationalParameter, radius);

compensationDepth = 70e3;
referenceDensity = 1570;
save("PrattModel/compensationDepth.mat","compensationDepth")
save("PrattModel/referenceDensity.mat","referenceDensity")

prattCrust = pratt(topography, referenceDensity, compensationDepth);
save("PrattModel/prattCrust.mat","prattCrust")

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];
lonLim =    [0.5 359.5 1];
height =    0; % LAMO is 375km

ModelSHCoefficients = segment_2layer_model(topography, ...
                    Model.l2.bound, -compensationDepth, prattCrust, 0, ...
                    compensationDepth / 5, Model);
data = model_SH_synthesis(lonLim, latLim, height, SHbounds, ...
    ModelSHCoefficients, Model);

gravityResidual = gravity - flip(data.vec.R, 1);
save("PrattModel/gravityResidual.mat","gravityResidual")

[n,modelDegreeVariance] = degreeVariance(ModelSHCoefficients);
modelDegreeVariance(1:2) = [0, 0];
[n,gravityDegreeVariance] = degreeVariance(gravitySHCoefficients);
gravityDegreeVariance(1:2) = [0, 0];

normalizedDegreeVariance = modelDegreeVariance ./ gravityDegreeVariance;

figure(1)
plot(normalizedDegreeVariance((SHbounds(1)+1):(SHbounds(2)+1)))
yline(1)
% set(gca,"Yscale", "log")
xticks(1:(SHbounds(2) - SHbounds(1)));
xticklabels(string(SHbounds(1):SHbounds(2)));
xlabel("Spherical Harmonics Degree [-]")
ylabel("Normalized Degree Variance [-]")

latitudeTicks = 0:30:180;
latitudeTickLabels = string(flip(-90:30:90));
longitudeTicks = 0:60:360;
longitudeTickLabels = string(-180:60:180);

figurePosition = get(groot, 'DefaultFigurePosition');
figurePosition(1) = figurePosition(1) - (2 * figurePosition(4) - ...
    figurePosition(3)) / 2;
figurePosition(3) = 2 * figurePosition(4);

figure('Position', figurePosition)
colormap("turbo")
imagesc(gravity * 1e5)
axis image
bar = colorbar;
bar.Label.String = "Gravitational Acceleration [mGal]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);

figure('Position', figurePosition)
colormap("turbo")
imagesc(prattCrust)
axis image
bar = colorbar;
bar.Label.String = "Crust Density [kg/m^3]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("PrattModel/PrattCrustDensity")
saveas(gcf, "PrattModel/PrattCrustDensity.png")

figure('Position', figurePosition)
colormap("turbo")
imagesc(flip(data.vec.R, 1) * 1e5)
axis image
bar = colorbar;
bar.Label.String = "Model Gravitational Acceleration [mGal]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("PrattModel/PrattGravity")
saveas(gcf, "PrattModel/PrattGravity.png")

figure('Position', figurePosition)
colormap("turbo")
imagesc(gravityResidual * 1e5)
axis image
bar = colorbar;
bar.Label.String = "Gravity Residual [mGal]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
savefig("PrattModel/PrattResidual")
saveas(gcf, "PrattModel/PrattResidual.png")