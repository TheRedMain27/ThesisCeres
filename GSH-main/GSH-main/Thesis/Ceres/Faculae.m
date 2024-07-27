clear;
close all;
clc;
clearvars;

prattCrust = matfile("PrattModel/prattCrust.mat").prattCrust;
faculae = readtable("..\..\..\..\Data\1-s2.0-S0019103517303627-mmc2.csv");
faculae = table2array(faculae(:,1:2));
faculae(faculae > 180) = faculae(faculae > 180);
faculae(:,2) = -faculae(:,2) + 90;

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
imagesc(prattCrust)
axis image
bar = colorbar;
bar.Label.String = "Crust Density [kg/m^3]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
hold on
scatter(faculae(:,1), faculae(:,2), "black", "filled")
hold off

sizeFaculae = size(faculae);

faculaeDensities = zeros(sizeFaculae(1),1);

for i = 1:sizeFaculae(1)
    faculaeDensities(i,1) = prattCrust(round(faculae(i,2)), round(faculae(i,1)));
end

disp(mean(faculaeDensities))