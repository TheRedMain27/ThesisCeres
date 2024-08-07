clear;
close all;
clc;
clearvars;

prattCrust = matfile("InvertedPrattModel/prattCrust.mat").prattCrust;
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
cbar = colorbar;
cbar.Label.String = "Crust Density [kg/m^3]";
xticks(longitudeTicks);
xticklabels(longitudeTickLabels);
yticks(latitudeTicks);
yticklabels(latitudeTickLabels);
hold on
scatter(faculae(:,1), faculae(:,2), "black", "filled")
yline([30, 150], "LineStyle","--","LineWidth",2)
hold off
savefig("Images/FaculaeMap")
saveas(gcf, "Images/PNG/FaculaeMap.png")

sizeFaculae = size(faculae);

faculaeDensities = zeros(sizeFaculae(1),1);

for i = 1:sizeFaculae(1)
    faculaeDensities(i,1) = prattCrust(round(faculae(i,2)), round(faculae(i,1)));
end

% disp(mean(faculaeDensities))
% disp(max(faculaeDensities))
% disp(min(faculaeDensities))

bucketEdges = 1430:10:1740;
sizeBucketEdges = size(bucketEdges);
histogram = zeros(sizeBucketEdges(2) - 1,1);

for i = 1:sizeBucketEdges(2) - 1
    sizeSelection = ...
        size(faculaeDensities((faculaeDensities > bucketEdges(i)) & ...
        (faculaeDensities < bucketEdges(i + 1))));
    histogram(i) = sizeSelection(1);
end

mu = mean(faculaeDensities);
sigma = std(faculaeDensities);
x = 0:31;
y = 10 * sizeFaculae(1) * exp(- 0.5 * ((bucketEdges - mu) / sigma) .^ 2) / ...
    (sigma * sqrt(2 * pi));

figure()
bar(histogram)
xticks(0.5:5:31.5);
xticklabels(string(1430:50:1740));
xlabel("Crust Density [kg/m^3]")
ylabel("Number of Faculae [-]")
hold on
plot(x, y, "LineWidth",2,"Color","red")
hold off
savefig("Images/FaculaeHistogram")
saveas(gcf, "Images/PNG/FaculaeHistogram.png")