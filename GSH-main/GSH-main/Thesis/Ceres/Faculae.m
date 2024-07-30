clear;
close all;
clc;
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pratt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
savefig("Images/PrattFaculaeMap")
saveas(gcf, "Images/PNG/PrattFaculaeMap.png")

sizeFaculae = size(faculae);

faculaeDensities = zeros(sizeFaculae(1),1);

for i = 1:sizeFaculae(1)
    faculaeDensities(i,1) = prattCrust(round(faculae(i,2)), round(faculae(i,1)));
end

disp(mean(faculaeDensities))
disp(max(faculaeDensities))
disp(min(faculaeDensities))

bucketEdges = 1430:10:1740;
sizeBucketEdges = size(bucketEdges);
histogram = zeros(sizeBucketEdges(2) - 1,1);

for i = 1:sizeBucketEdges(2) - 1
    sizeSelection = ...
        size(faculaeDensities((faculaeDensities > bucketEdges(i)) & ...
        (faculaeDensities < bucketEdges(i + 1))));
    histogram(i) = sizeSelection(1);
end

x = 0:31;
mu = mean(faculaeDensities);
sigma = std(faculaeDensities);
y1 = 10 * sizeFaculae(1) * exp(- 0.5 * ((bucketEdges - mu) / sigma) .^ 2) / ...
    (sigma * sqrt(2 * pi));

mu = mean(prattCrust, "all");
sigma = std(prattCrust, 0, "all");
y2 = 10 * sizeFaculae(1) * exp(- 0.5 * ((bucketEdges - mu) / sigma) .^ 2) / ...
    (sigma * sqrt(2 * pi));


figure()
bar(histogram)
xticks(0.5:5:31.5);
xticklabels(string(1430:50:1740));
xlabel("Crust Density [kg/m^3]")
ylabel("Number of Faculae [-]")
hold on
plot(x, y1, "LineWidth",2,"Color","red")
plot(x, y2, "LineWidth",2,"Color","green")
hold off
savefig("Images/PrattFaculaeHistogram")
saveas(gcf, "Images/PNG/PrattFaculaeHistogram.png")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Airy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crustDensity = matfile("InvertedAiryModel/crustDensity.mat").crustDensity;
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
imagesc(crustDensity)
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
savefig("Images/AiryFaculaeMap")
saveas(gcf, "Images/PNG/AiryFaculaeMap.png")

sizeFaculae = size(faculae);

faculaeDensities = zeros(sizeFaculae(1),1);

for i = 1:sizeFaculae(1)
    faculaeDensities(i,1) = crustDensity(round(faculae(i,2)), round(faculae(i,1)));
end

disp(mean(faculaeDensities))
disp(max(faculaeDensities))
disp(min(faculaeDensities))

bucketEdges = 1190:10:1510;
sizeBucketEdges = size(bucketEdges);
histogram = zeros(sizeBucketEdges(2) - 1,1);

for i = 1:sizeBucketEdges(2) - 1
    sizeSelection = ...
        size(faculaeDensities((faculaeDensities > bucketEdges(i)) & ...
        (faculaeDensities < bucketEdges(i + 1))));
    histogram(i) = sizeSelection(1);
end

x = 0:32;
mu = mean(faculaeDensities);
sigma = std(faculaeDensities);
y1 = 10 * sizeFaculae(1) * exp(- 0.5 * ((bucketEdges - mu) / sigma) .^ 2) / ...
    (sigma * sqrt(2 * pi));

mu = mean(crustDensity, "all");
sigma = std(crustDensity, 0, "all");
y2 = 10 * sizeFaculae(1) * exp(- 0.5 * ((bucketEdges - mu) / sigma) .^ 2) / ...
    (sigma * sqrt(2 * pi));

figure()
bar(histogram)
xticks(0.5:5:32.5);
xticklabels(string(1190:50:1510));
xlabel("Crust Density [kg/m^3]")
ylabel("Number of Faculae [-]")
hold on
plot(x, y1, "LineWidth",2,"Color","red")
plot(x, y2, "LineWidth",2,"Color","green")
hold off
savefig("Images/AiryFaculaeHistogram")
saveas(gcf, "Images/PNG/AiryFaculaeHistogram.png")