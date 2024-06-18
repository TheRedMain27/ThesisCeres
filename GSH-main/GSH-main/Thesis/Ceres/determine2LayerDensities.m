function [mantleDensity, crustDensity] = ...
    determine2LayerDensities(referenceDepth, radius, mass, momentOfInertia)
    % solves 2 layer densities constrained by mass and moment of inertia
    % innerRadius in meters
    innerRadius = radius - referenceDepth;
    A = [innerRadius ^ 3, radius ^ 3 - innerRadius ^ 3;
        innerRadius ^ 5, radius ^ 5 - innerRadius ^ 5];
    b = [mass / (4/3 * pi); momentOfInertia / (8/15 * pi)];
    densities = linsolve(A,b);
    mantleDensity = densities(1);
    crustDensity = densities(2);
end

