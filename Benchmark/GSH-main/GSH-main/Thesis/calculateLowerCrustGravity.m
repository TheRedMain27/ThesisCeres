function gravity = calculateLowerCrustGravity(radius, referenceDepth,...
    crustDensity, mass)
    % calculates the gravity at the crust-mantle boundary
    crustMass = (4 / 3) * pi * ...
        (radius ^ 3 - (radius - referenceDepth) ^ 3) * crustDensity;
    residualMass = mass - crustMass;
    gravity = 6.6743e-11 * residualMass / (radius - referenceDepth) ^ 2;
end

