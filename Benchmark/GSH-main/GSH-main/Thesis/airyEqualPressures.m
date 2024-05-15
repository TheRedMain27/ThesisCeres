function AiryCrust = airyEqualPressures(topography, crustDensity, ...
    mantleDensity, compensationDepth, g0, crustGravity)
    extraDepth = topography * crustDensity / ...
        (mantleDensity - crustDensity) *  (g0 / crustGravity);
    AiryCrust = compensationDepth + topography + extraDepth;
end