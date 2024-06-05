function PrattCrustDensity = pratt(topography, referenceDensity, ...
    compensationDepth)
    columnMass = referenceDensity * compensationDepth;

    thickness = topography + compensationDepth;

    PrattCrustDensity = columnMass ./ thickness;
end

