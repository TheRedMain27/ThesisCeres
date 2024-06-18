function correction = bouguerCorrection(topography, density, G)
    correction = 2 * pi * G * density * topography;
end

