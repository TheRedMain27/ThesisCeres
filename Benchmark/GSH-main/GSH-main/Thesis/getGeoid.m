function geoid = getGeoid()
    latitudes = reshape(repmat(flip(-89.5:89.5), 360, 1), [1, 64800]);
    longitudes = repmat(-179.5:179.5, 1, 180);
    
    % Use aerospace toolbox, needs geoid data for aetoolbox package
    % egm2008 is bugged, so use egm96
    geoid = geoidheight(latitudes, longitudes);
    geoid = transpose(reshape(geoid, [360, 180]));
end

