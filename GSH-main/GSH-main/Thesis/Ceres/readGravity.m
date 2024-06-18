function [gravity, SH_coefficients] = readGravity(filename, SHbounds, ...
    gravitationalParameter, radius)

    SH_coefficients = readtable(filename);
    SH_coefficients = table2array(SH_coefficients(:,1:4));
    SH_coefficients = [0 0 0 0; SH_coefficients];
    SH_coefficients = sortrows(SH_coefficients,2);
    SH_coefficients(3,3) = 0;
    SH_coefficients(5,3) = 0;
    
    latLim =    [-89.5 89.5 1];
    lonLim =    [0.5 359.5 1];
    lon = lonLim(1):lonLim(3):lonLim(2);
    lat = latLim(1):latLim(3):latLim(2);
    Lon = repmat(lon,length(lat),1);
    Lat = repmat(lat',1,length(lon));

    data = gravityModule(Lat, Lon, 470e3, SHbounds, ...
        SH_coefficients, radius, gravitationalParameter);
    gravity = flip(data.vec.R, 1);
end

