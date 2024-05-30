function gravityField = getEarthGravity()
    latitudes = reshape(repmat(flip(-89.5:89.5), 360, 1), [1, 64800]);
    longitudes = repmat(-179.5:179.5, 1, 180);
    % use radius of circular orbit 500km altitude
    r = repmat(6871e3, 1, 64800);

    [xECEF, yECEF, zECEF] = sph2cart(longitudes, latitudes, r);
    
    % Use aerospace toolbox
    [gx, gy, gz] = gravitysphericalharmonic(...
        transpose([xECEF; yECEF; zECEF]));
    
    [~, ~, radialGravity] = cart2sph(gx, gy, gz);

    gravityField = transpose(reshape(-radialGravity, [360, 180]));

    gravityCS = GSHA(gravityField, 179);
    gravitySC = cs2sc(gravityCS);
    gravitySC(1,180) = 0;
    gravitySC(3,180) = 0;
    % gravitySC(2,179:181) = zeros(1, 3);
    % gravitySC(3,178:182) = zeros(1, 5);
    % gravitySC(4,177:183) = zeros(1, 7);
    % gravitySC(5,176:184) = zeros(1, 9);
    % gravitySC(6,175:185) = zeros(1, 11);
    % gravitySC(7,174:186) = zeros(1, 13);
    % gravitySC(8,173:187) = zeros(1, 15);
    % gravitySC(9,172:188) = zeros(1, 17);
    % gravitySC(10,171:189) = zeros(1, 19);

    gravityField = GSHS(gravitySC, 0.5:359.5, 0.5:179.5, 179);
end
