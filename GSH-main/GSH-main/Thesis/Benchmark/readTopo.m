function topography = readTopo(filename)
    fileID = fopen(filename);
    topography = fread(fileID, "int16", "ieee-be");
    topography = rot90(reshape(topography, [4320, 2160])) + 6371e3;

    topographyCS = GSHA(topography, 179);
    topographySC = cs2sc(topographyCS);
    topographySC(1,180) = 0;
    topographySC(3,180) = 0;
    topography = GSHS(topographySC, 0.5:359.5, 0.5:179.5, 179);
end

