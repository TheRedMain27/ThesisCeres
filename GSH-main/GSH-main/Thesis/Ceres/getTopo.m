function topography = getTopo(filename)
    topography = double(imread(filename));

    topographyCS = GSHA(topography, 179);
    topographySC = cs2sc(topographyCS);
    topographySC(1,180) = 0;
    topographySC(3,180) = 0;
    topographySC(5,180) = 0;
    topography = GSHS(topographySC, 0.5:359.5, 0.5:179.5, 179);
    save("topography.mat", "topography");
end

