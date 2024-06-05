function getMask(filename)
    fileID = fopen(filename);
    mask = fread(fileID, "int16", "ieee-be");
    mask = rot90(reshape(mask, [21600, 10800]));
    mask = mask ~= 2;

    maskCS = GSHA(mask, 179);
    maskSC = cs2sc(maskCS);
    maskSC(1,180) = 0;
    maskSC(3,180) = 0;
    mask = GSHS(maskSC, 0.5:359.5, 0.5:179.5, 179);
    mask = mask > 0;
    save("mask360x180.mat", "mask")
end

