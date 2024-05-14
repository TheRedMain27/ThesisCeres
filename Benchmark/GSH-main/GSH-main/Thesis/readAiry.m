function AiryMatrix = readAiry(filename)
    fileID = fopen(filename);
    AiryMatrix = fread(fileID, "float64");
    AiryMatrix = transpose(reshape(AiryMatrix, [4320, 2160]));
end