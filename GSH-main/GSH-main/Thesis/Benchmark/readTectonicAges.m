function [fullCrustMap, mesozoicCenozoicMap, paleozoicMap, ...
    proterozoicMap, archeanMap] = readTectonicAges(filename)
    opts = detectImportOptions(filename);
    opts = setvartype(opts,'string');
    dataTable = readtable(filename,opts);
    dataTable = dataTable(:,"Thermo_tectonicAge");
    ageData = table2array(dataTable);
    ageMap = transpose(reshape(ageData, [360 180]));
    mesozoicCenozoicMap = ageMap == "Meso- and Cenozoic";
    paleozoicMap = ageMap == "Paleozoic";
    proterozoicMap = ageMap == "Late Proterozoic" | ...
        ageMap == "Middle Proterozoic" | ageMap == "Early Proterozoic";
    archeanMap = ageMap == "Archean";
    fullCrustMap = logical(mesozoicCenozoicMap + paleozoicMap + ...
        proterozoicMap + archeanMap);
end

