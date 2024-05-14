function Moho = readMoho(filename)
    Moho = table2array(readtable(filename, NumHeaderLines=34));
    Moho = transpose(reshape(Moho(:,3), [360 180]));
end

