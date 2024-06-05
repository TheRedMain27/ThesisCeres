function MohoUncertainty = readMohoUncertainty(filename)
    MohoUncertainty = table2array(readtable(filename, NumHeaderLines=34));
    MohoUncertainty = transpose(reshape(MohoUncertainty(:,4), [360 180]));
end

