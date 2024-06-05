function error = selectiveRMSE(model, measurement, mask)
    error = model - measurement;
    error = sqrt(mean(error(mask) .^ 2));
end

