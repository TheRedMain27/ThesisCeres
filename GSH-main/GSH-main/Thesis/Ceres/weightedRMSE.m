function rmse = weightedRMSE(model, measurement, weights)
    error = (model - measurement) .* weights;
    rmse = sqrt(mean(error .^ 2, "all"));
end

