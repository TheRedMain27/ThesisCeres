function [error, percentageError, RMS, minError, maxError] = ...
    characterizeError(measuredThickness, modelThickness)
    error = (modelThickness - measuredThickness);
    percentageError = abs(error) ./ measuredThickness * 100;
    RMS = rms(error, "all");
    minError = min(error, [], "all");
    maxError = max(error, [], "all");
end

