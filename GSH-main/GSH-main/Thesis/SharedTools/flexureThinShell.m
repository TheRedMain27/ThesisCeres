function FlexureCrust = flexureThinShell(airyCrust, bulkModulus,...
    shearModulus, elasticThickness, mantleDensity, crustDensity, g0, ...
    radius)
    
    Poisson = (3 * bulkModulus - 2 * shearModulus) /...
        (6 * bulkModulus + 2 * shearModulus);
    YoungsModulus = 3 * bulkModulus * (1 - 2 * Poisson);
    
    AiryCS = GSHA(airyCrust, 179);
    AirySC = cs2sc(AiryCS);
    
    n = 1:size(AirySC,1);
    
    flexuralRigidity = YoungsModulus * elasticThickness^3 / ...
        (12*(1-Poisson^2));
    
    commonFactor = flexuralRigidity / (mantleDensity - crustDensity) / g0;
    firstTerm = (n .* (n + 1) - 2) .^ 2 ./ (1 - (1 - Poisson) ./ ...
        (n .* (n + 1))) ./ radius .^ 4;
    secondTermFirstFactor = 12 .* (1 - Poisson .^ 2) ./ ...
        elasticThickness .^ 2 ./ radius .^ 2;
    secondTermSecondFactor = (1 - 2 ./ (n .* (n + 1))) ./ ...
        (1 - (1 - Poisson) ./ (n .* (n + 1)));
    
    PHI = (1 + commonFactor .* (firstTerm + secondTermFirstFactor .* ...
        secondTermSecondFactor)) .^ (-1);
    
    AirySC_Flexure = zeros(size(AirySC));
    
    for m = 1:size(AirySC,2)
        AirySC_Flexure(:,m) = AirySC(:,m).*PHI';
    end
    FlexureCrust = GSHS(AirySC_Flexure, 0.5:359.5, 0.5:179.5, 179);
end