function FlexureCrust = flexureInfinitePlate(airyCrust, bulkModulus,...
    shearModulus, elasticThickness, mantleDensity, crustDensity, g0, ...
    radius)
    
    addpath('..\Tools');
    
    Poisson = (3 * bulkModulus - 2 * shearModulus) /...
        (6 * bulkModulus + 2 * shearModulus);
    YoungsModulus = 3 * bulkModulus * (1 - 2 * Poisson);
    
    AiryCS = GSHA(airyCrust, 179);
    AirySC = cs2sc(AiryCS);
    
    n = 1:size(AirySC,1);
    
    flexuralRigidity = YoungsModulus * elasticThickness^3 / ...
        (12*(1-Poisson^2));
    PHI = (1 + flexuralRigidity / (mantleDensity - crustDensity) / g0 * ...
        ((n + 1) / radius) .^ 4) .^ (-1);
    
    AirySC_Flexure = zeros(size(AirySC));
    
    for m = 1:size(AirySC,2)
        AirySC_Flexure(:,m) = AirySC(:,m).*PHI';
    end
    FlexureCrust = GSHS(AirySC_Flexure, 0.5:359.5, 0.5:179.5, 179);
end