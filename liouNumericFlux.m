function FbarLhalf = liouNumericFlux(uL, uR, pL, pR, HL, HR, phiL, phiR, alpha, beta, gamma)

    aStarL = sqrt(2 * (gamma - 1) * HL / (gamma + 1)); 
    
    aTildeL = aStarL * min(1, aStarL / abs(uL)); 
    
    aStarR = sqrt(2 * (gamma - 1)  * HR / (gamma + 1)); 
    
    aTilderR = aStarR * min(1, aStarR / abs(uR)); 
    
    aHalfjp = min(aTildeL, aTilderR); 
    
    ML = uL / aHalfjp; 
    MR = uR / aHalfjp; 
    
    if abs(ML) >= 1
        ML_plus = (ML + abs(ML))/2; 
        pL_plus = (1 + sign(ML))/2; 
    else
        ML_plus = (ML + 1)^2/4 + beta * (ML^2 - 1)^2; 
        pL_plus = (ML + 1)^2 * (2 - ML)/4 + alpha * ML * (ML^2-1)^2; 
    end
    
    if abs(MR) >= 1
        MR_minus = (MR - abs(MR)) / 2; 
        pR_minus = (1 - sign(MR))/2; 
    else
        MR_minus = -(MR - 1)^2/4 - beta * (MR^2 - 1)^2; 
        pR_minus = (MR - 1)^2 * (2 + MR)/4 - alpha * MR * (MR^2-1)^2; %certp
    end
    
    MLhalf = ML_plus + MR_minus; 
    
    pLhalf = pL_plus * pL + pR_minus * pR; 
    
    MLhalf_plus = (MLhalf + abs(MLhalf)) / 2; 
    MLhalf_minus = (MLhalf - abs(MLhalf)) / 2; 
    
    FbarLhalf = aHalfjp * (MLhalf_plus * phiL + MLhalf_minus * phiR) + [0; pLhalf; 0]; 

end