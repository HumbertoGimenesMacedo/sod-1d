function FLhalf = roeNumericFlux(rhoL, rhoR, uL, uR, pL, pR, HL, HR, gamma)

    EL = pL / ((gamma - 1) * rhoL) + uL^2/2;
    ER = pR / ((gamma - 1) * rhoR) + uR^2/2;

    uROE = (sqrt(rhoL) * uL + sqrt(rhoR) * uR) / (sqrt(rhoL) + sqrt(rhoR));

    HROE = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR));

    aROE = sqrt((gamma - 1) * (HROE - uROE^2/2));

    qL = [rhoL; rhoL * uL; rhoL * EL];
    qR = [rhoR; rhoR * uR; rhoR * ER];

    FL =  [qL(2); ...
           (3 - gamma) * qL(2)^2/(2 * qL(1)) + (gamma - 1) * qL(3); ...
           gamma * qL(2) * qL(3)/qL(1) - (gamma - 1) * qL(2)^3/(2 * qL(1)^2)];

    FR = [qR(2); ...
          (3 - gamma) * qR(2)^2/(2 * qR(1)) + (gamma - 1) * qR(3); ...
          gamma * qR(2) * qR(3)/qR(1) - (gamma - 1) * qR(2)^3/(2 * qR(1)^2)];


    L1 = uROE - aROE;
    L2 = uROE;
    L3 = uROE + aROE;


    alpha2 = (gamma - 1) * ((HROE - uROE^2) * (rhoR - rhoL) + ...
                            uROE * (rhoR * uR - rhoL * uL) - ...
                            (rhoR * ER - rhoL * EL)) / aROE^2;

    alpha1 = (aROE * (rhoR - rhoL) - aROE * alpha2 - (rhoR * uR - rhoL * uL) + ...
             uROE * (rhoR - rhoL)) / (2 * aROE);

    alpha3 = (aROE * (rhoR - rhoL) - aROE * alpha2 + (rhoR * uR - rhoL * uL) - ...
             uROE * (rhoR - rhoL)) / (2 * aROE);

    r1 = [1; uROE - aROE; HROE - uROE * aROE];
    r2 = [1; uROE; uROE^2/2];
    r3 = [1; uROE + aROE; HROE + uROE * aROE];

    sum = alpha1 * abs(L1) * r1 + alpha2 * abs(L2) * r2 + alpha3 * abs(L3) * r3;

    FLhalf = (FL + FR - sum) / 2;
end