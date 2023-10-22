function [R, alpha, nu] = getValues(qL, qR, gamma, lamb)

    alpha = zeros(3, 1);
    nu = zeros(3, 1);

    [uROE, HROE, aROE] = roeAverageProperties(qL, qR, gamma);
    
    % Calcula a_{j + 1/2}^k para k = 1, 2, 3
    a1 = uROE - aROE;
    a2 = uROE;
    a3 = uROE + aROE;

    % Calcula R_{j + 1/2}^k para k = 1, 2, 3
    r1 = [1; uROE - aROE; HROE - uROE * aROE];
    r2 = [1; uROE; uROE^2/2];
    r3 = [1; uROE + aROE; HROE + uROE * aROE];

    R = [r1, r2, r3];

    % Calcula alpha_{j + 1/2}^k para k = 1, 2, 3
    C1 = (gamma - 1) * ((qR(3) - qL(3)) + uROE^2 * (qR(1) - qL(1)) / 2 - uROE * (qR(2) - qL(2)))/aROE^2;
    C2 = ((qR(2) - qL(2)) - uROE * (qR(1) - qL(1))) / aROE;

    alpha(1) = (C1 - C2)/2;
    alpha(2) = (qR(1) - qL(1)) - C1;
    alpha(3) = (C1 + C2)/2;

    % Calcula nu_{j + 1/2}^k para k = 1, 2, 3
    nu(1) = lamb * a1;
    nu(2) = lamb * a2;
    nu(3) = lamb * a3;
end