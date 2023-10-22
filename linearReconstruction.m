function [VL, VR, qL, qR, aL, aR]  = linearReconstruction(V, L, R, gamma)
    
    VL = V(L, :); %+ (V(L, :)- V(L - 1, :)) / 2;

    rhoL = VL(1);
    uL = VL(2);
    pL = VL(3);
    EL = pL / ((gamma - 1) * rhoL) + uL^2/2;

    VR = V(R, :);% - (V(R + 1, :) - V(R, :)) / 2;

    rhoR = VR(1);
    uR = VR(2);
    pR = VR(3);
    ER = pR / ((gamma - 1) * rhoR) + uR^2/2;

    qL = [rhoL; rhoL * uL; rhoL * EL];
    qR = [rhoR; rhoR * uR; rhoR * ER];

    aL = sqrt(gamma * pL / rhoL);
    aR = sqrt(gamma * pR / rhoR);
end