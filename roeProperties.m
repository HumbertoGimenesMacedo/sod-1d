function [VL, VR, qL, qR, aL, aR]  = roeProperties(V, L, R, gamma, fluxLimiter)

    rho = V(:, 1);
    u = V(:, 2);
    p = V(:, 3);

    Z_plus = [fluxLimiter((rho(R) - rho(L)) / (rho(L) - rho(L - 1))),...
              fluxLimiter((u(R) - u(L)) / (u(L) - u(L - 1))),...
              fluxLimiter((p(R) - p(L)) / (p(L) - p(L - 1)))];
    
    Z_minus = [fluxLimiter((rho(R) - rho(L)) / (rho(R + 1) - rho(R))),...
               fluxLimiter((u(R) - u(L)) / (u(R + 1) - u(R))),...
               fluxLimiter((p(R) - p(L)) / (p(R + 1) - p(R)))];
    
    Z_plus(isnan(Z_plus)) = 1;
    Z_minus(isnan(Z_minus)) = 1;
    
    VL = V(L, :) + Z_plus .* (V(L, :)- V(L - 1, :)) / 2;
    
    rhoL = VL(1);
    uL = VL(2);
    pL = VL(3);
    EL = pL / ((gamma - 1) * rhoL) + uL^2/2;
    
    VR = V(R, :) - Z_minus .* (V(R + 1, :) - V(R, :)) / 2;
    
    rhoR = VR(1);
    uR = VR(2);
    pR = VR(3);
    ER = pR / ((gamma - 1) * rhoR) + uR^2/2;
    
    qL = [rhoL; rhoL * uL; rhoL * EL];
    qR = [rhoR; rhoR * uR; rhoR * ER];
    
    aL = sqrt(gamma * pL / rhoL);
    aR = sqrt(gamma * pR / rhoR);
end