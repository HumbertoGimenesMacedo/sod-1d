function [F_minus, F_plus] = vanLeerNumericFlux(q, gamma)

    % calcula as variáveis primitivas
    rho = q(1);
    a = sqrt(gamma * (gamma - 1) * (q(3) ./ q(1) - q(2).^2 ./ (2 * q(1).^2))); 

    % determina o número de Mach local
    M = q(2) / (q(1) * a);

    % vetor de fluxo (Euler, 1D) em função das variáveis conservadas
    F =  [q(2); ...
          (3 - gamma) * q(2)^2/(2 * q(1)) + (gamma - 1) * q(3); ...
          gamma * q(2) * q(3)/q(1) - (gamma - 1) * q(2)^3/(2 * q(1)^2)];
    
    % define os vetores de fluxo "mais" e "menos" de van Leer
    f_mass_plus = rho * a * ((M + 1) / 2)^2;
    f_mass_minus = - rho * a * ((M - 1) / 2)^2;
    
    F_plus_vanLeer =  f_mass_plus * ...
                      [1; 
                      ((gamma - 1) * q(2) / q(1) + 2 * a) / gamma; 
                      ((gamma - 1) * q(2) / q(1) + 2 * a)^2 / (2 * (gamma^2 - 1))];
    
    F_minus_vanLeer = f_mass_minus * ...
                      [1;
                       ((gamma - 1) * q(2) / q(1) - 2 * a) / gamma;
                       ((gamma - 1) * q(2) / q(1) - 2 * a)^2 / (2 * (gamma^2 - 1))];

    % determina o vetor de fluxo no ponto j
    if M >= 1
        F_plus = F;
        F_minus = 0;
    elseif M <= -1
        F_plus = 0;
        F_minus = F;
    else
        F_plus = F_plus_vanLeer;
        F_minus = F_minus_vanLeer;
    end
end