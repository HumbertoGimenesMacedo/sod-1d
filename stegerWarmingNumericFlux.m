function [F_minus, F_plus] = stegerWarmingNumericFlux(q, gamma)

    % calcula as variáveis primitivas
    rho = q(1);
    u = q(2) / q(1);
    a = sqrt(gamma * (gamma - 1) * (q(3) ./ q(1) - q(2).^2 ./ (2 * q(1).^2))); 

    % determina o número de Mach local
    M = u / a;

    factor = rho / (2 * gamma);

    % vetor de fluxo (Euler, 1D) em função das variáveis conservadas
    F =  [q(2); ...
          (3 - gamma) * q(2)^2/(2 * q(1)) + (gamma - 1) * q(3); ...
          gamma * q(2) * q(3)/q(1) - (gamma - 1) * q(2)^3/(2 * q(1)^2)];

    % determina o vetor de fluxo no ponto j
    if M > 1
        F_plus = F;
        F_minus = 0;
    elseif M < -1
        F_plus = 0;
        F_minus = F;
    else
        if M >= 0 && M <= 1

            F_plus_stegerWarming = factor * [2 * gamma * u + a - u; ...
                                            2 * (gamma - 1) * u^2 + (u + a)^2;...
                                            (gamma - 1) * u^3 + (u + a)^3/2 + (3 - gamma) * (u + a) * a^2 / (2 * (gamma - 1))];

            F_minus_stegerWarming = factor * [u - a; ...
                                              (u - a)^2; ...
                                              (u - a)^3 / 2 + (3 - gamma) * (u - a) * a^2 / (2 * (gamma - 1))];

            F_plus = F_plus_stegerWarming;
            F_minus = F_minus_stegerWarming;
        else

             F_plus_stegerWarming = factor * [u + a; ...
                                              (u + a)^2; ...
                                              (u + a)^3/2 + (3 - gamma) * (u + a) * a^2 / (2 * (gamma - 1))];

             F_minus_stegerWarming = factor * [2 * (gamma - 1) * u + u - a; ...
                                               2 * (gamma - 1) * u^2 + (u - a)^2; ...
                                               (gamma - 1) * u^3 + (u - a)^3/2 + (3 - gamma) * (u - a) * a^2 / (2 * (gamma - 1))];

            F_plus = F_plus_stegerWarming;
            F_minus = F_minus_stegerWarming;
        end

    end
end