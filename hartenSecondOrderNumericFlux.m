function FbarMjphalf = hartenSecondOrderNumericFlux(qLm, qL, qR, qRp, Fj, Fjp, gamma, lamb, epsValues)

    % Obtém os valores de R, alpha e nu na interface entre L e R para k = 1, 
    % 2, 3
    [R, alpha, nu] = getValues(qL, qR, gamma, lamb);

    % Obtém os valores de g nos pontos j para k = 1, 2, 3
    gj = calcG(qLm, qL, qR, gamma, lamb, epsValues);

    % Obtém os valores de g no ponto j + 1 para k = 1, 2, 3
    gjp = calcG(qL, qR, qRp, gamma, lamb, epsValues);

    % Obtém os valores de gamma na interface entre L e R para k = 1, 2, 3
    gam = zeros(1, 3);
    for k= 1:3
        if alpha(k) ~=0
            gam(k) = (gjp(k) - gj(k)) / alpha(k);
        end
    end

    % Calcula a soma que está na definição do fluxo numérico
    sum = 0;
    for k = 1:3
          sum = sum + R(:, k) * (gj(k) + gjp(k) - getPsi(nu(k) + gam(k), epsValues(k)) * alpha(k)) / lamb;
    end

    % Calcula o fluxo numérico na interface entre L e R
    FbarMjphalf = (Fj + Fjp + sum) / 2;
end