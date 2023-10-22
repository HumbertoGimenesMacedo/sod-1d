function FbarMjphalf = hartenSecondOrderNumericFluxImp(qLm, qL, qR, qRp, Fj, Fjp, gamma, lamb, epsValues)

    % Obtém os valores de alpha e R na interface entre L e R para k = 1, 2, 3
    [R, alpha, nu] = getValues(qL, qR, gamma, lamb);
    a = nu / lamb;
    % Obtém os valores de g nos pontos j e j + 1
    gj = calcGImp(qLm, qL, qR, gamma, lamb, epsValues);
    gjp = calcGImp(qL, qR, qRp, gamma, lamb, epsValues);

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
          sum = sum + R(:, k) * (gj(k) + gjp(k) - getPsi(a(k) + gam(k), epsValues(k)) * alpha(k));
    end

    % Calcula o fluxo numérifo na interface entre L e R
    FbarMjphalf = (Fj + Fjp + sum) / 2;
end