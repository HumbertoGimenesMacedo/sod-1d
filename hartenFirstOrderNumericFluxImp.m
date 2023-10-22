function FbarMjphalf = hartenFirstOrderNumericFluxImp(qL, qR, Fj, Fjp, gamma, lamb, epsValues)

    % Obtém os valores de alpha e R na interface entre L e R para k = 1, 2, 3
    [R, alpha, nu] = getValues(qL, qR, gamma, lamb);
    a = nu / lamb;

    % Calcula a soma que está na definição do fluxo numérico
    sum = 0;
    for k = 1:3
          sum = sum - R(:, k) * getPsi(a(k), epsValues(k)) * alpha(k);
    end

    % Calcula o fluxo numérifo na interface entre L e R
    FbarMjphalf = (Fj + Fjp + sum) / 2;
end