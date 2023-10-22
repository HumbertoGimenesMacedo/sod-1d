function FbarMjphalf = hartenNumericFlux(qL, qR, Fj, Fjp, gamma, lamb, epsValues)

    [R, alpha, nu] = getValues(qL, qR, gamma, lamb);

    sum = 0;
    for k= 1:3
        sum = sum - R(:, k) * getPsi(nu(k), epsValues(k)) * alpha(k) / lamb;
    end

    FbarMjphalf = (Fj + Fjp + sum) / 2;
end