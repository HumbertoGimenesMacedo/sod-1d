function J = getJ(qLm, qL, qR, qRp, gamma, lamb, epsValues, type)

    [R, alpha, nu] = getValues(qL, qR, gamma, lamb);
    a = nu / lamb;

    gj = calcGImp(qLm, qL, qR, gamma, lamb, epsValues);
    gjp = calcGImp(qL, qR, qRp, gamma, lamb, epsValues); 

    gam = zeros(1, 3);
    D = zeros(3, 3);
    C = zeros(1, 3);

    for k= 1:3
        if alpha(k) ~=0
            gam(k) = (gjp(k) - gj(k))/alpha(k);
        end

        if type == 1 % plus
            C(k) = (getPsi(a(k) + gam(k), epsValues(k)) + a(k) + gam(k)) / 2;
        elseif type == 2 % minus
            C(k) = (getPsi(a(k) + gam(k), epsValues(k)) - a(k) - gam(k)) / 2;
        end
        D(k, k) = C(k);
    end

    [uROE, HROE, aROE] = roeAverageProperties(qL, qR, gamma);
    J = R * D * getRinv(uROE, aROE, HROE);
end