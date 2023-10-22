function J = getJFirstOrder(qL, qR, gamma, lamb, epsValues, type)

    [R, ~, nu] = getValues(qL, qR, gamma, lamb);
    a = nu / lamb;

    D = zeros(3, 3);
    C = zeros(1, 3);

    for k= 1:3
        if type == 1 % plus
            C(k) = (getPsi(a(k), epsValues(k)) + a(k)) / 2;
        elseif type == 2 % minus
            C(k) = (getPsi(a(k) , epsValues(k)) - a(k)) / 2;
        end
        D(k, k) = C(k);
    end

    [uROE, HROE, aROE] = roeAverageProperties(qL, qR, gamma);
    J = R * D * getRinv(uROE, aROE, HROE);
end