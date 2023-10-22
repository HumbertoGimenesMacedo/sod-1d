function gt = gTilde(qL, qR, gamma, lamb, epsValues)

    gt = zeros(1, 3);
    [~, alpha_jpHalf, nu_jpHalf] = getValues(qL, qR, gamma, lamb);
    for k = 1:3
        gt(k) = (getPsi(nu_jpHalf(k), epsValues(k)) - nu_jpHalf(k)^2) * alpha_jpHalf(k) / 2;
    end

    % [uROE, ~, aROE] = roeAverageProperties(qL, qR, gamma);
    % 
    % % Calcula a_{j + 1/2}^k para k = 1, 2, 3
    % a1 = uROE - aROE;
    % a2 = uROE;
    % a3 = uROE + aROE;
    % 
    % % Calcula alpha_{j + 1/2}^k para k = 1, 2, 3
    % C1 = (gamma - 1) * ((qR(3) - qL(3)) + uROE^2 * (qR(1) - qL(1)) / 2 - uROE * (qR(2) - qL(2)))/aROE^2;
    % C2 = ((qR(2) - qL(2)) - uROE * (qR(1) - qL(1))) / aROE;
    % 
    % alpha1 = (C1 - C2)/2;
    % alpha2 = (qR(1) - qL(1)) - C1;
    % alpha3 = (C1 + C2)/2;
    % 
    % % Calcula nu_{j + 1/2}^k para k = 1, 2, 3
    % nu1 = lamb * a1;
    % nu2 = lamb * a2;
    % nu3 = lamb * a3;

    

    % Calcula g_tilde_{j + 1}^k para k = 1, 2, 3
    
    % gt1 = (getPsi(nu1, epsValues(1)) - nu1^2) * alpha1 / 2;
    % gt2 = (getPsi(nu2, epsValues(2)) - nu2^2) * alpha2 / 2;
    % gt3 = (getPsi(nu3, epsValues(3)) - nu3^2) * alpha3 / 2;

    %gt = [gt1, gt2, gt3];
end