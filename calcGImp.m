function [gj] =  calcGImp(qLm, qL, qR, gamma, lamb, epsValues)

    gj = zeros(1, 3);

    % Cálculo de alpha na interface entre L - 1 e L para k = 1, 2, 3
    [~, alpha_jm, nu_jm] = getValues(qLm, qL, gamma, lamb);
    a_jm = nu_jm / lamb;

    % Cálculo do alpha na interface entre L e R para k = 1, 2, 3
    [~, alpha_j, nu_j] = getValues(qL, qR, gamma, lamb);
    a_j = nu_j / lamb;

    % Cálculo de s na interface entre L e R para k = 1, 2, 3
    s = sign(alpha_j);

    % Cálculo de sigma nas interfaces em (L - 1,  L) e (L, R) para k = 1, 2, 3
    sigma_jm = zeros(1, 3);
    sigma_j = zeros(1, 3);
    for k = 1:3
        sigma_jm(k) = getSigma(a_jm(k), epsValues(k), lamb);
        sigma_j(k) = getSigma(a_j(k), epsValues(k), lamb);
    end

    % % Calcula theta no ponto j para k = 1, 2, 3
    % theta = zeros(1, 3);
    % for k = 1:3
    %     if abs(alpha_j(k)) + abs(alpha_jm(k)) ~= 0
    %         theta(k) = abs(alpha_j(k) - alpha_jm(k)) / (abs(alpha_j(k)) + abs(alpha_jm(k)));
    %     end
    % end
        
    % Calcula g no ponto j para k = 1 , 2, 3
    %w = [0 0 0];
    for k = 1:3 
        %gj(k) = (1 + w(k) * theta(k)) * s(k) * max(0, min(sigma_j(k) * abs(alpha_j(k)), s(k) * sigma_jm(k) * alpha_jm(k)));
        gj(k) = s(k) * max(0, min(sigma_j(k) * abs(alpha_j(k)), s(k) * sigma_jm(k) * alpha_jm(k)));
    end
end