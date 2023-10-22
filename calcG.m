function [g] =  calcG(qLm, qL, qR, gamma, lamb, epsValues)

    % Vetor que armazena os valores de g no ponto j para k = 1, 2, 3
    g = zeros(3, 1);

    % Obtém os valores de g̃ na interface entre L - 1 e L para k = 1, 2, 3
    gt_jmHalf = gTilde(qLm, qL, gamma, lamb, epsValues);
    
    % Obtém os valores de g̃ na interface entre L e R para k = 1, 2, 3
    gt_jpHalf = gTilde(qL, qR, gamma, lamb, epsValues);
    
    % Cálculo de s na interface entre L e R para k = 1, 2, 3
    s = sign(gt_jpHalf);

    for k = 1:3
        g(k) = s(k) * max(0, min(abs(gt_jpHalf(k)), gt_jmHalf(k) * s(k)));
    end
end