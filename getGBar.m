function [gbb, theta] = getGBar(s, alpha_jmHalf, alpha_jpHalf, nu_jmHalf, nu_jpHalf, epsValues)
    gbb = zeros(1, 3);
    theta = zeros(1, 3);
    for k = 1:3
        sigma_jmHalf = (1 - getPsi(nu_jmHalf(k), epsValues(k))) / 2;
        sigma_jpHalf = (1 - getPsi(nu_jpHalf(k), epsValues(k))) / 2;

        gbb(k) = s(k) * max(0, min(s(k) * sigma_jmHalf * alpha_jmHalf(k), sigma_jpHalf * abs(alpha_jpHalf(k))));
    
        if abs(alpha_jpHalf(k)) + abs(alpha_jmHalf(k)) ~= 0
            theta(k) = abs(alpha_jpHalf(k) - alpha_jmHalf(k)) / (abs(alpha_jpHalf(k)) + abs(alpha_jmHalf(k)));
        end
    end
end