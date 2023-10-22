function [u, H, a] = averageProperties(qL, qR, gamma)
    uL = qL(2) / qL(1);
    uR = qR(2) / qR(1);

    u = (uL + uR) / 2;

    aL = sqrt(gamma * (gamma - 1) * (qL(3) ./ qL(1) - qL(2).^2 ./ (2 * qL(1).^2)));
    aR = sqrt(gamma * (gamma - 1) * (qR(3) ./ qR(1) - qR(2).^2 ./ (2 * qR(1).^2)));

    a = (aL + aR) / 2;

    HL = (qL(2) / qL(1))^2 / 2 + aL^2 / (gamma - 1);
    HR = (qR(2) / qR(1))^2 / 2 + aR^2 / (gamma - 1);

    H = (HL + HR) / 2;
end