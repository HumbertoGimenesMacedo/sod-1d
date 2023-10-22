function [uROE, HROE, aROE] = roeAverageProperties(qL, qR, gamma)
  %  global averageType;
    
   % if averageType == 1
        % rhoL = qL(1);
        % uL = qL(2) / qL(1);
        % 
        % rhoR = qR(1);
        % uR = qR(2) / qR(1);
        % 
        % aL = sqrt(gamma * (gamma - 1) * (qL(3) ./ qL(1) - qL(2).^2 ./ (2 * qL(1).^2)));
        % aR = sqrt(gamma * (gamma - 1) * (qR(3) ./ qR(1) - qR(2).^2 ./ (2 * qR(1).^2)));
        % 
        % HL = (qL(2) / qL(1))^2 / 2 + aL^2 / (gamma - 1);
        % HR = (qR(2) / qR(1))^2 / 2 + aR^2 / (gamma - 1);
        % 
        % uROE = (sqrt(rhoL) * uL + sqrt(rhoR) * uR) / (sqrt(rhoL) + sqrt(rhoR));
        % 
        % HROE = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR));
        % 
        % aROE = sqrt((gamma - 1) * (HROE - uROE^2/2));
   % elseif averageType == 2

        uL = qL(2) / qL(1);
        uR = qR(2) / qR(1);

        uROE = (uL + uR) / 2;

        aL = sqrt(gamma * (gamma - 1) * (qL(3) ./ qL(1) - qL(2).^2 ./ (2 * qL(1).^2)));
        aR = sqrt(gamma * (gamma - 1) * (qR(3) ./ qR(1) - qR(2).^2 ./ (2 * qR(1).^2)));

        aROE = (aL + aR) / 2;

        HL = (qL(2) / qL(1))^2 / 2 + aL^2 / (gamma - 1);
        HR = (qR(2) / qR(1))^2 / 2 + aR^2 / (gamma - 1);

        HROE = (HL + HR) / 2;
    %end

% roeAverageProperties(VL, VR, qL, qR, aL, aR, H, gamma)
    % HL = H(qL, aL); 
    % HR = H(qR, aR);
    % 
    % rhoL = VL(1);
    % uL = VL(2);
    % 
    % rhoR = VR(1);
    % uR = VR(2);
    % 
    % uROE = (sqrt(rhoL) * uL + sqrt(rhoR) * uR) / (sqrt(rhoL) + sqrt(rhoR));
    % 
    % HROE = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR));
    % 
    % aROE = sqrt((gamma - 1) * (HROE - uROE^2/2));
end