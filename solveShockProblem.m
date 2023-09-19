function [x, rho, u, E, p, totalTime] = solveShockProblem(METHOD, order, rho_l, u_l, p_l, rho_r, u_r, p_r)
    %% Definição das variáveis

    startTime = tic;
    
    % variáveis que identificam o esquema numérico
    LAX_WENDROFF = 1;
    MAC_COMARCK = 2;
    EXP_BEAM_WARMING = 3;
    IMP_BEAM_WARMING = 4;
    EXP_STEGER_WARMING = 5;
    AUSM_PLUS = 6;
    VAN_LEER = 7;
    IMP_STEGER_WARMING = 8;
    ROE = 9;

    % número de pontos interiores da malha
    M = 2998;
    
    % comprimento do tubo
    L_ref = 1;
    L = 2000 / L_ref;
    
    % espaçamento da malha
    dx = L / (M + 1);

    % malha
    x = (0:M + 1)' * dx;

    % razão de calores específicos
    gamma = 1.4;

    % instante de tempo
    t = 0;

    % instante final de tempo
    tEnd = 1;

    % define o limitador de fluxo
    fluxLimiter = @(r) max(0, min(1, r));

    % matriz contendo a solução em cada instante de tempo
    q = zeros(M + 2, 3);

    %% Condição de contorno + inicial

    % define as grandezas de referência
    rho_atm = 1.225;
    p_atm = 101325;

    rho_ref = rho_atm;
    u_ref = sqrt(gamma * p_atm/rho_atm);
    p_ref = rho_ref * u_ref^2;

    pressureRatio = round(p_l/p_r);

    % adimensionaliza as grandezas utilizando as referências
    rho_r = rho_r / rho_ref;
    u_r = u_r / u_ref;
    p_r = p_r / p_ref;
    
    rho_l = rho_l / rho_ref;
    u_l = u_l / u_ref;
    p_l = p_l / p_ref;
    
    % condições de contorno à esquerda
    q(1, 1) = rho_l; % densidade
    q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
    q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total

    % condições de contorno à direita
    q(M + 2, 1) = rho_r; % densidade
    q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
    q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total

    % condição inicial (primeira metade)
    for j = 2:ceil(M/2)
        q(j, 1) = rho_l; % densidade
        q(j, 2) = rho_l * u_l; % quantidade de movimento na direção x
        q(j, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
    end
    
    % condição inicial (segunda metade)
    for j = ceil(M/2) + 1:M + 1
        q(j, 1) = rho_r; % densidade
        q(j, 2) = rho_r * u_r; % quantidade de movimento na direção x
        q(j, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
    end

    % velocidade do som empregando as variáveis conservadas (n = 0)
    a0 = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2)));
    
    % calcula a velocidade na direção x (n = 0)
    u0 = q(:, 2) ./ q(:, 1);

    %% Aplicação do esquema numérico

    switch(METHOD)

        case ROE

            % entalpia total
            H = @(q, a) (q(2) / q(1))^2 / 2 + a^2 / (gamma - 1);

            % determina o passo no tempo
            if order == 1
                % determina o passo temporal
                if ismember(pressureRatio, [5, 10])
                    CFL = 0.3;
                elseif pressureRatio == 20
                    CFL = 0.1;
                elseif pressureRatio == 50
                    CFL = 0.05;
                elseif pressureRatio == 100
                    CFL = 0.01;
                end
            elseif order == 2
                 % determina o passo temporal
                if ismember(pressureRatio, [5, 10, 20])
                    CFL = 0.1;
                elseif pressureRatio == 50
                    CFL = 0.05;
                elseif pressureRatio == 100
                    CFL = 0.01;
                end
            end
            eigenvalues = [u0, u0 + a0, u0 - a0];
            dt = CFL * dx / max(eigenvalues(:));

            % coeficiente do esquema numérico
            lamb = dt / dx;

            while 1

                % efetua um incremento temporal
                t = t + dt;
   
                if t * (L_ref / u_ref) > 1 
                    break
                end

                % mostra a solução corrente
                %show(x, q, gamma);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                % calcula as variáveis primitivas
                rho = q(:, 1);
                u = q(:, 2) ./ rho;
                p = (gamma - 1) * (q(:, 3) - q(:, 2).^2 ./ (2 * q(:, 1)));
                a = sqrt(gamma * p ./ rho);

                % constrói o vetor de variáveis primitivas
                V = [rho, u, p];

                % calcula a solução no próximo instante de tempo
                for j = 3:M 
                    if order == 1
                        % calcula o fluxo numérico de Liou na interface j + 1/2
                        Fbarjphalf = roeNumericFlux(rho(j), rho(j + 1), ...
                                                    u(j), u(j + 1), ...
                                                    p(j), p(j + 1), ...
                                                    H(q(j, :), a(j)), H(q(j + 1, :), a(j + 1)), gamma);

                        % calcula o fluxo numérico de Liou na interface j - 1/2
                        Fbarjmhalf = roeNumericFlux(rho(j - 1), rho(j), ...
                                                    u(j - 1), u(j), ...
                                                    p(j - 1), p(j), ...
                                                    H(q(j - 1, :), a(j - 1)), H(q(j, :), a(j)), gamma);
                    elseif order == 2
                        
                        % calcula o fluxo numérico de Liou na interface j - 1/2
                        [VL, VR, qL, qR, aL, aR] = roeProperties(V, j - 1, j, gamma, fluxLimiter);
                        Fbarjmhalf = roeNumericFlux(VL(1), VR(1), ...
                                                    VL(2), VR(2), ...
                                                    VL(3), VR(3), ...
                                                    H(qL, aL), H(qR, aR), gamma);

                        % calcula o fluxo numérico de Liou na interface j + 1/2
                        [VL, VR, qL, qR, aL, aR] = roeProperties(V, j, j + 1, gamma, fluxLimiter);
                        Fbarjphalf = roeNumericFlux(VL(1), VR(1), ...
                                                    VL(2), VR(2), ...
                                                    VL(3), VR(3), ...
                                                    H(qL, aL), H(qR, aR), gamma);
        
                    end

                    next_q(j, :) = q(j, :) - lamb * (Fbarjphalf - Fbarjmhalf)';
                end

                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                % aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
        
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
            end

        case VAN_LEER

            % determina o passo no tempo
            CFL = 0.1;
            eigenvalues = [u0, u0 + a0, u0 - a0];
            dt = CFL * dx / max(eigenvalues(:));

            % coeficiente do esquema numérico
            lamb = dt / dx;

            while 1

                % efetua um incremento temporal
                t = t + dt;
   
                if t * (L_ref / u_ref) > 1 
                    break
                end

                % mostra a solução corrente
                %show(x, q, gamma);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                % armazena o valor da densidade em todos os pontos da malha
                rho = q(:, 1);
                 
                % calcula a solução no próximo instante de tempo
                for j = 3:M 

                    % calcula o vetor de fluxo numérico de van Leer 
                    % em quatro pontos consecutivos da malha
                    [~, Fjmm_plus] = vanLeerNumericFlux(q(j - 2, :), gamma);
                    [~, Fjm_plus] = vanLeerNumericFlux(q(j - 1, :), gamma);
                    [Fj_minus, Fj_plus] = vanLeerNumericFlux(q(j , :), gamma);
                    [Fjp_minus, ~] = vanLeerNumericFlux(q(j + 1, :), gamma);
                    [Fjpp_minus, ~] = vanLeerNumericFlux(q(j + 2, :), gamma);
                   
                    % verifica em qual ordem o esquema deve ser aplicado
                    if order == 1
                        F_bar_jphalf_plus = Fj_plus;
                        F_bar_jmhalf_plus = Fjm_plus ;
                        F_bar_jphalf_minus = Fjp_minus;
                        F_bar_jmhalf_minus = Fj_minus;
                    elseif order == 2

                        % calcula os valores do limitador
                        Z_plus_1 = fluxLimiter((rho(j + 1) - rho(j)) / (rho(j) - rho(j - 1)));
                        Z_minus_1 = fluxLimiter((rho(j + 1) - rho(j)) / (rho(j + 2) - rho(j + 1)));
                        
                        Z_plus_1(isnan(Z_plus_1)) = 1;
                        Z_minus_1(isnan(Z_minus_1)) = 1;
                        
                        Z_plus_2 = fluxLimiter((rho(j) - rho(j - 1)) / (rho(j - 1) - rho(j - 2)));
                        Z_minus_2 = fluxLimiter((rho(j) - rho(j - 1)) / (rho(j + 1) - rho(j)));
                        
                        Z_plus_2(isnan(Z_plus_2)) = 1;
                        Z_minus_2(isnan(Z_minus_2)) = 1;
            
                        F_bar_jphalf_plus = Fj_plus + Z_plus_1 .* (Fj_plus - Fjm_plus) / 2;
                        F_bar_jphalf_minus = Fjp_minus - Z_minus_1 .* (Fjpp_minus - Fjp_minus) / 2; 
                        F_bar_jmhalf_plus = Fjm_plus + Z_plus_2 .* (Fjm_plus - Fjmm_plus) / 2;
                        F_bar_jmhalf_minus = Fj_minus - Z_minus_2 .* (Fjp_minus - Fj_minus) / 2;
                    end

                    next_q(j, :) = q(j, :) - lamb * (F_bar_jphalf_plus - F_bar_jmhalf_plus)'  ...
                                - lamb * (F_bar_jphalf_minus - F_bar_jmhalf_minus)';

                end

                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                % aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
        
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
            end

        case AUSM_PLUS

            % entalpia total
            H = @(q, a) (q(2) / q(1))^2 / 2 + a^2 / (gamma - 1);

            % vetor de escalares passivos 
            phi = @(q, a) [q(1); q(2); q(1) * H(q, a)];

            if order == 1
                % determina o passo temporal
                if ismember(pressureRatio, [5, 10])
                    CFL = 0.3;
                elseif pressureRatio == 20
                    CFL = 0.1;
                elseif pressureRatio == 50
                    CFL = 0.05;
                elseif pressureRatio == 100
                    CFL = 0.01;
                end
            elseif order == 2
                 % determina o passo temporal
                if ismember(pressureRatio, [5, 10, 20])
                    CFL = 0.1;
                elseif pressureRatio == 50
                    CFL = 0.05;
                elseif pressureRatio == 100
                    CFL = 0.01;
                end
            end

            dt = CFL * dx / max(abs(u0) + a0);

            % parâmetros do esquema de Liou
            alpha = 3 / 16;
            beta = 1/8;

            % coeficiente do esquema numérico
            lamb = dt / dx;

            while 1

                % efetua um incremento temporal
                t = t + dt;

                if t * (L_ref / u_ref) > 1 
                    break
                end

                % mostra a solução corrente
                %show(x, q, gamma);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                % calcula as variáveis primitivas
                rho = q(:, 1);
                u = q(:, 2) ./ rho;
                p = (gamma - 1) * (q(:, 3) - q(:, 2).^2 ./ (2 * q(:, 1)));
                a = sqrt(gamma * p ./ rho);

                % constrói o vetor de variáveis primitivas
                V = [rho, u, p];

                for j = 3:M

                    if order == 1

                        % calcula o fluxo numérico de Liou na interface j + 1/2
                        Fbarjphalf = liouNumericFlux(u(j), u(j + 1), ...
                                                     p(j), p(j + 1), ...
                                                     H(q(j, :), a(j)), H(q(j + 1, :), a(j + 1)), ...
                                                     phi(q(j, :), a(j)), phi(q(j + 1, :), a(j + 1)), ...
                                                     alpha, beta, gamma);
    
                         % calcula o fluxo numérico de Liou na interface j - 1/2
                        Fbarjmhalf = liouNumericFlux(u(j - 1), u(j), ...
                                                     p(j - 1), p(j), ...
                                                     H(q(j - 1, :), a(j - 1)), H(q(j, :), a(j)), ...
                                                     phi(q(j - 1, :), a(j - 1)), phi(q(j, :), a(j)), ...
                                                     alpha, beta, gamma);
                    elseif order == 2
                        
                        % calcula o fluxo numérico de Liou na interface j - 1/2
                        [uL, uR, pL, pR, qL, qR, aL, aR] = liouProperties(V, j - 1, j, gamma, fluxLimiter);
                        Fbarjmhalf = liouNumericFlux(uL, uR, ...
                                                     pL, pR, ...
                                                     H(qL, aL), H(qR, aR), ...
                                                     phi(qL, aL), phi(qR, aR), ...
                                                     alpha, beta, gamma);

                        % calcula o fluxo numérico de Liou na interface j + 1/2
                        [uL, uR, pL, pR, qL, qR, aL, aR] = liouProperties(V, j, j + 1, gamma, fluxLimiter);
                        Fbarjphalf = liouNumericFlux(uL, uR, ...
                                                     pL, pR, ...
                                                     H(qL, aL), H(qR, aR), ...
                                                     phi(qL, aL), phi(qR, aR), ...
                                                     alpha, beta, gamma);  
                    end

                    next_q(j, :) = q(j, :) - lamb * (Fbarjphalf - Fbarjmhalf)';
                end

                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                % aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
        
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
            end
        
        case IMP_STEGER_WARMING

            % matriz jacobiana de fluxo em função dos autovalores
            A = @(u, a, lambda) reshape([-lambda(1).*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0)-(1.0./a.^2.*lambda(2).*(a.*u-(u.^2.*(gamma-1.0))./2.0))./2.0+(1.0./a.^2.*lambda(3).*(a.*u+(u.^2.*(gamma-1.0))./2.0))./2.0,-lambda(1).*u.*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0)-lambda(2).*((1.0./a.^2.*u)./2.0+1.0./(a.*2.0)).*(a.*u-(u.^2.*(gamma-1.0))./2.0)+lambda(3).*((1.0./a.^2.*u)./2.0-1.0./(a.*2.0)).*(a.*u+(u.^2.*(gamma-1.0))./2.0),lambda(1).*u.^2.*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0).*(-1.0./2.0)-lambda(2).*(a.*u-(u.^2.*(gamma-1.0))./2.0).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+lambda(3).*(a.*u+(u.^2.*(gamma-1.0))./2.0).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0),(1.0./a.^2.*lambda(2).*(a-u.*(gamma-1.0)))./2.0-(1.0./a.^2.*lambda(3).*(a+u.*(gamma-1.0)))./2.0+1.0./a.^2.*lambda(1).*u.*(gamma-1.0),1.0./a.^2.*lambda(1).*u.^2.*(gamma-1.0)+(1.0./a.^2.*lambda(2).*(a+u).*(a+u-gamma.*u))./2.0+(1.0./a.^2.*lambda(3).*(a-u).*(a-u+gamma.*u))./2.0,lambda(2).*(a-u.*(gamma-1.0)).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)-lambda(3).*(a+u.*(gamma-1.0)).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+(1.0./a.^2.*lambda(1).*u.^3.*(gamma-1.0))./2.0,(1.0./a.^2.*(gamma-1.0).*(lambda(1).*-2.0+lambda(2)+lambda(3)))./2.0,(1.0./a.^2.*(gamma-1.0).*(a.*lambda(2)-a.*lambda(3)-lambda(1).*u.*2.0+lambda(2).*u+lambda(3).*u))./2.0,lambda(2).*(gamma-1.0).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+lambda(3).*(gamma-1.0).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)-(1.0./a.^2.*lambda(1).*u.^2.*(gamma-1.0))./2.0],[3,3]);
            
            % vetor de fluxo em função dos autovalores
            F = @(q, lambda) A(q(2)/q(1), sqrt(gamma * (gamma - 1) * (q(3) ./ q(1) - q(2).^2 ./ (2 * q(1).^2))), lambda) * q';

            % determina o passo no tempo
            if pressureRatio == 50 || pressureRatio == 100
                CFL = 0.6;
            else
                CFL = 0.3;
            end

            eigenvalues = [u0, u0 + a0, u0 - a0];
            dt =  CFL * dx / max(eigenvalues(:));

            % coeficiente do esquema numérico
            lamb = dt / dx;

            % parâmetros que determinam o método de marcha
            % no tempo
            theta = 1;
            xi = 0;
            
            % variação inicial das propriedades
            dq = zeros(3 * M - 6, 1);

            while 1

                % efetua um incremento temporal
                t = t + dt;

                if t * (L_ref / u_ref) > 1 
                    break
                end

                % mostra a solução corrente
                %show(x, q, gamma);
               
                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                % calcula os três autovalores + e - em todos os 
                % pontos da malha
                u = q(:, 2) ./ q(:, 1);
                a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2)));
                
                lp = [u, u + a, u - a];
                lm = lp;

                lp = max(lp, 0);
                lm = min(lm, 0);

                % aloca as matrizes e vetores para os sistemas da versão
                % implícita do esquema de Steger-Warming
                B_1 = zeros(3 * M - 6, 3 * M - 6);
                B_2 = zeros(3 * M - 6, 3 * M - 6);
                f = zeros(3 * M - 6, 1);
                i = 1;

                for j = 3:M
                    % monta a matriz de coeficientes para o primeiro
                    % sistema de equações lineares
                    if j == 3
                        B_1(1:3, 1:3) = eye(3) + theta * lamb * A(u(j), a(j), lp(j,:)) / (1 + xi);
                    else
                        B_1(3 * (i - 1) + 1:3*i, (1:3) + 3 * (i - 2)) = - theta * lamb * A(u(j - 1), a(j - 1), lp(j - 1,:))  / (1 + xi);
                        B_1(3 * (i - 1) + 1:3*i, (1:3) + 3 * (i - 1)) = eye(3) + theta * lamb * A(u(j), a(j), lp(j,:))  / (1 + xi);
                    end

                    % monta a matriz de coeficientes para o segundo
                    % sistema de equações lineares
                    if j == M
                        B_2(3*M - 8:3*M - 6, 3*M - 8:3*M - 6) = eye(3) - theta * lamb * A(u(j), a(j), lm(j,:)) / (1 + xi);
                    else
                        B_2(3 * (i - 1) + 1:3*i, (1:3) + 3 * (i - 1)) = eye(3) - theta * lamb * A(u(j), a(j), lm(j,:)) / (1 + xi);
                        B_2(3 * (i - 1) + 1:3*i, (1:3) + 3 * i) = theta * lamb * A(u(j + 1), a(j + 1),  lm(j + 1,:)) / (1 + xi);
                    end
            
                    % monta o vetor do lado direito para o primeiro sistema
                    % de equações lineares
                    f(3 * (i - 1) + 1:3*i) =  3 * (F(q(j, :), lp(j, :)) - F(q(j, :), lm(j,:))) ...
                                           - 4 * (F(q(j - 1, :), lp(j - 1, :)) - F(q(j + 1, :), lm(j + 1, :))) ...
                                           + 1 * (F(q(j - 2, :), lp(j - 2, :)) - F(q(j + 2, :), lm(j + 2, :)));

                    i = i + 1;
                end

                % monta o vetor do lado direito para o primeiro sistema
                % de equações lineares
                f = - lamb * f / (2 * (xi + 1));

                % resolve o primeiro sistema
                dq_star = linsolve(B_1, f + xi * dq / (1 + xi));

                % resolve o segundo sistema
                dq = linsolve(B_2, dq_star);

                % atualiza a solução corrente
                for j = 3:M
                    indexes = 3 * (j - 3) + 1:3 * (j - 2);
                    next_q(j, :) = q(j, :) + dq(indexes(1):indexes(3))';
                end

                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                % aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
    
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
            end

        case EXP_STEGER_WARMING

            % matriz jacobiana de fluxo em função dos autovalores
            A = @(u, a, lambda) reshape([-lambda(1).*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0)-(1.0./a.^2.*lambda(2).*(a.*u-(u.^2.*(gamma-1.0))./2.0))./2.0+(1.0./a.^2.*lambda(3).*(a.*u+(u.^2.*(gamma-1.0))./2.0))./2.0,-lambda(1).*u.*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0)-lambda(2).*((1.0./a.^2.*u)./2.0+1.0./(a.*2.0)).*(a.*u-(u.^2.*(gamma-1.0))./2.0)+lambda(3).*((1.0./a.^2.*u)./2.0-1.0./(a.*2.0)).*(a.*u+(u.^2.*(gamma-1.0))./2.0),lambda(1).*u.^2.*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0).*(-1.0./2.0)-lambda(2).*(a.*u-(u.^2.*(gamma-1.0))./2.0).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+lambda(3).*(a.*u+(u.^2.*(gamma-1.0))./2.0).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0),(1.0./a.^2.*lambda(2).*(a-u.*(gamma-1.0)))./2.0-(1.0./a.^2.*lambda(3).*(a+u.*(gamma-1.0)))./2.0+1.0./a.^2.*lambda(1).*u.*(gamma-1.0),1.0./a.^2.*lambda(1).*u.^2.*(gamma-1.0)+(1.0./a.^2.*lambda(2).*(a+u).*(a+u-gamma.*u))./2.0+(1.0./a.^2.*lambda(3).*(a-u).*(a-u+gamma.*u))./2.0,lambda(2).*(a-u.*(gamma-1.0)).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)-lambda(3).*(a+u.*(gamma-1.0)).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+(1.0./a.^2.*lambda(1).*u.^3.*(gamma-1.0))./2.0,(1.0./a.^2.*(gamma-1.0).*(lambda(1).*-2.0+lambda(2)+lambda(3)))./2.0,(1.0./a.^2.*(gamma-1.0).*(a.*lambda(2)-a.*lambda(3)-lambda(1).*u.*2.0+lambda(2).*u+lambda(3).*u))./2.0,lambda(2).*(gamma-1.0).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+lambda(3).*(gamma-1.0).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)-(1.0./a.^2.*lambda(1).*u.^2.*(gamma-1.0))./2.0],[3,3]);
            
            % vetor de fluxo em função dos autovalores
            F = @(q, lambda) A(q(2)/q(1), sqrt(gamma * (gamma - 1) * (q(3) ./ q(1) - q(2).^2 ./ (2 * q(1).^2))), lambda) * q';

            % determina o passo no tempo
            if order == 1
                CFL = 0.3;
            elseif order == 2
                CFL = 0.1;
            end
            eigenvalues = [u0, u0 + a0, u0 - a0];
            dt = CFL * dx / max(eigenvalues(:));

            % coeficiente do esquema numérico
            lamb = dt / dx;

            while 1

                % efetua um incremento temporal
                t = t + dt;
                
                if t * (L_ref / u_ref) > 1 
                    break
                end
                
                % mostra a solução corrente
                %show(x, q, gamma);
                
                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);
                
                % calcula os três autovalores + e - em todos os 
                % pontos da malha
                u = q(:, 2) ./ q(:, 1);
                a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 
                
                lp = [u, u + a, u - a];
                lm = lp;
                
                lp = max(lp, 0);
                lm = min(lm, 0);

                % calcula a solução no próximo instante de tempo
                for j = 3:M
                    % verifica em qual ordem o esquema deve ser aplicado
                    if order == 1
                        next_q(j, :) = q(j, :) - lamb  * (F(q(j, :), lp(j, :)) - F(q(j - 1, :), lp(j - 1, :)))' ...
                                               - lamb  * (F(q(j + 1, :), lm(j + 1, :)) - F(q(j, :), lm(j, :)))';
                    elseif order == 2
                        next_q_bar_jm = q(j - 1, :) ...
                                   - lamb * (F(q(j - 1, :), lp(j - 1, :)) - F(q(j - 2, :), lp(j - 2, :)))' ...
                                   - lamb * (F(q(j, :), lm(j, :)) - F(q(j - 1, :), lm(j - 1, :)))';
                    
                        next_q_bar_j = q(j, :) ...
                                  - lamb * (F(q(j, :), lp(j, :)) - F(q(j - 1, :), lp(j - 1, :)))' ...
                                  - lamb * (F(q(j + 1, :), lm(j + 1, :)) - F(q(j, :), lm(j, :)))';
                    
                        next_q_bar_jp = q(j + 1, :) ...
                                  - lamb * (F(q(j + 1, :), lp(j + 1, :)) - F(q(j, :), lp(j, :)))' ...
                                  - lamb * (F(q(j + 2, :), lm(j + 2, :)) - F(q(j + 1, :), lm(j + 1, :)))';

                        
                        next_q(j, :) = (q(j, :) + next_q_bar_j) / 2 ...
                                     - lamb * (F(q(j, :), lp(j, :)) - 2 * F(q(j - 1, :), lp(j - 1, :)) + F(q(j - 2, :), lp(j - 2, :)) + F(next_q_bar_j, lp(j, :)) - F(next_q_bar_jm, lp(j - 1, :)))' / 2 ...
                                     + lamb * (F(q(j + 2, :), lm(j + 2, :)) - 2 * F(q(j + 1, :), lm(j + 1, :)) + F(q(j, :), lm(j, :)) - (F(next_q_bar_jp, lm(j + 1, :)) - F(next_q_bar_j, lm(j, :))))' / 2;
                    end
                end

                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                %aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
        
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
             end


        case {LAX_WENDROFF, MAC_COMARCK}

            % matriz jacobiana de fluxo em função das variáveis conservadas
            A = @(q) [0, 1, 0; 
                      (gamma - 3) * q(2)^2/(2 * q(1)^2), (3 - gamma) * q(2)/q(1), gamma - 1;
                      -gamma * q(2) * q(3)/q(1)^2+(gamma-1) * q(2)^3/q(1)^3, gamma * q(3)/q(1) - 3 * (gamma - 1) * q(2)^2/(2*q(1)^2), gamma * q(2)/q(1)];

            % vetor de fluxo em função das variáveis conservadas
            F = @(q) [q(2); ...
                (3 - gamma) * q(2)^2/(2 * q(1)) + (gamma - 1) * q(3); ...
                gamma * q(2) * q(3)/q(1) - (gamma - 1) * q(2)^3/(2 * q(1)^2)];


            % determina o passo no tempo
            CFL = 0.2;
            dt = CFL * dx / max(abs(u0) + a0);

            % coeficiente dos esquemas
            lamb = dt / dx;

            % parâmetros da dissipação artificial linear (segunda e quarta)
            eps_e_4 = 0;
            eps_e_2 = 1 / 15;

            while 1

                % efetua um incremento temporal
                t = t + dt;

                if t * (L_ref / u_ref) > 1 
                    break
                end

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);
    
                 % verifica qual é o método explícito que deve ser aplicado
                if METHOD == LAX_WENDROFF

                    % aplica o esquema de LW em cada ponto interior da malha
                    for j = 3:M
                        A_half_plus = A((q(j + 1, :) + q(j, :))/2);
                        A_half_minus = A((q(j, :) + q(j - 1, :))/2);

                        next_q(j, :) = q(j, :) - lamb * (F(q(j + 1, :)) - F(q(j - 1,:)))' / 2 + ...
                                        lamb^2 * (A_half_plus * (F(q(j + 1,:)) - F(q(j, :))) - A_half_minus * (F(q(j, :)) - F(q(j - 1, :))))'/2 - eps_e_4 * (q(j - 2, :) - 4 * q(j - 1, :) + 6 * q(j, :) - 4 * q(j + 1, :) + q(j + 2, :)) + eps_e_2 * (q(j + 1, :) - 2 * q(j, :) + q(j - 1, :));
                    end
    
                elseif METHOD == MAC_COMARCK

                    % aplica o esquema de MC em cada ponto interior da
                    % malha
                    for j = 3:M
                        % preditor
                        q_bar_j = q(j, :) - lamb * (F(q(j + 1, :)) - F(q(j, :)))';
                        
                        % corretor
                        q_bar_j_minus = q(j - 1, :) - lamb * (F(q(j, :)) - F(q(j - 1, :)))';
                        next_q(j, :) = (q(j, :) + q_bar_j) / 2 - lamb * (F(q_bar_j) - F(q_bar_j_minus))' / 2 - eps_e_4 * (q(j - 2, :) - 4 * q(j - 1, :) + 6 * q(j, :) - 4 * q(j + 1, :) + q(j + 2, :)) + eps_e_2 * (q(j + 1, :) - 2 * q(j, :) + q(j - 1, :));
                    end
                end

                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                % aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
        
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
            end

        case {EXP_BEAM_WARMING, IMP_BEAM_WARMING} 

            % matriz jacobiana de fluxo em função das variáveis conservadas
            A = @(q) [0, 1, 0; 
                      (gamma - 3) * q(2)^2/(2 * q(1)^2), (3 - gamma) * q(2)/q(1), gamma - 1;
                      -gamma * q(2) * q(3)/q(1)^2+(gamma-1) * q(2)^3/q(1)^3, gamma * q(3)/q(1) - 3 * (gamma - 1) * q(2)^2/(2*q(1)^2), gamma * q(2)/q(1)];

            % vetor de fluxo em função das variáveis conservadas
            F = @(q) [q(2); ...
                (3 - gamma) * q(2)^2/(2 * q(1)) + (gamma - 1) * q(3); ...
                gamma * q(2) * q(3)/q(1) - (gamma - 1) * q(2)^3/(2 * q(1)^2)];


            if METHOD == IMP_BEAM_WARMING

                % estipula um valor de CFL
                CFL = 0.088;

                % parâmetro de dissipação explícita (quarta ordem)
                if(pressureRatio < 100)
                     eps_e = 1 / 20;
                else
                     eps_e = 1 / 100;
                end
               
                % parâmetro de dissipação implícita (segunda ordem)
                eps_i = 3  * eps_e;

            elseif METHOD == EXP_BEAM_WARMING

                % estipula um valor de CFL
                CFL = 0.2;

                % parâmetros da dissipação artificial linear (segunda e
                % quarta)
                eps_e_2 = 1 / 5;
                eps_e_4 = 0;
            end

            % determina o passo no tempo
            dt = CFL * dx / max(abs(u0) + a0);
        
            % coeficiente do sistema tridiagonal em bloco
            lamb = dt / dx;
            
            % laço para efetuar a marcha no tempo
            while 1

                % efetua um incremento temporal
                t = t + dt;

                if t * (L_ref / u_ref) > tEnd
                    break
                end
    
                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);
    
                % verifica se o tratamento deve ser implícito ou explícito
                if METHOD == IMP_BEAM_WARMING
    
                    % monta a matriz de coeficiente do sistema tridiagonal
                    % em blocos, assim como o vetor do lado direito
                    B = zeros(3 * M - 6, 3 * M - 6);
                    f = zeros(3 * M - 6, 1);
                    g = zeros(3 * M - 6, 1);
                    k = 1;
                    i = 1;
                    for j = 3:M
                        if j == 3
                            B(1:3, 1:3) = (1 + 2 * eps_i) * eye(3); 
                            B(1:3, 4:6) = (lamb * A(q(j + 1, :)) / 4 - eps_i * eye(3));
                        elseif j == M
                            B(3*M - 8:3*M - 6, (3*M - 11):(3*M-9)) = (-lamb * A(q(j - 1, :))/4 - eps_i * eye(3));
                            B(3*M - 8:3*M - 6, 3*M - 8:3*M - 6) = (1+ 2 * eps_i) * eye(3);
                        else
                            B(3 * (i - 1) + 1:3*i, (1:3) + 3 * (k - 1)) = (-lamb * A(q(j - 1,:))/4 - eps_i * eye(3));
                            B(3 * (i - 1) + 1:3*i, (4:6) + 3 * (k - 1)) = (1 + 2 * eps_i) *  eye(3);
                            B(3 * (i - 1) + 1:3*i, (7:9) + 3 * (k - 1)) = (lamb * A(q(j + 1, :))/4 - eps_i * eye(3));
                            k = k + 1;
                        end
    
                        f(3 * (i - 1) + 1:3*i) =  (F(q(j + 1, :)) - F(q(j - 1, :)));
                        g(3 * (i - 1) + 1:3*i) = (q(j - 2, :) - 4 * q(j - 1, :) + 6 * q(j, :) - 4 * q(j + 1, :) + q(j + 2, :))';
                        
                        i = i + 1;
                    end
    
                    % coloca o fator multiplicativo nos vetores do lado
                    % direito do sistema tridiagonal em blocos
                    f = - lamb * f/2;
                    g = - eps_e * g;
            
                    % resolve o sistema tridiagonal em blocos para encontrar
                    % a correção que leva a solução no instante corrente para 
                    % o próximo instante
                    dq = linsolve(B, f + g);

                    for j = 3:M
                        indexes = 3 * (j - 3) + 1:3 * (j - 2);
                        next_q(j, :) = q(j, :) + dq(indexes(1):indexes(3))';
                    end
    
                elseif METHOD == EXP_BEAM_WARMING
                    for j = 3:M
                        next_q(j, :) = q(j, :) - lamb * (F(q(j + 1, :)) - F(q(j - 1, :)))' / 2 + ... 
                                    - eps_e_4 * (q(j - 2, :) - 4 * q(j - 1, :) + 6 * q(j, :) - 4 * q(j + 1, :) + q(j + 2, :)) + eps_e_2 * (q(j - 1, :)  - 2 * q(j, :) + q(j + 1, :));
                    end   
                end
    
                % aplica as condições de contorno à esquerda
                next_q(1, 1) = rho_l; % densidade
                next_q(1, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(1, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                next_q(2, 1) = rho_l; % densidade
                next_q(2, 2) = rho_l * u_l; % quantidade de movimento na direção x
                next_q(2, 3) = p_l / (gamma - 1) + rho_l * u_l^2 / 2; % energia total
        
                % aplica as condições de contorno à direita
                next_q(M + 1, 1) = rho_r; % densidade
                next_q(M + 1, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 1, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
            
                next_q(M + 2, 1) = rho_r; % densidade
                next_q(M + 2, 2) = rho_r * u_r; % quantidade de movimento na direção x
                next_q(M + 2, 3) = p_r / (gamma - 1) + rho_r * u_r^2 / 2; % energia total
    
                % a solução no próximo instante se torna a solução corrente
                q = next_q;
            end
    end

    % obtém a solução em termos das variáveis primitivas (t = 1)
    rho = q(:, 1) * rho_ref;
    u = q(:, 2) * u_ref ./ q(:, 1);
    E = q(:, 3) * u_ref^2 ./ q(:, 1);
    p = (gamma - 1) * (q(:, 3) - q(:, 2).^2 ./ (2 * q(:, 1))) * p_ref;

    totalTime = toc(startTime);
end