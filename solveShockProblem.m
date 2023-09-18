function [x, rho, u, E, p] = solveShockProblem(METHOD, rho_l, u_l, p_l, rho_r, u_r, p_r)
    %% Definição das variáveis

    % variáveis que identificam o esquema numérico
    LAX_WENDROFF = 1;
    MAC_COMARCK = 2;
    EXP_BEAM_WARMING = 3;
    IMP_BEAM_WARMING = 4;
    EXP_STEGER_WARMING = 5;
    IMP_STEGER_WARMING = 8;
    AUSM_PLUS = 6;
    VAN_LEER = 7;
    AUSM_PLUS_SEC = 9;

    % número de pontos interiores da malha
    M = 1998;
    
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

    % vetor de fluxo em função das variáveis conservadas
    F = @(q) [q(2); ...
                (3 - gamma) * q(2)^2/(2 * q(1)) + (gamma - 1) * q(3); ...
                gamma * q(2) * q(3)/q(1) - (gamma - 1) * q(2)^3/(2 * q(1)^2)];


    % matriz jacobiana de fluxo em função das variáveis conservadas
    % A = @(q) [0, 1, 0; 
    %           (gamma - 3) * q(2)^2/(2 * q(1)^2), (3 - gamma) * q(2)/q(1), gamma - 1;
    %           -gamma * q(2) * q(3)/q(1)^2+(gamma-1) * q(2)^3/q(1)^3, gamma * q(3)/q(1) - 3 * (gamma - 1) * q(2)^2/(2*q(1)^2), gamma * q(2)/q(1)];

    % 0 <= u <= c
    Fpp = @(q, a) (q(1)/(2*gamma))*[2 * (gamma - 1) * q(2)/q(1) + q(2)/q(1) + a;
                    2 * (gamma - 1) * (q(2)/q(1))^2 + (q(2)/q(1)+a)^2;
                    (gamma - 1) * (q(2)/q(1))^3 + (q(2)/q(1) + a)^3/2 + (3 - gamma) * (q(2)/q(1) + a) * a^2/(2*(gamma - 1))];

    Fnp = @(q, a) (q(1)/(2*gamma)) * [q(2)/q(1) - a; (q(2)/q(1) - a)^2; (q(2)/q(1) - a)^3/2 + (3-gamma) * (q(2)/q(1) - a)*a^2/(2*(gamma - 1))]; 

    % -c <= u <= 0
    Fpn = @(q, a) (q(1)/(2*gamma)) * [q(2)/q(1) + a; (q(2)/q(1) + a)^2; (q(2)/q(1) + a)^3/2 + (3 - gamma) * (q(2)/q(1) + a)*a^2/(2 *(gamma - 1))];
    
    Fnn = @(q, a) (q(1)/(2*gamma))*[2 * (gamma - 1) * q(2)/q(1) + q(2)/q(1) - a;
                    2 * (gamma - 1) * (q(2)/q(1))^2 + (q(2)/q(1)-a)^2;
                    (gamma - 1) * (q(2)/q(1))^3 + (q(2)/q(1) - a)^3/2 + (3 - gamma) * (q(2)/q(1) - a) * a^2/(2*(gamma - 1))];

    %H = @(p, rho, u) p * gamma / ((gamma - 1) * rho) + u^2/2;
    
    %H = @(a, u) u^2/2 + a^2/(gamma - 1);
    %phi = @(rho, u, Hv) [rho; rho * u; rho * Hv];
%fluxLimiter = @(r) max([0, min(2*r, 1), min(r, 2)]);
 %   fluxLimiter = @(r) (r + abs(r))/(1 + abs(r));
%   fluxLimiter = @(r) max(0, min(1, r)); %mel
   % fluxLimiter = @(r) r * (3 * r + 1) / (r + 1)^2 .* (r > 0) ;
fluxLimiter = @(r) 2 * (r + abs(r)) / (r + 3); % bom pra 20
 % fluxLimiter = @(r) abs((r^2+ r)/(r^2 + 1));
%fluxLimiter = @(r) max(0, min([2*r, 0.5*(1+r), 2]));
%fluxLimiter = @(r)max(0, min(2 * r, min((1 + 2 * r)/3, 2)));
H = @(q, a) (q(2)/q(1))^2/2 + a^2 / (gamma - 1);
    phi = @(q,a) [q(1); q(2); q(1) * H(q, a)];

    f_mass_plus = @(rho, a, M) rho * a * ((M + 1)/2)^2;
    f_mass_minus = @(rho, a, M) -rho * a * ((M - 1)/2)^2;

    F_plus_vanLeer = @(q, a) f_mass_plus(q(1), a, q(2) / (q(1) * a)) * ...
                      [1; 
                      ((gamma - 1) * q(2) / q(1) + 2 * a) / gamma; 
                      ((gamma - 1) * q(2) / q(1) + 2 * a)^2 / (2 * (gamma^2 - 1))];

    F_minus_vanLeer = @(q, a) f_mass_minus(q(1), a, q(2) / (q(1) * a)) * ...
                              [1;
                               ((gamma - 1) * q(2) / q(1) - 2 * a) / gamma;
                               ((gamma - 1) * q(2) / q(1) - 2 * a)^2 / (2 * (gamma^2 - 1))];
                        

    % F_plus = @(rho, a, u, M) [f_mass_positive(rho, a, M); 
    %             f_mass_positive(rho, a, M) * ((gamma - 1) * u + 2 * a)/gamma;
    %             f_mass_positive(rho, a, M) * ((gamma - 1) * u + 2 * a)^2/(2*(gamma^2 - 1))];
    % 
    % F_minus = @(rho, a, u, M) [f_mass_negative(rho, a, M); 
    %             f_mass_negative(rho, a, M) * ((gamma - 1) * u - 2 * a)/gamma;
    %             f_mass_negative(rho, a, M) * ((gamma - 1) * u - 2 * a)^2/(2*(gamma^2 - 1))];
    
    %phi = @(rho, u, Hv, a) [rho*a; rho * u*a; rho * Hv*a];


    % X = @(rho, u, a) [ 1, 1 / (2 * a^2), 1 / (2 * a^2); 0, 1 / (2 * rho * a), - 1 / (2 * rho * a); 0, 1/2, 1/2];
    % invX = @(rho, u , a) [1, 0, -1/a^2; 0, rho * a, 1; 0, - rho * a, 1];
    % 
    % Mm = @(rho, u) [1 0 0; u rho 0; u^2/2, rho * u, 1/(gamma - 1)];
    % invMm = @(rho, u) [1, 0 ,0; -u/rho, 1/rho, 0; (gamma - 1) * u^2/2, - (gamma - 1) * u, gamma - 1];
    % 
    % invT = @(rho, u , a) invX(rho, u, a) * invMm(rho, u);
    % T = @(rho, u , a) Mm(rho, u) * X(rho, u, a);



    A = @(u,a,lambda) reshape([-lambda(1).*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0)-(1.0./a.^2.*lambda(2).*(a.*u-(u.^2.*(gamma-1.0))./2.0))./2.0+(1.0./a.^2.*lambda(3).*(a.*u+(u.^2.*(gamma-1.0))./2.0))./2.0,-lambda(1).*u.*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0)-lambda(2).*((1.0./a.^2.*u)./2.0+1.0./(a.*2.0)).*(a.*u-(u.^2.*(gamma-1.0))./2.0)+lambda(3).*((1.0./a.^2.*u)./2.0-1.0./(a.*2.0)).*(a.*u+(u.^2.*(gamma-1.0))./2.0),lambda(1).*u.^2.*((1.0./a.^2.*u.^2.*(gamma-1.0))./2.0-1.0).*(-1.0./2.0)-lambda(2).*(a.*u-(u.^2.*(gamma-1.0))./2.0).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+lambda(3).*(a.*u+(u.^2.*(gamma-1.0))./2.0).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0),(1.0./a.^2.*lambda(2).*(a-u.*(gamma-1.0)))./2.0-(1.0./a.^2.*lambda(3).*(a+u.*(gamma-1.0)))./2.0+1.0./a.^2.*lambda(1).*u.*(gamma-1.0),1.0./a.^2.*lambda(1).*u.^2.*(gamma-1.0)+(1.0./a.^2.*lambda(2).*(a+u).*(a+u-gamma.*u))./2.0+(1.0./a.^2.*lambda(3).*(a-u).*(a-u+gamma.*u))./2.0,lambda(2).*(a-u.*(gamma-1.0)).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)-lambda(3).*(a+u.*(gamma-1.0)).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+(1.0./a.^2.*lambda(1).*u.^3.*(gamma-1.0))./2.0,(1.0./a.^2.*(gamma-1.0).*(lambda(1).*-2.0+lambda(2)+lambda(3)))./2.0,(1.0./a.^2.*(gamma-1.0).*(a.*lambda(2)-a.*lambda(3)-lambda(1).*u.*2.0+lambda(2).*u+lambda(3).*u))./2.0,lambda(2).*(gamma-1.0).*(1.0./(gamma.*2.0-2.0)+u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)+lambda(3).*(gamma-1.0).*(1.0./(gamma.*2.0-2.0)-u./(a.*2.0)+(1.0./a.^2.*u.^2)./4.0)-(1.0./a.^2.*lambda(1).*u.^2.*(gamma-1.0))./2.0],[3,3]);
    F = @(q, lambda) A(q(2)/q(1),sqrt(gamma * (gamma - 1) * (q(3) ./ q(1) - q(2).^2 ./ (2 * q(1).^2))),lambda) * q';




    % 
    % T = @(u, a) [1, 1 / (2 * a^2), 1 / (2 * a^2);
    %              u, (u + a)/(2 * a^2), (u - a)/(2 * a^2);
    %              u^2 / 2, (u^2/2 + u * a + a^2/(gamma - 1))/(2 * a^2), (u^2/2 - u * a + a^2/(gamma - 1))/(2 * a^2)];
    % 
    % invT = @(u, a) [(a^2 - (gamma - 1) * u^2/2)/a^2, (gamma - 1) * u/a^2, -(gamma - 1)/a^2;
    %                  (gamma - 1) * u^2/2 - u * a, (u + a) - gamma * u, gamma-1;
    %                  (gamma - 1)*u^2/2 + u * a, (u - a) - gamma * u, gamma - 1];
    % 
    % absLamb = @(u, a) [abs(u), 0 ,0; 0, abs(u + a), 0; 0, 0, abs(u - a)];
    % 
    % absA = @(u, a) T(u, a) * absLamb(u, a) * invT(u, a);

    % L_p = @(u, a) [(u + abs(u))/2, 0 ,0; 0, (u + a + abs(u + a)) / 2, 0; 0, 0, (u - a + abs(u - a))/2]; 
    % L_m = @(u, a) [(u - abs(u))/2, 0 ,0; 0, (u + a - abs(u + a)) / 2, 0; 0, 0, (u - a - abs(u - a))/2]; 

   

    %   A_p = @(u, a) T( u, a) * L_p(u, a) * invT( u, a);
    % A_n = @(u, a) T( u, a) * L_m(u, a) * invT( u, a);

    % A_p = @(rho, u, a) T(rho, u, a) * Lp* invT(rho, u, a);
    % A_n = @(rho, u, a) T(rho, u, a) * Lm* invT(rho, u, a);
    % 
    % Fplus = @(q, a ) A_p(q(1), q(2)/q(1), a) * q'; 
    % Fminus = @(q, a, Lm) A_n(q(1), q(2)/q(1), a) * q'; 

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

        case VAN_LEER

            CFL = 0.1;%0.0001;
            dt = CFL * dx /max(abs(u0) + a0);

                 % coeficiente dos esquemas
            lamb = dt / dx;
            while 1

                % efetua um incremento temporal
                t = t + dt;
   
                if t * (L_ref / u_ref) > 1 
                    break
                end

                  subplot(2,2, 1);
                plot(x, q(:, 1), 'k');

                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
                ax = gca;
                ax.FontSize = 16; 
                subplot(2,2, 2);
                plot(x, q(:, 2) ./ q(:, 1), 'r');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$u(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 3);
                plot(x, (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1))), 'b');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$p(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 4);
                plot(x, q(:, 3)./q(:, 1), 'g');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$E(x, t)$', 'Interpreter', 'latex');
                pause(0.0001);


                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                rho = q(:, 1);
                u = q(:, 2) ./ rho;
             a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 
               p = (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1)));
                 
           %     a = sqrt(gamma * p ./ rho);
                for j = 3:M 
                    Mj = u(j)/a(j);
                    if Mj >= 1
                        Fj_positive = F(q(j, :));
                        Fjm_positive = F(q(j - 1, :));
                        Fjmm_positive = F(q(j - 2, :));

                        Fj_negative = 0;
                        Fjp_negative = 0;
                        Fjpp_negative = 0;
                    elseif Mj <= -1

                        Fj_positive = 0;
                        Fjm_positive = 0;
                        Fjmm_positive = 0;

                        Fj_negative = F(q(j, :));
                        Fjp_negative = F(q(j + 1, :));
                        Fjpp_negative = F(q(j + 2, :));
                    else

                        Fj_positive = F_plus_vanLeer(q(j, :), a(j));
                        Fjm_positive = F_plus_vanLeer(q(j - 1, :), a(j - 1));
                        Fjmm_positive = F_plus_vanLeer(q(j - 2, :), a(j - 2));

                        Fj_negative = F_minus_vanLeer(q(j, :), a(j));
                        Fjp_negative = F_minus_vanLeer(q(j + 1, :), a(j + 1));
                        Fjpp_negative = F_minus_vanLeer(q(j + 2, :), a(j + 2));
                    end

                   Z_plus_1 = [fluxLimiter((rho(j + 1) - rho(j)) / (rho(j) - rho(j - 1)));...
                                     fluxLimiter((u(j + 1) - u(j)) / (u(j) - u(j - 1)));...
                                     fluxLimiter((p(j + 1) - p(j)) / (p(j) - p(j - 1)))];
                   Z_minus_1 = [fluxLimiter((rho(j + 1) - rho(j)) / (rho(j + 2) - rho(j + 1)));...
                                 fluxLimiter((u(j + 1) - u(j)) / (u(j + 2) - u(j + 1)));...
                                 fluxLimiter((p(j + 1) - p(j)) / (p(j + 2) - p(j + 1)))];
                    
                    Z_plus_1(isnan(Z_plus_1)) = 1;
                    Z_minus_1(isnan(Z_minus_1)) = 1;
                    % 
                    %  Z_plus_1(:) = 1;
                    % Z_minus_1(:) = 1;
                    

                   Z_plus_2 = [fluxLimiter((rho(j) - rho(j - 1)) / (rho(j - 1) - rho(j - 2)));...
                                     fluxLimiter((u(j) - u(j - 1)) / (u(j - 1) - u(j - 2)));...
                                     fluxLimiter((p(j) - p(j - 1)) / (p(j - 1) - p(j - 2)))];

                   Z_minus_2 = [fluxLimiter((rho(j) - rho(j - 1)) / (rho(j + 1) - rho(j)));...
                                 fluxLimiter((u(j) - u(j - 1)) / (u(j + 1) - u(j)));...
                                 fluxLimiter((p(j) - p(j - 1)) / (p(j + 1) - p(j)))];

                    Z_plus_2(isnan(Z_plus_2)) = 1;
                    Z_minus_2(isnan(Z_minus_2)) = 1;

                    %  Z_plus_2(:) = 1;
                   %  Z_minus_2(:) = 1;
                   

                    F_bar_jphalf_positive = Fj_positive + Z_plus_1(1).*(Fj_positive - Fjm_positive)/2;%Fj_positive;%Fj_positive + (Fj_positive - Fjm_positive)/2;
                    
                    
                    F_bar_jmhalf_positive = Fjm_positive + Z_plus_2(1).*(Fjm_positive - Fjmm_positive)/2;%Fjm_positive;%Fjm_positive + (Fjm_positive - Fjmm_positive)/2;
    
                    F_bar_jphalf_negative =Fjp_negative - Z_minus_1(1).*(Fjpp_negative - Fjp_negative)/2; %Fjp_negative;%Fjp_negative - (Fjpp_negative - Fjp_negative)/2;
                    
                    F_bar_jmhalf_negative = Fj_negative - Z_minus_2(1).*(Fjp_negative - Fj_negative)/2;%Fj_negative;%Fj_negative - (Fjp_negative - Fj_negative)/2;
    
                    next_q(j, :) = q(j, :) - lamb * (F_bar_jphalf_positive - F_bar_jmhalf_positive)'  ...
                            - lamb * (F_bar_jphalf_negative - F_bar_jmhalf_negative)';
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

            % 
            %     rho = q(:, 1);
            %     u = q(:, 2) ./ rho;
            %     a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2)));
            % 
            %        CFL = 0.05;
            % 
            % dt = CFL * dx /max(abs(u0) + a0);
            % 
            %      % coeficiente dos esquemas
            % lamb = dt / dx;
            end


        case AUSM_PLUS_SEC

            CFL = 0.05;

            dt = CFL * dx /max(abs(u0) + a0);

            alpha = 3 / 16;
            beta = 1/8;

                 % coeficiente dos esquemas
            lamb = dt / dx;

            while 1

                % efetua um incremento temporal
                t = t + dt;

                if t * (L_ref / u_ref) > 1 
                    break
                end


                subplot(2,2, 1);
                plot(x, q(:, 1), 'k');

                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
                ax = gca;
                ax.FontSize = 16; 
                subplot(2,2, 2);
                plot(x, q(:, 2) ./ q(:, 1), 'r');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$u(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 3);
                plot(x, (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1))), 'b');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$p(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 4);
                plot(x, q(:, 3)./q(:, 1), 'g');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$E(x, t)$', 'Interpreter', 'latex');
                pause(0.0001);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                rho = q(:, 1);
                u = q(:, 2) ./ rho;
               % a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 
                p = (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1)));
                a = sqrt(gamma * p ./ rho);
                V = [rho, u, p];

                for j = 3:M

                   Z_plus = [fluxLimiter((rho(j + 1) - rho(j)) / (rho(j) - rho(j - 1))),...
                                     fluxLimiter((u(j + 1) - u(j)) / (u(j) - u(j - 1))),...
                                     fluxLimiter((p(j + 1) - p(j)) / (p(j) - p(j - 1)))];
                   Z_minus = [fluxLimiter((rho(j + 1) - rho(j)) / (rho(j + 2) - rho(j + 1))),...
                                 fluxLimiter((u(j + 1) - u(j)) / (u(j + 2) - u(j + 1))),...
                                 fluxLimiter((p(j + 1) - p(j)) / (p(j + 2) - p(j + 1)))];
                    
                    Z_plus(isnan(Z_plus)) = 1;
                    Z_minus(isnan(Z_minus)) = 1;
                    % 
                    % Z_plus(:) = 1;
                    % Z_minus(:) = 1;

                    VL = V(j, :) +  Z_plus.*(V(j, :)- V(j - 1, :))/2;

                    rhoL = VL(1);
                    uL = VL(2);
                    pL = VL(3);
                    EL = pL / ((gamma - 1) * rhoL) + uL^2/2;

                    VR = V(j + 1, :) - Z_minus.*(V(j + 2, :) - V(j + 1, :))/2;

                    rhoR = VR(1);
                    uR = VR(2);
                    pR = VR(3);
                    ER = pR / ((gamma - 1) * rhoR) + uR^2/2;

                    qL = [rhoL; rhoL * uL; rhoL * EL];
                    qR = [rhoR; rhoR * uR; rhoR * ER];

                    aL = sqrt(gamma * pL / rhoL);
                    aR = sqrt(gamma * pR / rhoR);
                    
                    aStarL = sqrt(2 * (gamma - 1) * H(qL, aL) / (gamma + 1)); % certo
                    
                    aTildeL = aStarL * min(1, aStarL / abs(uL)); % certo

                    aStarR = sqrt(2 * (gamma - 1)  * H(qR, aR) / (gamma + 1)); % certo

                    aTildeR = aStarR * min(1, aStarR / abs(uR)); % certo

                    aHalfjp = min(aTildeL, aTildeR); %certo

                    Mj = uL / aHalfjp; % certo
                    Mjp = uR / aHalfjp; % certo

                    if abs(Mj) >= 1
                        Mj_plus = (Mj + abs(Mj))/2; % certo
                        pj_plus = (1 + sign(Mj))/2; % certo
                    else
                        Mj_plus = (Mj + 1)^2/4 + beta * (Mj^2 - 1)^2; % certo
                        pj_plus = (Mj + 1)^2 * (2 - Mj)/4 + alpha * Mj * (Mj^2-1)^2; % certo
                    end

                    if abs(Mjp) >= 1
                        Mjp_negative = (Mjp - abs(Mjp)) / 2; % certo
                        pjp_negative = (1 - sign(Mjp))/2; % certo
                    else
                        Mjp_negative = -(Mjp - 1)^2/4 - beta * (Mjp^2 - 1)^2; % certo
                        pjp_negative = (Mjp - 1)^2 * (2 + Mjp)/4 - alpha * Mjp * (Mjp^2-1)^2; %certp
                    end

                    Mjhalf = Mj_plus + Mjp_negative; % certo

                    pjhalf = pj_plus * pL + pjp_negative * pR; %certo

                    Mjhalf_plus = (Mjhalf + abs(Mjhalf)) / 2; %certo
                    Mjhalf_negative = (Mjhalf - abs(Mjhalf)) / 2; % certo

                    phij = phi(qL, aL);  %certo
                    phijp = phi(qR, aR);  % certo

                    Fbarjhalf = aHalfjp * (Mjhalf_plus * phij + Mjhalf_negative * phijp) + [0; pjhalf; 0]; % certo


                    %% calculando Fbarjmhalf

                   Z_plus = [fluxLimiter((rho(j) - rho(j - 1)) / (rho(j - 1) - rho(j - 2))),...
                                     fluxLimiter((u(j) - u(j - 1)) / (u(j - 1) - u(j - 2))),...
                                     fluxLimiter((p(j) - p(j - 1)) / (p(j - 1) - p(j - 2)))];

                   Z_minus = [fluxLimiter((rho(j) - rho(j - 1)) / (rho(j + 1) - rho(j))),...
                                 fluxLimiter((u(j) - u(j - 1)) / (u(j + 1) - u(j))),...
                                 fluxLimiter((p(j) - p(j - 1)) / (p(j + 1) - p(j)))];

                    Z_plus(isnan(Z_plus)) = 1;
                    Z_minus(isnan(Z_minus)) = 1;

                    %      Z_plus(:) = 1;
                    % Z_minus(:) = 1;

                    VL = V(j - 1, :) +  Z_plus.*(V(j - 1, :)- V(j - 2, :))/2;

                    rhoL = VL(1);
                    uL = VL(2);
                    pL = VL(3);
                    EL = pL / ((gamma - 1) * rhoL) + uL^2/2;

                    VR = V(j, :) -  Z_minus.*(V(j + 1, :) - V(j, :))/2;

                    rhoR = VR(1);
                    uR = VR(2);
                    pR = VR(3);
                    ER = pR / ((gamma - 1) * rhoR) + uR^2/2;

                    qL = [rhoL; rhoL * uL; rhoL * EL];
                    qR = [rhoR; rhoR * uR; rhoR * ER];

                    aL = sqrt(gamma * pL / rhoL);
                    aR = sqrt(gamma * pR / rhoR);

                    aStarL = sqrt(2 * (gamma - 1) * H(qL, aL) / (gamma + 1));

                    aTildeL = aStarL * min(1, aStarL / abs(uL));

                    aStarR = sqrt(2 * (gamma - 1) * H(qR, aR) / (gamma + 1));

                    aTildeR = aStarR * min(1, aStarR / abs(uR));

                    aHalfjm = min(aTildeL, aTildeR);

                    Mjm = uL / aHalfjm;
                    Mj = uR/aHalfjm;

                    if abs(Mjm) >= 1
                        Mjm_positive = (Mjm + abs(Mjm)) / 2; % certo
                        pjm_positive = (1 + sign(Mjm)) / 2; % certo
                    else
                        Mjm_positive = (Mjm + 1)^2/4 + beta * (Mjm^2 - 1)^2; % certo
                        pjm_positive = (Mjm + 1)^2 * (2 - Mjm) / 4 + alpha * Mjm * (Mjm^2 - 1)^2; % certo
                    end

                    if abs(Mj) >= 1
                        Mj_negative = (Mj - abs(Mj))/2; % certo
                        pj_negative = (1 - sign(Mj)) / 2; %certo

                    else
                        Mj_negative = -(Mj - 1)^2 /4 - beta * (Mj^2 - 1)^2; % certo
                        pj_negative = (Mj - 1)^2 * (2 + Mj) / 4 - alpha * Mj*(Mj^2 - 1)^2;%certo
                    end

                    Mjmhalf = Mjm_positive + Mj_negative; %certo
                    pjmhalf = pjm_positive * pL + pj_negative * pR; % CERTO

                    Mjmhalf_positive = (Mjmhalf + abs(Mjmhalf)) / 2;
                    Mjmhalf_negative = (Mjmhalf - abs(Mjmhalf)) / 2;

                    phijm = phi(qL, aL); 
                    phij = phi(qR, aR);

                    Fbarjmhalf = aHalfjm * (Mjmhalf_positive* phijm + Mjmhalf_negative * phij) + [0;pjmhalf;0];

                    next_q(j, :) = q(j, :) - lamb * (Fbarjhalf - Fbarjmhalf)';
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

            CFL = 0.01;

            dt = CFL * dx /max(abs(u0) + a0);

            alpha = 3 / 16;
            beta = 1/8;

                 % coeficiente dos esquemas
            lamb = dt / dx;

            while 1

                % efetua um incremento temporal
                t = t + dt;

                if t * (L_ref / u_ref) > 1 
                    break
                end


                subplot(2,2, 1);
                plot(x, q(:, 1), 'k');

                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
                ax = gca;
                ax.FontSize = 16; 
                subplot(2,2, 2);
                plot(x, q(:, 2) ./ q(:, 1), 'r');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$u(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 3);
                plot(x, (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1))), 'b');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$p(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 4);
                plot(x, q(:, 3)./q(:, 1), 'g');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$E(x, t)$', 'Interpreter', 'latex');
                pause(0.0001);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                rho = q(:, 1);
                u = q(:, 2) ./ rho;
               % a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 
                p = (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1)));
                a = sqrt(gamma * p ./ rho);
                for j = 3:M
                    
                    aStarj = sqrt(2 * (gamma - 1) * H(q(j, :), a(j)) / (gamma + 1)); % certo
                    
                    aTildej = aStarj * min(1, aStarj / abs(u(j))); % certo

                    aStarjp = sqrt(2 * (gamma - 1)  * H(q(j + 1, :), a(j+1)) / (gamma + 1)); % certo

                    aTildejp = aStarjp * min(1, aStarjp / abs(u(j + 1))); % certo

                    aHalfjp = min(aTildej, aTildejp); %certo

                    Mj = u(j) / aHalfjp; % certo
                    Mjp = u(j+1) / aHalfjp; % certo

                    if abs(Mj) >= 1
                        Mj_plus = (Mj + abs(Mj))/2; % certo
                        pj_plus = (1 + sign(Mj))/2; % certo
                    else
                        Mj_plus = (Mj + 1)^2/4 + beta * (Mj^2 - 1)^2; % certo
                        pj_plus = (Mj + 1)^2 * (2 - Mj)/4 + alpha * Mj * (Mj^2-1)^2; % certo
                    end

                    if abs(Mjp) >= 1
                        Mjp_negative = (Mjp - abs(Mjp)) / 2; % certo
                        pjp_negative = (1 - sign(Mjp))/2; % certo
                    else
                        Mjp_negative = -(Mjp - 1)^2/4 - beta * (Mjp^2 - 1)^2; % certo
                        pjp_negative = (Mjp - 1)^2 * (2 + Mjp)/4 - alpha * Mjp * (Mjp^2-1)^2; %certp
                    end

                    Mjhalf = Mj_plus + Mjp_negative; % certo

                    pjhalf = pj_plus * p(j) + pjp_negative * p(j + 1); %certo

                    Mjhalf_plus = (Mjhalf + abs(Mjhalf)) / 2; %certo
                    Mjhalf_negative = (Mjhalf - abs(Mjhalf)) / 2; % certo

                    phij = phi(q(j, :), a(j));  %certo
                    phijp = phi(q(j+1, :), a(j+1));  % certo

                    Fbarjhalf = aHalfjp * (Mjhalf_plus * phij + Mjhalf_negative * phijp) + [0; pjhalf; 0]; % certo


                    %% calculando Fbarjmhalf

                    aStarjm = sqrt(2 * (gamma - 1) * H(q(j - 1, :), a(j - 1)) / (gamma + 1));

                    aTildejm = aStarjm * min(1, aStarjm / abs(u(j - 1)));

                    aHalfjm = min(aTildejm, aTildej);

                    Mjm = u(j - 1) / aHalfjm;
                    Mj = u(j)/aHalfjm;

                    if abs(Mjm) >= 1
                        Mjm_positive = (Mjm + abs(Mjm)) / 2; % certo
                        pjm_positive = (1 + sign(Mjm)) / 2; % certo
                    else
                        Mjm_positive = (Mjm + 1)^2/4 + beta * (Mjm^2 - 1)^2; % certo
                        pjm_positive = (Mjm + 1)^2 * (2 - Mjm) / 4 + alpha * Mjm * (Mjm^2 - 1)^2; % certo
                    end

                    if abs(Mj) >= 1
                        Mj_negative = (Mj - abs(Mj))/2; % certo
                        pj_negative = (1 - sign(Mj)) / 2; %certo

                    else
                        Mj_negative = -(Mj - 1)^2 /4 - beta * (Mj^2 - 1)^2; % certo
                        pj_negative = (Mj - 1)^2 * (2 + Mj) / 4 - alpha * Mj*(Mj^2 - 1)^2;%certo
                    end

                    Mjmhalf = Mjm_positive + Mj_negative; %certo
                    pjmhalf = pjm_positive * p(j - 1) + pj_negative * p(j); % CERTO

                    Mjmhalf_positive = (Mjmhalf + abs(Mjmhalf)) / 2;
                    Mjmhalf_negative = (Mjmhalf - abs(Mjmhalf)) / 2;

                    phijm = phi(q(j - 1, :), a(j-1)); 

                    Fbarjmhalf = aHalfjm * (Mjmhalf_positive* phijm + Mjmhalf_negative * phij) + [0;pjmhalf;0];

                    next_q(j, :) = q(j, :) - lamb * (Fbarjhalf - Fbarjmhalf)';
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


            CFL = 0.1;
            eigenvalues = [u0, u0 + a0, u0 - a0];
            dt =  CFL * dx / max(eigenvalues(:));




                 % coeficiente dos esquemas
            lamb = dt / dx;

            theta = 1/2;
            xi = 0;

            aux = q(:);
            dq = zeros(3 * M - 6, 1);%aux(7:3*M);%zeros(3 * M - 6, 1);%
            while 1

                % efetua um incremento temporal
                t = t + dt;

               
                if t * (L_ref / u_ref) > 1 
                    break
                end

                  subplot(2,2, 1);
                plot(x, q(:, 1), 'k');

                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
                ax = gca;
                ax.FontSize = 16; 
                subplot(2,2, 2);
                plot(x, q(:, 2) ./ q(:, 1), 'r');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$u(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 3);
                plot(x, (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1))), 'b');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$p(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 4);
                plot(x, q(:, 3)./q(:, 1), 'g');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$E(x, t)$', 'Interpreter', 'latex');
                pause(0.0001);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                rho = q(:, 1);
                u = q(:, 2) ./ rho;
                a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 

                B_1 = zeros(3 * M - 6, 3 * M - 6);
                B_2 = zeros(3 * M - 6, 3 * M - 6);
                f = zeros(3 * M - 6, 1);
                i = 1;
                k = 1;
                for j = 3:M

                     
                    % if u(j) >= 0 && u(j) <= a(j)
                    %     F_plus_jmm = Fpp(q(j - 2, :), a(j - 2));
                    %     F_plus_jm =  Fpp(q(j - 1, :), a(j - 1));
                    %     F_plus_j = Fpp(q(j, :), a(j));
                    % 
                    %     F_minus_jmm = Fnp(q(j - 2, :), a(j-2));
                    %     F_minus_jm = Fnp(q(j - 1, :), a(j-1));
                    %     F_minus_j = Fnp(q(j, :), a(j));
                    % 
                    % elseif u(j) >= -a(j) && u(j) < 0
                    % 
                    %     F_plus_jmm = Fpn(q(j - 2, :), a(j - 2));
                    %     F_plus_jm =  Fpn(q(j - 1, :), a(j - 1));
                    %     F_plus_j = Fpn(q(j, :), a(j));
                    % 
                    %     F_minus_jmm = Fnn(q(j-2, :), a(j-2));
                    %        F_minus_jm =Fnn(q(j-1, :), a(j-1));
                    %     F_minus_j = Fnn(q(j, :), a(j));
                    % 
                    % 
                    % elseif u(j) > a(j)
                    %     F_plus_jmm = F(q(j - 2, :));
                    %     F_plus_jm = F(q(j - 1, :));
                    %     F_plus_j = F(q(j, :));
                    % 
                    % 
                    %     F_minus_jmm = 0;
                    %     F_minus_jm = 0;
                    %     F_minus_j = 0;
                    % 
                    % 
                    % elseif u(j) < -a(j)
                    %         F_plus_jmm = 0;
                    %       F_plus_jm = 0;
                    %     F_plus_j = 0;
                    % 
                    % 
                    %     F_minus_jmm = F(q(j - 2, :));
                    %     F_minus_jm = F(q(j - 1, :));
                    %     F_minus_j = F(q(j, :));
                    % 
                    % end

                    % l1_plus_j = 0;
                    % l2_plus_j = 0;
                    % l3_plus_j = 0;
                    % 
                    % l1_minus_j = 0;
                    % l2_minus_j = 0;
                    % l3_minus_j = 0;
                    % 
                    % if u(j) > 0
                    %     l1_plus = u(j);
                    % else
                    %     l1_minus = u(j);
                    % end
                    % 
                    % if u(j) + a(j) > 0
                    %     l2_plus = u(j) + a(j);
                    % else
                    %     l2_minus = u(j) + a(j);
                    % end
                    % 
                    %  if u(j) - a(j) > 0
                    %     l3_plus = u(j) - a(j);
                    % else
                    %     l3_minus = u(j) - a(j);
                    % end
                    % 
                    % lp = [l1_plus, l2_plus, l3_plus];
                    % lm = [l1_minus, l2_minus, l3_minus];

                    lp = [u, u + a, u - a];
                    lm = lp;

                    lp = max(lp, 0);
                    lm = min(lm, 0);

                    if j == M
                        B_2(3*M - 8:3*M - 6, 3*M - 8:3*M - 6) = eye(3) - theta * lamb * A(u(j), a(j), lm(j,:)) / (1 + xi);
                    else
                        B_2(3 * (i - 1) + 1:3*i, (1:3) + 3 * (i - 1)) = eye(3) - theta * lamb * A(u(j), a(j), lm(j,:)) / (1 + xi);
                        B_2(3 * (i - 1) + 1:3*i, (1:3) + 3 * i) = theta * lamb * A(u(j + 1), a(j + 1),  lm(j + 1,:)) / (1 + xi);
                    end

                    if j == 3
                        B_1(1:3, 1:3) = eye(3) + theta * lamb * A(u(j), a(j), lp(j,:)) / (1 + xi);
                    else
                        B_1(3 * (i - 1) + 1:3*i, (1:3) + 3 * (i - 2)) = - theta * lamb * A(u(j - 1), a(j - 1), lp(j - 1,:))  / (1 + xi);
                        B_1(3 * (i - 1) + 1:3*i, (1:3) + 3 * (i - 1)) = eye(3) + theta * lamb * A(u(j), a(j), lp(j,:))  / (1 + xi);
                    end
                    % if j == 3
                    %     B_1(1:3, 1:3) = eye(3) + theta * lamb * A_p(rho(j),u(j), a(j)) / (1 + xi);
                    % elseif j == M
                    %     B_1(3*M - 8:3*M - 6, (3*M - 11):(3*M-9)) = - theta * lamb * A_p(rho(j-1),u(j - 1), a(j - 1)) / (1 + xi);
                    %     B_1(3*M - 8:3*M - 6, 3*M - 8:3*M - 6) = eye(3) + theta * lamb * A_p(rho(j),u(j), a(j)) / (1 + xi);
                    % else
                    %     B_1(3 * (i - 1) + 1:3*i, (1:3) + 3 * (k - 1)) = - theta * lamb * A_p(rho(j-1),u(j - 1), a(j - 1)) / (1 + xi);
                    %     B_1(3 * (i - 1) + 1:3*i, (4:6) + 3 * (k - 1)) = eye(3) + theta * lamb * A_p(rho(j),u(j), a(j)) / (1 + xi);
                    %     k = k + 1;
                    % end
    
                    %f(3 * (i - 1) + 1:3*i) =  3 * (F_plus_j - F_minus_j) - 4 * (F_plus_jm - F_minus_jm) + (F_plus_jmm - F_minus_jmm);
                  %  f(3 * (i - 1) + 1:3*i) =  3 * (Fplus(q(j, :), a(j)) - Fminus(q(j, :), a(j))) - 4 * (Fplus(q(j - 1, :), a(j - 1)) - Fminus(q(j - 1, :), a(j - 1))) + (Fplus(q(j - 2, :), a(j - 2)) - Fminus(q(j - 2, :), a(j - 2)));
                   %F = @(u, a, lambda, q) A(u,a,lambda) * q;    


                   f(3 * (i - 1) + 1:3*i) =  3 * (F(q(j, :), lp(j, :)) - F(q(j, :), lm(j,:))) ...
                                           - 4 * (F(q(j - 1, :), lp(j - 1, :)) - F(q(j + 1, :), lm(j + 1, :))) ...
                                           + 1 * (F(q(j - 2, :), lp(j - 2, :)) - F(q(j + 2, :), lm(j + 2, :)));

                    i = i + 1;
                end
                f = - lamb * f / (2 * (xi + 1));

                dq_star = linsolve(B_1, f);
                dq = linsolve(B_2, dq_star);

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

               %  u = q(:, 2) ./ rho;
               %  a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 
               % 
               % eigenvalues = [u, u + a, u - a];
               %  dt =  CFL * dx / max(eigenvalues(:));
               % 
               % 
               % 
               % 
               %   coeficiente dos esquemas
               %  lamb = dt / dx;



            end



        case EXP_STEGER_WARMING

              % determina o passo no tempo
            CFL = 0.1;

            eigenvalues = [u0, u0 + a0, u0 - a0];
            % lm = lp;
            % 
            % lp = max(lp, 0);
            % lp = lp(:);
            % lm = min(lm, 0);
            % 
            % lm = lm(:);
            % eigenvalues = abs([lp; lm]);

            % arr = cat(3, u0 + a0, u0 - a0, u0);
             dt = CFL * dx / max(eigenvalues(:));

            % coeficiente dos esquemas
            lamb = dt / dx;

                  k4 = 1/200;
            k2 = 1/4;

             while 1

                % efetua um incremento temporal
                t = t + dt;


                if t * (L_ref / u_ref) > 1 
                    break
                end


                subplot(2,2, 1);
                plot(x, q(:, 1), 'k');

                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
                ax = gca;
                ax.FontSize = 16; 
                subplot(2,2, 2);
                plot(x, q(:, 2) ./ q(:, 1), 'r');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$u(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 3);
                plot(x, (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1))), 'b');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$p(x, t)$', 'Interpreter', 'latex');
                subplot(2,2, 4);
                plot(x, q(:, 3)./q(:, 1), 'g');
                ax = gca;
                ax.FontSize = 16; 
                grid('on');
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$E(x, t)$', 'Interpreter', 'latex');
                pause(0.0001);

                % matriz contendo a solução em cada instante de tempo
                next_q = zeros(M + 2, 3);

                rho = q(:, 1);
                u = q(:, 2) ./ rho;
                a = sqrt(gamma * (gamma - 1) * (q(:, 3) ./ q(:, 1) - q(:, 2).^2 ./ (2 * q(:, 1).^2))); 
                p = (gamma - 1) * (q(:, 3) - q(:, 2).^2./(2*q(:, 1)));
               % a = sqrt(gamma*p./rho);  
                
                for j = 3:M

                     
                   %  if u(j) >= 0 && u(j) <= a(j)
                   %      F_plus_jmm = Fpp(q(j - 2, :), a(j - 2));
                   %      F_plus_jm =  Fpp(q(j - 1, :), a(j - 1));
                   %      F_plus_j = Fpp(q(j, :), a(j));
                   %      F_plus_jp = Fpp(q(j+1, :), a(j+1));
                   %      F_plus_jpp =Fpp(q(j+2, :), a(j+2));
                   % 
                   %      F_minus_jm = Fnp(q(j - 1, :), a(j-1));
                   %      F_minus_j = Fnp(q(j, :), a(j));
                   %      F_minus_jp = Fnp(q(j + 1, :), a(j + 1));
                   %      F_minus_jpp = Fnp(q(j + 2, :), a(j + 2));
                   % 
                   %  elseif u(j) >= -a(j) && u(j) < 0
                   %      F_plus_jmm = Fpn(q(j - 2, :), a(j - 2));
                   %      F_plus_jm =  Fpn(q(j - 1, :), a(j - 1));
                   %      F_plus_j = Fpn(q(j, :), a(j));
                   %      F_plus_jp=Fpn(q(j+1, :), a(j+1));
                   %      F_plus_jpp= Fpn(q(j+2, :), a(j+2));
                   % 
                   %         F_minus_jm =Fnn(q(j-1, :), a(j-1));
                   %      F_minus_j = Fnn(q(j, :), a(j));
                   %      F_minus_jp = Fnn(q(j + 1, :), a(j + 1));
                   %      F_minus_jpp = Fnn(q(j + 2, :), a(j + 2));
                   % 
                   %  elseif u(j) > a(j)
                   %      F_plus_jmm = F(q(j - 2, :));
                   %      F_plus_jm = F(q(j - 1, :));
                   %      F_plus_j = F(q(j, :));
                   %      F_plus_jp = F(q(j + 1, :));
                   %      F_plus_jpp = F(q(j + 2, :));
                   % 
                   %      F_minus_jm = 0;
                   %      F_minus_j = 0;
                   %      F_minus_jp = 0;
                   %      F_minus_jpp = 0;
                   % 
                   %  elseif u(j) < -a(j)
                   %          F_plus_jmm = 0;
                   %        F_plus_jm = 0;
                   %      F_plus_j = 0;
                   %      F_plus_jp  =0;
                   %      F_plus_jpp = 0;
                   % 
                   %      F_minus_jm = F(q(j - 1, :));
                   %      F_minus_j = F(q(j, :));
                   %      F_minus_jp = F(q(j+1, :));
                   %      F_minus_jpp = F(q(j+2, :));
                   %  end
                   % 
                   %  %Fjm = absA(u(j - 1), a(j - 1)) * (q(j - 1, :))';
                   %  %Fj  = absA(u(j), a(j)) * (q(j, :))';
                   %  %Fjp = absA(u(j + 1), a(j + 1)) * (q(j + 1, :))';
                   % 
                   % 
                   % 
                   %  next_q_bar_jp = q(j + 1, :) - lamb * (F_plus_jpp - F_plus_jp)' - lamb * (F_minus_jp - F_minus_j)';
                   % next_q_bar_j = q(j, :) - lamb * (F_plus_j - F_plus_jm)' - lamb * (F_minus_jp - F_minus_j)';
                   % next_q_bar_jm = q(j - 1, :) - lamb * (F_plus_jm - F_plus_jmm)' - lamb * (F_minus_j - F_minus_jm)';


              


                   % 
                   % if u(j) >= 0 && u(j) <= a(j)
                   %      % F_plus_bar_jmm = Fpp(q(j - 2, :), a(j - 2));
                   %      F_plus_bar_jm =  Fpp(next_q_bar_jm, a(j - 1));
                   %      F_plus_bar_j = Fpp(next_q_bar_j, a(j));
                   % 
                   %      % F_minus_bar_jm = Fnp(next_q_bar_jm, a(j-1));
                   %      F_minus_bar_j = Fnp(next_q_bar_j, a(j));
                   %      F_minus_bar_jp = Fnp(next_q_bar_jp, a(j + 1));
                   % 
                   %  elseif u(j) >= -a(j) && u(j) < 0
                   %      % F_plus_bar_jmm = Fpn(q(j - 2, :), a(j - 2));
                   %      F_plus_bar_jm =  Fpn(next_q_bar_jm, a(j - 1));
                   %      F_plus_bar_j = Fpn(next_q_bar_j, a(j));
                   % 
                   %         % F_minus_bar_jm =Fnn(next_q_bar_jm, a(j-1));
                   %      F_minus_bar_j = Fnn(next_q_bar_j, a(j));
                   %      F_minus_bar_jp = Fnn(next_q_bar_jp, a(j + 1));
                   % 
                   %  elseif u(j) > a(j)
                   %    %  F_plus_bar_jmm = F(q(j - 2, :));
                   %      F_plus_bar_jm = F(next_q_bar_jm);
                   %      F_plus_bar_j = F(next_q_bar_j);
                   % 
                   %     % F_minus_bar_jm = 0;
                   %      F_minus_bar_j = 0;
                   %      F_minus_bar_jp = 0;
                   % 
                   %  elseif u(j) < -a(j)
                   %      %    F_plus_bar_jmm = 0;
                   %        F_plus_bar_jm = 0;
                   %      F_plus_bar_j = 0;
                   % 
                   %      %F_minus_bar_jm = F(next_q_bar_jm);
                   %      F_minus_bar_j = F(next_q_bar_j);
                   %      F_minus_bar_jp = F(q(j+1, :));
                   % end
                   % 
                   %  Tjm = abs(p(j) - 2 * p(j-1) + p(j - 2))/abs(p(j)+ 2*p(j-1)+ p(j -2));
                   %  Tj = abs(p(j + 1) - 2 * p(j) + p(j - 1))/abs(p(j + 1)+ 2*p(j)+ p(j - 1));
                   %  Tjp = abs(p(j + 2) - 2 * p(j+1) + p(j))/abs(p(j + 2)+ 2*p(j+1)+ p(j));
                   % 
                   %  r2 = k2 * dt  * max([Tjp, Tj, Tjm]);
                   %  r4 = max(0, k4 * dt - r2);
                   % 
                   %  % longe da descontinuidade
                   %  if r4 ~= 0
                   %       % 2 ordem
                   %      next_q(j, :) = (q(j, :) + next_q_bar_j)/2 - lamb * (F_plus_j - 2 * F_plus_jm + F_plus_jmm + F_plus_bar_j - F_plus_bar_jm)'/2 + lamb *(F_minus_jpp - 2 * F_minus_jp + F_minus_j - (F_minus_bar_jp - F_minus_bar_j))'/2;
                   % 
                   %  else % perto da descontinuidade
                   % 
                   % 
                   %      %1 ordem
                   %  next_q(j, :) = q(j, :) - lamb  * (F_plus_j - F_plus_jm)' - lamb * (F_minus_jp - F_minus_j)';
                   % 
                   %  end


                   %next_q(j, :) = q(j, :) - lamb  * (Fplus(q(j, :), a(j)) - Fplus(q(j - 1, :), a(j - 1)))' - lamb * (Fminus(q(j + 1, :), a(j + 1)) - Fminus(q(j, :), a(j)))';

                         lp = [u, u + a, u - a];
                    lm = lp;

                    lp = max(lp, 0);
                    lm = min(lm, 0);

                   

                    next_q_bar_jm = q(j - 1, :) ...
                                   - lamb * (F(q(j - 1, :), lp(j - 1, :)) - F(q(j - 2, :), lp(j - 2, :)))' ...
                                   - lamb * (F(q(j, :), lm(j, :)) - F(q(j - 1, :), lm(j - 1, :)))';
                    
                    next_q_bar_j = q(j, :) ...
                                  - lamb * (F(q(j, :), lp(j, :)) - F(q(j - 1, :), lp(j - 1, :)))' ...
                                  - lamb * (F(q(j + 1, :), lm(j + 1, :)) - F(q(j, :), lm(j, :)))';
                    
                    next_q_bar_jp = q(j + 1, :) ...
                                  - lamb * (F(q(j + 1, :), lp(j + 1, :)) - F(q(j, :), lp(j, :)))' ...
                                  - lamb * (F(q(j + 2, :), lm(j + 2, :)) - F(q(j + 1, :), lm(j + 1, :)))';

    
                    %     Tjm = abs(p(j) - 2 * p(j-1) + p(j - 2))/abs(p(j)+ 2*p(j-1)+ p(j -2));
                    % Tj = abs(p(j + 1) - 2 * p(j) + p(j - 1))/abs(p(j + 1)+ 2*p(j)+ p(j - 1));
                    % Tjp = abs(p(j + 2) - 2 * p(j+1) + p(j))/abs(p(j + 2)+ 2*p(j+1)+ p(j));
                    % 
                    % r2 = k2 * dt  * max([Tjp, Tj, Tjm]);
                    % r4 = max(0, k4 * dt - r2);

 %                      % longe da descontinuidade
 %                    if r4 ~= 0
 % 
 %                    %2 ordem
 %                    next_q(j, :) = (q(j, :) + next_q_bar_j) / 2 ...
 %                                    - lamb * (F(q(j, :), lp(j, :)) - 2 * F(q(j - 1, :), lp(j - 1, :)) + F(q(j - 2, :), lp(j - 2, :)) + F(next_q_bar_j, lp(j, :)) - F(next_q_bar_jm, lp(j - 1, :)))'/2 ...
 %                                    + lamb *(F(q(j + 2, :), lm(j + 2, :)) - 2 * F(q(j + 1, :), lm(j + 1, :)) + F(q(j, :), lm(j, :)) - (F(next_q_bar_jp, lm(j + 1, :)) - F(next_q_bar_j, lm(j, :))))'/2;
 % 
 % 
 %                    else % perto da descontinuidade
 % 
 % 
 % 
 % 1 ordem
                  next_q(j, :) = q(j, :) - lamb  * (F(q(j, :), lp(j, :)) - F(q(j - 1, :), lp(j - 1, :)))' - lamb * (F(q(j + 1, :), lm(j + 1, :)) - F(q(j, :), lm(j, :)))';
 % 
 %                    end

                    % 1 ORDEM
               %     next_q(j, :) = q(j, :) - 0.5 * lamb * (3 * F(q(j, :), lp(j, :)) - 4 * F(q(j - 1, :), lp(j - 1, :)) + F(q(j - 2, :), lp(j - 2, :)))' ...
              %                             - 0.5 * lamb * (- 3 * F(q(j, :), lm(j, :)) + 4 * F(q(j + 1, :), lm(j + 1, :)) - F(q(j + 2, :), lm(j + 2, :)))';
                   
                   
                   
                   % 2 ordem
                    % next_q(j, :) = (q(j, :) + next_q_bar_j) / 2 ...
                    %                - lamb * (F(q(j, :), lp(j, :)) - 2 * F(q(j - 1, :), lp(j - 1, :)) + F(q(j - 2, :), lp(j - 2, :)) + F(next_q_bar_j, lp(j, :)) - F(next_q_bar_jm, lp(j - 1, :)))'/2 ...
                    %                + lamb *(F(q(j + 2, :), lm(j + 2, :)) - 2 * F(q(j + 1, :), lm(j + 1, :)) + F(q(j, :), lm(j, :)) - (F(next_q_bar_jp, lm(j + 1, :)) - F(next_q_bar_j, lm(j, :))))'/2;
                    % 
                   
                   %F(u(j), a(j), lp(j,:), q(j, :))



                   % next_q(j, :) = q(j, :) - lamb * (F(q(j + 1, :)) - F(q(j - 1, :)))' / 2 + ...
                      %  lamb * (Fjp - 2 * Fj + Fjm)' / 2;
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
end