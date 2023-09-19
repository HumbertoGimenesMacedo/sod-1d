clear, clc
close all

% variáveis que identificam o esquema numérico
LAX_WENDROFF = 1;
MAC_COMARCK = 2;
EXP_BEAM_WARMING = 3;
IMP_BEAM_WARMING = 4;
EXP_STEGER_WARMING = 5;
IMP_STEGER_WARMING = 8;
AUSM_PLUS = 6;
VAN_LEER = 7;
ROE = 9;
gamma = 1.4;

% estipula o método
METHOD = ROE;

% ordem do método
order = 2;

scheme_name = "";
switch (METHOD)
    case LAX_WENDROFF
        scheme_name = "Lax-Wendroff";
    case MAC_COMARCK
        scheme_name = "Mac Comarck";
    case EXP_BEAM_WARMING
        scheme_name = "Explicit Beam-Warming";
    case IMP_BEAM_WARMING
        scheme_name = "Implicit Beam-Warming";
    case EXP_STEGER_WARMING
         if order == 1
            scheme_name = "First Order Explicit Steger-Warming";
        elseif order == 2
            scheme_name = "Second Order Explicit Steger-Warming";
         end
    case IMP_STEGER_WARMING
        scheme_name = "Second Order Implicit Steger-Warming";
    case VAN_LEER
        if order == 1
            scheme_name = "First Order Van Leer";
        elseif order == 2
            scheme_name = "Second Order Van Leer";
        end
    case AUSM_PLUS
        if order == 1
            scheme_name = "First Order Liou Scheme";
        elseif order == 2
            scheme_name = "Second Order Liou Scheme";
        end   
    case ROE
        if order == 1
            scheme_name = "First Order Roe Scheme";
        elseif order == 2
            scheme_name = "Second Order Roe Scheme";
        end   
end

% razão entre as pressões em ambos os lados do diafragma
pressureRatio = 100;

% constante dos gases ideais
R = 287.052874;

% estado do gás à direita do diafragma (atmosfera padrão)
rho_r = 1.225;
u_r = 0;
p_r = 101325;
T_r = 288.15;

% estado do gás à esquerda do diafragma
p_l = pressureRatio * p_r;
T_l = T_r;
u_l = 0;
rho_l = p_l / (R * T_l);

% resolve de forma numérica
[num_x, num_rho, num_u, num_E, num_p, totalTime] = solveShockProblem(METHOD, order, rho_l, u_l, ...
                                    p_l, rho_r, u_r, p_r);

% Create a table with the data and variable names
T = table(num_x, num_rho, num_u, num_E, num_p, ones(length(num_x), 1) * totalTime, 'VariableNames', {'x', 'rho', 'u', 'E', 'p', 'time'});
% Write data to text file
writetable(T, strcat(scheme_name, '_', string(pressureRatio), '.txt'))

% obtém solução analítica (razão 5)
data = readmatrix(strcat('3000_d_pressure_ratio_', string(pressureRatio), '.txt'));
ana_x = data(:, 1);
ana_rho = data(:, 2);
ana_u = data(:, 3);
ana_E = data(:, 4) ./ (data(:, 2) * (gamma - 1)) + data(:, 3).^2 / 2;
ana_p = data(:, 4);

sgtitle(strcat(scheme_name, sprintf(' ($t = %0.2f$)', 1.0)), 'Interpreter', 'latex'); 
set(gcf, 'WindowState', 'maximized');

% densidade
subplot(2, 2, 1);
plot(num_x, num_rho, 'r', 'LineWidth', 1);
hold('on');
plot(ana_x, ana_rho, 'b', 'LineWidth', 1);
grid('on');
axis('tight')
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 16; 

% velocidade
subplot(2, 2, 2);
plot(num_x, num_u, 'r', 'LineWidth', 1);
hold('on');
plot(ana_x, ana_u, 'b', 'LineWidth', 1);
grid('on');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u(x, t)$', 'Interpreter', 'latex');
legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 16; 

% pressão
subplot(2,2, 3);
plot(num_x, num_p, 'r', 'LineWidth', 1);
hold('on');
plot(ana_x, ana_p, 'b', 'LineWidth', 1);
grid('on');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$p(x, t)$', 'Interpreter', 'latex');
legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 16; 

% energia
subplot(2,2, 4);
plot(num_x, num_E, 'r', 'LineWidth', 1);
hold('on');
plot(ana_x, ana_E, 'b', 'LineWidth', 1);
grid('on');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$E(x, t)$', 'Interpreter', 'latex');
legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 16; 