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
AUSM_PLUS_SEC = 9;
gamma = 1.4;

% estipula o método
METHOD = EXP_STEGER_WARMING;

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
        scheme_name = "Explicit Steger-Warming";
    case IMP_STEGER_WARMING
        scheme_name = "Implicit Steger-Warming";
end

% razão entre as pressões em ambos os lados do diafragma
pressureRatio = 50;

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
[num_x, num_rho, num_u, num_E, num_p] = solveShockProblem(METHOD, rho_l, u_l, ...
                                    p_l, rho_r, u_r, p_r);

% Create a table with the data and variable names
T = table(num_x, num_rho, num_u, num_E, num_p, 'VariableNames', {'x', 'rho', 'u', 'E', 'p'});
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