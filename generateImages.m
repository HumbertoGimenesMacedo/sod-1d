function generateImages(METHOD, scheme_name, order, pressureRatio, gamma)
    % 
    % % variáveis que identificam o esquema numérico
    % LAX_WENDROFF = 1;
    % MAC_COMARCK = 2;
    % EXP_BEAM_WARMING = 3;
    % IMP_BEAM_WARMING = 4;
    % EXP_STEGER_WARMING = 5;
    % IMP_STEGER_WARMING = 8;
    % AUSM_PLUS = 6;
    % VAN_LEER = 7;
    % ROE = 9;
    % HARTEN = 10;
    % 
    % switch (METHOD)
    %     case LAX_WENDROFF
    %         scheme_name = "Second Order Lax-Wendroff";
    %     case MAC_COMARCK
    %         scheme_name = "Second Order MacComarck";
    %     case EXP_BEAM_WARMING
    %         scheme_name = "Second Order Explicit Beam-Warming";
    %     case IMP_BEAM_WARMING
    %         scheme_name = "Second Order Implicit Beam-Warming";
    %     case EXP_STEGER_WARMING
    %          if order == 1
    %             scheme_name = "First Order Explicit Steger-Warming";
    %         elseif order == 2
    %             scheme_name = "Second Order Explicit Steger-Warming";
    %          end
    %     case IMP_STEGER_WARMING
    %         scheme_name = "Second Order Implicit Steger-Warming";
    %     case VAN_LEER
    %         if order == 1
    %             scheme_name = "First Order Van Leer";
    %         elseif order == 2
    %             scheme_name = "Second Order Van Leer";
    %         end
    %     case AUSM_PLUS
    %         if order == 1
    %             scheme_name = "First Order Liou Scheme";
    %         elseif order == 2
    %             scheme_name = "Second Order Liou Scheme";
    %         end   
    %     case ROE
    %         if order == 1
    %             scheme_name = "First Order Roe Scheme";
    %         elseif order == 2
    %             scheme_name = "Second Order Roe Scheme";
    %         end   
    %     case HARTEN
    %         if order == 1
    %             scheme_name = "First Order Explicit Harten Scheme";
    %         elseif order == 2
                 scheme_name = "Second Order Implicit Harten Scheme";
                 METHOD = 10;
                 order = 2;
                 pressureRatio=100;
    %         end
    % end
    
    [num_x, num_rho, num_u, num_E, num_p] = loadNumericSolution(METHOD, scheme_name, order, pressureRatio);
    
    
    % obtém solução analítica (razão 5)
    data = readmatrix(strcat('Analytic Solution/3000_d_pressure_ratio_', string(pressureRatio), '.txt'));
    ana_x = data(:, 1);
    ana_rho = data(:, 2);
    ana_u = data(:, 3);
    ana_E = data(:, 4) ./ (data(:, 2) * (gamma - 1)) + data(:, 3).^2 / 2;
    ana_p = data(:, 4);

    %figure();
    figure('units','normalized','outerposition',[0 0 1 1], 'visible','off')
    %figure;
    sgtitle(strcat(scheme_name, sprintf(' ($t = %0.2f$)', 1.0)), 'Interpreter', 'latex'); 
    set(gcf, 'WindowState', 'maximized');
    
    % densidade
    subplot(2, 2, 1);
    plot(num_x, num_rho, 'r', 'LineWidth', 1);
    hold('on');
    plot(ana_x, ana_rho, 'b--', 'LineWidth', 1);
    grid('on');
    axis('tight')
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
    legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 16; 
    
    
    % create a new pair of axes inside current figure
    axes('units', 'normalized', 'position',[.15 .62 0.08 0.08])
    box on % put box around new pair of axes
    indexOfInterest = (num_x < 2000) & (num_x > 500); % range of t near perturbation
    plot(num_x(indexOfInterest),num_rho(indexOfInterest), 'r', 'LineWidth',1) % plot on new axes
    % xlabel('$x$', 'Interpreter', 'latex');
    %ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
    set(gca,'ytick',[])
        ax = gca;
    ax.FontSize = 13; 
    axis tight
    
    % velocidade
    subplot(2, 2, 2);
    plot(num_x, num_u, 'r', 'LineWidth', 1);
    hold('on');
    plot(ana_x, ana_u, 'b--', 'LineWidth', 1);
    grid('on');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$u(x, t)$', 'Interpreter', 'latex');
    legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex', 'Location', 'northwest');
    ax = gca;
    ax.FontSize = 16; 


    % create a new pair of axes inside current figure
    axes('units', 'normalized', 'position',[.73 .62 0.08 0.08])
    box on % put box around new pair of axes
    indexOfInterest = (num_x < 2000) & (num_x > 500); % range of t near perturbation
    plot(num_x(indexOfInterest),num_u(indexOfInterest), 'r', 'LineWidth',1) % plot on new axes
    % xlabel('$x$', 'Interpreter', 'latex');
    %ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 13; 
    set(gca,'ytick',[])
    axis tight
    
    % pressão
    subplot(2,2, 3);
    plot(num_x, num_p, 'r', 'LineWidth', 1);
    hold('on');
    plot(ana_x, ana_p, 'b--', 'LineWidth', 1);
    grid('on');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$p(x, t)$', 'Interpreter', 'latex');
    legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 16; 

    % create a new pair of axes inside current figure
    axes('units', 'normalized', 'position',[.15 .17 0.08 0.08])
    box on % put box around new pair of axes
    indexOfInterest = (num_x < 2000) & (num_x > 400); % range of t near perturbation
    plot(num_x(indexOfInterest),num_p(indexOfInterest), 'r', 'LineWidth',1) % plot on new axes
    % xlabel('$x$', 'Interpreter', 'latex');
    %ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 13; 
    set(gca,'ytick',[])
    axis tight



    % energia
    subplot(2,2, 4);
    plot(num_x, num_E, 'r', 'LineWidth', 1);
    hold('on');
    plot(ana_x, ana_E, 'b--', 'LineWidth', 1);
    grid('on');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$E(x, t)$', 'Interpreter', 'latex');
    legend({"Num\'erica", "Exata"}, 'Interpreter', 'latex', 'Location', 'northwest');
    ax = gca;
    ax.FontSize = 16; 
    %set(gcf,'PaperSize',[10 10]);
    set(gcf, 'Color', 'w')


    % create a new pair of axes inside current figure
    axes('units', 'normalized', 'position',[.65 .25 0.08 0.08])
    box on % put box around new pair of axes
    indexOfInterest = (num_x < 2000) & (num_x > 400); % range of t near perturbation
    plot(num_x(indexOfInterest),num_E(indexOfInterest), 'r', 'LineWidth',1) % plot on new axes
    % xlabel('$x$', 'Interpreter', 'latex');
    %ylabel('$\rho (x,t)$', 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 13; 
    set(gca,'ytick',[])
    axis tight
   

    if order == 1
         export_fig(gcf, '-m3', strcat('First Order (Outputs)/', scheme_name, '_', string(pressureRatio),'_z', '.pdf'))
    elseif order == 2
         export_fig(gcf, '-m3', strcat('Second Order (Outputs)/', scheme_name, '_', string(pressureRatio), '_z', '.pdf'))
    end
end