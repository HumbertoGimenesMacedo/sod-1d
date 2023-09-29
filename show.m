function show(x, q, gamma)
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
end