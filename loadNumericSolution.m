function [num_x, num_rho, num_u, num_E, num_p] = loadNumericSolution(METHOD, scheme_name, order, pressureRatio) 
    % 
    % scheme_name = "";
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
    %             scheme_name = "Second Order Explicit Harten Scheme";
    %         end
    % end

    if order == 1
        S = readmatrix(strcat('First Order (Outputs)/', scheme_name, '_', string(pressureRatio), '.txt'));
    elseif order == 2
        S = readmatrix(strcat('Second Order (Outputs)/', scheme_name, '_', string(pressureRatio), '.txt'));
    end

    num_x = S(:, 1);
    num_rho = S(:, 2);
    num_u = S(:, 3);
    num_E = S(:, 4);
    num_p = S(:, 5);
end