function closed_loop()
    clear; 
    clc; close all;
    SetGraphics(); 

    % Parameters
    Y_xs = 0.4;   % g/g
    b = 0.2;  % h^(-1)
    P_m  = 50;    % g/L
    K_i  = 22;    % g/L
    a = 2.2; % g/
    mu_m = 0.48;  % h^(-1)
    K_m  = 1.2;   % g/L
    Sf_0 = 20;    % g/L
    D_0  = 0.202; % h^(-1) nominal point

    % plot parameters
    sample_time       = 0.04;   
    n        = 4000;  
    plt_every     = 15;     
    pause_t     = 0.05;
    settle_time = 0;      

    % Measurement noise
    rng(0);
    % Non-Linear Model
    mu = @(S,P) (mu_m * (1 - P / P_m) .* S) ./ (K_m + S + (S.^2) / K_i);

    f  = @(x,u,Sf_par) [ ...
        -u * x(1) + mu(x(2), x(3)) * x(1);             % dX/dt
         u * (Sf_par - x(2)) - (1 / Y_xs) * mu(x(2), x(3)) * x(1); % dS/dt
        -u * x(3) + (a * mu(x(2), x(3)) + b) * x(1)  % dP/dt
    ];
    X_steady = 5.9956 ; S_steady = 5.0109; P_steady = 19.1267;
    x_steady = [X_steady; S_steady; P_steady];
    fprintf('SS at D0=%.6f, Sf=%.3f -> X=%.6f, S=%.6f, P=%.6f\n', ...
        D_0, Sf_0, X_steady, S_steady, P_steady);

    % Choose Control Variable
    fprintf('\nSelect controlled variable (CV):\n');
    fprintf('  1 = Biomass X\n');
    fprintf('  2 = Substrate S\n');
    fprintf('  3 = Product P\n');
    User_choice = input('Choice (1/2/3): ');
    switch User_choice
        case 1
            y_index  = 1;
            cv_name  = 'X';
            cv_title = 'Biomass X';
            cv_ss    = X_steady;

            % FOPTD parameters for D -> X
            Kp    = -36.92;   
            Tau_P =  5.598;  
            Tau_D =  0.0799;    

        case 2
            y_index  = 2;
            cv_name  = 'S';
            cv_title = 'Substrate Conc. S';
            cv_ss    = S_steady;

            % FOPTD parameters for D -> S 
            Kp    =  93.50;   
            Tau_P =  5.598;    
            Tau_D =  0.0799;  

        case 3
            y_index  = 3;
            cv_name  = 'P';
            cv_title = 'Product Conc. P';
            cv_ss    = P_steady;

            % FOPTD parameters for D -> P 
            Kp    = -147.9;  
            Tau_P =  6.775;    
            Tau_D =  0.0799;  

        otherwise
            error('Invalid user choice. Must be 1, 2, or 3.');
    end

    fprintf('\nUsing FOPTD (D -> %s): Kp=%.5g, Tau_P=%.5g h, Tau_D=%.5g h\n', ...
        cv_name, Kp, Tau_P, Tau_D);

    % Direct-Syntheis Tuning
    Tau = Tau_P / 2;   
    K_c    = (1 / Kp) * (Tau_P / (Tau + Tau_D)); 
    Tau_I  = Tau_P;      
    T_deri = 0.9;            

    fprintf('\nSelect mode:\n  1 = Regulatory \n  2 = Servo \n', cv_name);
    Choice1 = input('Choice (1/2): ');
    fprintf('\nSelect controller:\n  1 = P\n  2 = PI\n  3 = PID\n');
    Choice2 = input('Choice (1/2/3): ');

    % Time vector
    vec_time = (1:n)' * sample_time;

    % Servo setpoint 
    setpts = cv_ss * ones(n,1);
    time_sample  = [300 700 1100 1500 1900 2300 2700];
    sp_offsets  = [+0.5  -0.3  +0.8  -0.6  +0.4  -0.8   0.0];   % offsets (g/L)
    for i = 1:numel(time_sample)
        setpts(time_sample(i):end) = cv_ss + sp_offsets(i);
    end

    % Regulatory disturbance 
    Sf_sample = Sf_0 * ones(n,1);
    reg_step   = [500  1500  2500];
    Sf_offsets = [+2   -2    +1];   

    
    for ii = 1:numel(reg_step)
        Sf_sample(reg_step(ii):end) = max(0.1, Sf_0 + Sf_offsets(ii));
    end

    % Initializing the values
    y_tru  = zeros(n,1);   % true CV
    y_setp     = zeros(n,1);   % setpoint
    D_hist   = zeros(n,1);   % MV trajectory
    Sf_hist  = zeros(n,1);   % disturbance trajectory
    err_hist = zeros(n,1);   % control error

    % Live plotting setup 
    fig = figure('Name','Closed-loop');
    ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
    h_y  = plot(ax1, vec_time(1), cv_ss, 'k-');      % plant output
    h_sp = plot(ax1, vec_time(1), cv_ss, 'r-');      % setpoint
    xlabel(ax1,'Time (h)'); ylabel(ax1, sprintf('%s (g/L)', cv_name));
    title(ax1, sprintf('CV (%s) vs Setpoint', cv_title));
    legend(ax1, {'output','desired sp'}, 'Location','best');

    ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
    h_D  = plot(ax2, vec_time(1), D_0, 'g-');
    xlabel(ax2,'Time (h)'); ylabel(ax2,'D (h^{-1})');
    title(ax2,'MV');

    ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
    h_Sf = plot(ax3, vec_time(1), Sf_0, 'k-');
    xlabel(ax3,'Time (h)'); ylabel(ax3,'S_f (g/L)');
    title(ax3,'DV');

    xlim(ax1,[0 vec_time(end)]);
    xlim(ax2,[0 vec_time(end)]);
    xlim(ax3,[0 vec_time(end)]);
    drawnow;

    % CLOSED LOOP SIMULATION
    xcl   = x_steady(:);       
    Dcl  = D_0;            
    D_ss = D_0;            

    
    D_min = -1;  
    D_max = 1;

    
    Iaccum   = 0;        % integral accumulator
    prev_y = cv_ss;    % previous measurement for derivative on PV

    for k = 1:n

        if k <= settle_time
            % Initial steady-state hold
            y_setp(k) = cv_ss;
            Sf_k    = Sf_0;

            [~, xtra] = ode45(@(t,xx) f(xx, D_ss, Sf_k), [0, sample_time], xcl);
            xcl = xtra(end,:).';

            y_tru(k)  = xcl(y_index);
            D_hist(k)  = D_ss;
            Sf_hist(k) = Sf_k;
            err_hist(k)= 0;

            if (mod(k,plt_every) == 0) || (k <= settle_time) || (k == n)
                set(h_y,  'XData', vec_time(1:k), 'YData', y_tru(1:k));
                set(h_sp, 'XData', vec_time(1:k), 'YData', y_setp(1:k));
                set(h_D,  'XData', vec_time(1:k), 'YData', D_hist(1:k));
                set(h_Sf, 'XData', vec_time(1:k), 'YData', Sf_hist(1:k));
                drawnow limitrate;
                if pause_t > 0
                    pause(pause_t);
                end
            end
            continue;
        end

        % Mode-dependent setpoint and disturbance
        if Choice1 == 1   % Regulatory
            y_setp(k) = cv_ss;
            Sf_k    = Sf_sample(k);
        else              % Servo
            y_setp(k) = setpts(k);
            Sf_k    = Sf_0;
        end

        % Plant propagation with current MV
        [~, xtra] = ode45(@(t,xx) f(xx, Dcl, Sf_k), [0, sample_time], xcl);
        xcl = xtra(end,:).';

        yk_true    = xcl(y_index);           % true CV
        yk         = yk_true;       % measured CV
        y_tru(k)  = yk_true;
        D_hist(k)  = Dcl;
        Sf_hist(k) = Sf_k;

        % Controllers
        err = y_setp(k) - yk;
        err_hist(k) = err;

        dy = 0;   % derivative term (on PV), default 0

        switch Choice2
            case 1  % P controller
                Dk_cmd = D_ss + K_c * err;

            case 2  % PI controller
                Iaccum   = Iaccum + (sample_time / Tau_I) * err;
                Dk_cmd = D_ss + K_c * (err + Iaccum);

            case 3  % PID controller (D on measurement/PV) 
                dy     = (prev_y - yk) / sample_time;   % d(-y)/dt
                Iaccum   = Iaccum + (sample_time / Tau_I) * err;
                Dk_cmd = D_ss + K_c * (err + Iaccum + T_deri * dy);
        end

        % update previous measurement for next iteration 
        prev_y = yk;

        % Saturation + anti-windup 
        Dk_unsaturate = Dk_cmd;
        Dcl = min(max(Dk_cmd, D_min), D_max);

        if Dcl ~= Dk_unsaturate && Choice2 >= 2 && abs(K_c) > 1e-12
            % Back-calculate Iacc to be consistent with saturated Dk
            if Choice2 == 2         % PI
                Iaccum = (Dcl - D_ss) / K_c - err;
            elseif Choice2 == 3     % PID with D on PV
                Iaccum = (Dcl - D_ss) / K_c - err - T_deri * dy;
            end
        end
        % Live plot update
        if (mod(k,plt_every) == 0) || (k <= settle_time) || (k == n)
            set(h_y,  'XData', vec_time(1:k), 'YData', y_tru(1:k));
            set(h_sp, 'XData', vec_time(1:k), 'YData', y_setp(1:k));
            set(h_D,  'XData', vec_time(1:k), 'YData', D_hist(1:k));
            set(h_Sf, 'XData', vec_time(1:k), 'YData', Sf_hist(1:k));
            drawnow limitrate;
            if pause_t > 0
                pause(pause_t);
            end
        end
    end
end
function SetGraphics()
    set(0,'DefaultLineLineWidth', 2)
    set(0,'DefaultaxesLineWidth', 2)
    set(0,'DefaultLineMarkerSize',10)
    set(0,'DefaultaxesFontSize', 18)
    set(0,'DefaultTextFontSize', 18)
    set(0,'DefaultaxesFontName', 'arial')
end
