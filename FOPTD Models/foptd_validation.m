
clear; clc; close all;

% Loading identification and simulation data
id  = load('bio_foptd_identification_all.mat');
sim = load('bioreactor_data.mat');

t          = sim.t_vec(:);
x_nl       = sim.x_nl;
x_nl_noisy = sim.x_nl_noisy;
x_ss       = sim.x_ss(:);
D_0         = sim.D0;
delta_D     = sim.deltaD;

Kp    = id.K_p;
Tau_p  = id.Tau_p;
Theta_vec = id.Theta;

n_states    = 3;
state_names = {'X (biomass)','S (substrate)','P (product)'};

D_vec = D_0 * ones(size(t));
D_vec(t >= 100) = D_0 - 1*delta_D;
D_vec(t >= 200) = D_0 - 2*delta_D;
D_vec(t >= 300) = D_0 - 3*delta_D;
u = D_vec - D_0;

Ts = t(2) - t(1);

figure(1); clf;

for y_index = 1:n_states

    y_data = x_nl_noisy(:, y_index);
    y_ss   = x_ss(y_index);

    K_p    = Kp(y_index);
    Tau   = Tau_p(y_index);
    Theta = Theta_vec(y_index);

    % FOPTD simulation for this state
    Nd = max(0, round(Theta_vec / Ts));
    u_delayed = zeros(size(u));
    if Nd < length(u)
        u_delayed(Nd+1:end) = u(1:end-Nd);
    end

    y_dev = zeros(size(t));
    for k = 1:length(t)-1
        y_dev(k+1) = y_dev(k) + (Ts/Tau)*(-y_dev(k) + K_p*u_delayed(k));
    end

    y_model = y_ss + y_dev;

    % Metrics 
    err  = y_data - y_model;
    SSE  = sum(err.^2);
    ybar = mean(y_data);
    SST  = sum((y_data - ybar).^2);
    R2   = 1 - SSE/SST;

    fprintf('\n Validation: D to %s (state %d) \n', ...
            state_names{y_index}, y_index);
    fprintf('  Kp    = %.6g\n', K_p);
    fprintf('  Tau_p  = %.6g\n', Tau);
    fprintf('  Theta = %.6g\n', Theta_vec);
    fprintf('  SSE       = %.6g\n', SSE);
    fprintf('  R_square        = %.6g\n', R2);

    % Subplot for this state
    subplot(3,1,y_index);
   
    plot(t, y_data, 'r.', 'MarkerSize', 5);hold on;
    plot(t, y_model, ' b--', 'LineWidth', 1.5);
    yline(y_ss, ':k', 'LineWidth', 1);
    xlabel('Time (h)');
    ylabel('State value (g/L)');
    title(sprintf('D to %s  (R^2 = %.3f)', state_names{y_index}, R2));
    legend('Noisy','FOPTD','SS','Location','best');

    grid on;

end
