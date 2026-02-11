clear; 
clc; 
close all;

% Given Parameters
Y_xs = 0.4;    % g/g
beta = 0.2;    % h^-1
P_m = 50;      % g/L
K_i = 22;      % g/L
alpha = 2.2;   % g/g
mu_m = 0.48;   % h^-1
K_m = 1.2;     % g/L
S_f = 20;      % g/L
d0 = 0.202;    % dilution rate at nominal operating point

%% MODEL EQUATIONS 
mu = @(S,P) (mu_m*(1 - P/P_m) * S) ./ (K_m + S + (S.^2)/K_i);
func = @(t,f,d,Sf)[ ...
    -d*f(1) + mu(f(2),f(3))*f(1); ...
     d*(Sf-f(2)) - (1/Y_xs)*mu(f(2),f(3))*f(1); ...
    -d*f(3) + (alpha*mu(f(2),f(3))+beta)*f(1)];

% Finding steady-state by simulating for long time
f_i = [0.1, 2, 0];  % initial guess
tspan = [0 500];    
[tsim, xsim] = ode45(@(t,x) func(t,x,d0,S_f), tspan, f_i);

% Plot transient response
figure('Name','Bioreactor Transient Simulation','NumberTitle','off');
plot(tsim, xsim(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(tsim, xsim(:,2), 'r-', 'LineWidth', 1.5);
plot(tsim, xsim(:,3), 'g-', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Concentration (g/L)');
title('Bioreactor Transient Response (Approach to Steady-State)');
legend('X (biomass)','S (substrate)','P (product)','Location','best');
grid on; box on; hold off;

% Steady state
x_ss = xsim(end,:)';
X_ss = x_ss(1); S_ss = x_ss(2); P_ss = x_ss(3);
fprintf('Steady-state concentrations:\nX_ss = %.5f, S_ss = %.5f, P_ss = %.5f\n', X_ss,S_ss,P_ss);

%% LINEARIZATION 
syms X S P d Sf

mu_sym = mu_m * (1 - P/P_m) * (S / (K_m + S + (S^2)/K_i));
f1 = -d*X + mu_sym*X;
f2 = d*(Sf - S) - (1/Y_xs)*mu_sym*X;
f3 = -d*P + (alpha*mu_sym + beta)*X;
f_vec = [f1; f2; f3];

% Jacobian (A matrix)
J_sym = jacobian(f_vec, [X S P]);
J = double(subs(J_sym, {X, S, P, d, Sf}, {X_ss, S_ss, P_ss, d0, S_f}));

% Input matrices
I_d_sym = diff(f_vec, d);   % effect of dilution rate
I_Sf_sym = diff(f_vec, Sf); % effect of feed substrate

I_d = double(subs(I_d_sym, {X, S, P, d, Sf}, {X_ss, S_ss, P_ss, d0, S_f}));
I_Sf = double(subs(I_Sf_sym, {X, S, P, d, Sf}, {X_ss, S_ss, P_ss, d0, S_f}));

fprintf('\nLinearization done.\n');
disp('Jacobian J ='); disp(J);
disp('Input vector I_d (D input) ='); disp(I_d);
disp('Input vector I_Sf (S_f input) ='); disp(I_Sf);
disp('Eigenvalues of J:'); disp(eig(J));

%% TRANSFER FUNCTIONS
% State-space with two inputs: [D, S_f]
A = J;
B = [I_d I_Sf];
C = eye(3);
D_mat = zeros(3,2);

sys = ss(A,B,C,D_mat);
tf_model = tf(sys);
tf_model_min = minreal(tf_model, 1e-6);

tf_model_min = minreal(tf_model, 1e-6);
state_names = {'X','S','P'};  % define before first use

disp('Transfer function model (both inputs, clean display):');
for i = 1:size(tf_model_min,1)
    for j = 1:size(tf_model_min,2)
        fprintf('\nFrom Input %d to Output %s:\n', j, state_names{i});
        display(tf_model_min(i,j)); % clean MATLAB-style output
    end
end


%%  Transfer Functions from D input 
[num_D, den_D] = ss2tf(A, I_d, C, zeros(3,1));
state_names = {'X','S','P'};

for i = 1:3
    fprintf('\nTransfer Function from D to %s:\n', state_names{i});
    G_i = tf(num_D(i,:), den_D);
    disp(G_i)
    
    p = pole(G_i);
    z = zero(G_i);
    fprintf('Poles:\n'); disp(p);
    fprintf('Zeros:\n'); disp(z);
    if all(real(p) < 0)
        fprintf('System is STABLE (all poles have negative real parts)\n');
    elseif any(real(p) > 0)
        fprintf('System is UNSTABLE (some poles have positive real parts)\n');
    else
        fprintf('System is MARGINALLY STABLE\n');
    end
end

%% Transfer Functions from S_f input
[num_Sf, den_Sf] = ss2tf(A, I_Sf, C, zeros(3,1));

for i = 1:3
    fprintf('\nTransfer Function from S_f to %s:\n', state_names{i});
    G_i = tf(num_Sf(i,:), den_Sf);
    disp(G_i)

    p = pole(G_i);
    z = zero(G_i);
    fprintf('Poles:\n'); disp(p);
    fprintf('Zeros:\n'); disp(z);
    if all(real(p) < 0)
        fprintf('System is STABLE (all poles have negative real parts)\n');
    elseif any(real(p) > 0)
        fprintf('System is UNSTABLE (some poles have positive real parts)\n');
    else
        fprintf('System is MARGINALLY STABLE\n');
    end
end

%% STEP RESPONSE & NOISE for D
t_points = linspace(0, 400, 2000);
D_step = d0*ones(size(t_points));

% Step changes in D
D_step(t_points>=50 & t_points<150)  = d0 + 0.05;
D_step(t_points>=150 & t_points<250) = d0 - 0.03;
D_step(t_points>=250 & t_points<350) = d0 + 0.1;

% Linear response for D input
x_d = lsim(sys(:,1), D_step'-d0, t_points);
x_lin = x_d + repmat(x_ss', length(t_points), 1);

% Nonlinear response
[t_non_lin, x_non_lin] = ode45(@(t,x) func(t,x,interp1(t_points,D_step,t),S_f), t_points, x_ss);

% Add Gaussian noise (1% std)
noise_std = [0.01*X_ss, 0.01*S_ss, 0.01*P_ss];
rng(0);
x_lin_noisy     = x_lin     + randn(size(x_lin))     .* noise_std;
x_non_lin_noisy = x_non_lin + randn(size(x_non_lin)) .* noise_std;

fprintf('\nNoise added to simulated outputs:\n');
fprintf('X: %.4f, S: %.4f, P: %.4f\n', noise_std(1), noise_std(2), noise_std(3));

%% PLOTTING
figure('Name','Bioreactor Step Responses with Noise','NumberTitle','off');
states = {'D','X','S','P'};
y_labels = {'D (1/h)','X (g/L)','S (g/L)','P (g/L)'};
titles = {'Step Changes in Dilution Rate','Biomass Response (Noisy)',...
          'Substrate Response (Noisy)','Product Response (Noisy)'};

subplot(4,1,1);
plot(t_points, D_step, 'm', 'LineWidth', 2);
xlabel('Time (min)'); ylabel(y_labels{1});
title(titles{1}); grid on; box on;

for i = 1:3
    subplot(4,1,i+1);
    plot(t_points, x_lin_noisy(:,i), 'r-', 'LineWidth', 1.8); hold on;
    plot(t_non_lin, x_non_lin_noisy(:,i), 'k-', 'LineWidth', 1.5);
    yline(x_ss(i), ':b', 'LineWidth', 1.2);
    ylabel(y_labels{i+1});
    title(titles{i+1});
    if i == 3, xlabel('Time (min)'); end
    legend('Linear (Noisy)','Nonlinear (Noisy)','Steady-State','Location','best');
    grid on; box on;
end

%% STEP RESPONSE FOR S_f INPUT
t_points = linspace(0, 400, 2000);
Sf_step = S_f * ones(size(t_points));  % base value

% Step changes in S_f (feed substrate)
Sf_step(t_points>=50 & t_points<150)  = S_f + 5;   % increase feed
Sf_step(t_points>=150 & t_points<250) = S_f - 3;   % decrease feed
Sf_step(t_points>=250 & t_points<350) = S_f + 8;   % larger increase

% Linear response for S_f input (2nd column of sys)
x_Sf = lsim(sys(:,2), Sf_step' - S_f, t_points);
x_lin_Sf = x_Sf + repmat(x_ss', length(t_points), 1);

% Nonlinear response (simulate full model)
[t_non_lin_Sf, x_non_lin_Sf] = ode45(@(t,x) func(t,x,d0,interp1(t_points, Sf_step, t)), ...
                                     t_points, x_ss);

% Add Gaussian noise (same as before)
x_lin_Sf_noisy     = x_lin_Sf     + randn(size(x_lin_Sf))     .* noise_std;
x_non_lin_Sf_noisy = x_non_lin_Sf + randn(size(x_non_lin_Sf)) .* noise_std;

fprintf('\nNoise added to S_f step response simulation.\n');

%% PLOTTING 
figure('Name','Bioreactor Step Responses (Feed Substrate S_f)','NumberTitle','off');
states = {'S_f','X','S','P'};
y_labels = {'S_f (g/L)','X (g/L)','S (g/L)','P (g/L)'};
titles = {'Step Changes in Feed Substrate Concentration','Biomass Response (Noisy)',...
          'Substrate Response (Noisy)','Product Response (Noisy)'};

subplot(4,1,1);
plot(t_points, Sf_step, 'm', 'LineWidth', 2);
xlabel('Time (min)'); ylabel(y_labels{1});
title(titles{1}); grid on; box on;

for i = 1:3
    subplot(4,1,i+1);
    plot(t_points, x_lin_Sf_noisy(:,i), 'r-', 'LineWidth', 1.8); hold on;
    plot(t_non_lin_Sf, x_non_lin_Sf_noisy(:,i), 'k-', 'LineWidth', 1.5);
    yline(x_ss(i), ':b', 'LineWidth', 1.2);
    ylabel(y_labels{i+1});
    title(titles{i+1});
    if i == 3, xlabel('Time (min)'); end
    legend('Linear (Noisy)','Nonlinear (Noisy)','Steady-State','Location','best');
    grid on; box on;
end

set(gcf,'Position',[100 100 850 900]);
sgtitle('Bioreactor: Linear vs Nonlinear Step Responses for Feed Substrate (S_f)','FontSize',13,'FontWeight','bold');


set(gcf,'Position',[100 100 850 900]);
sgtitle('Bioreactor: Linear vs Nonlinear Step Responses','FontSize',13,'FontWeight','bold');


save('bio_reactor_data.mat', 't_points', 'D_step', 'x_non_lin_noisy', 'x_ss', 'd0');
disp('Bioreactor data saved to bioreactor_data.mat');

