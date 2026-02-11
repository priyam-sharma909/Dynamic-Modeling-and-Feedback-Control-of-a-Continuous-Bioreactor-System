
clear; clc; close all;

% Loading step response data from bioreactor simulation
sim = load("bioreactor_data.mat");

t        = sim.t_vec(:);
x_nl     = sim.x_nl;
x_nl_noisy = sim.x_nl_noisy;
x_ss     = sim.x_ss(:);
deltaD   = sim.deltaD;

Ts = t(2) - t(1);

n_states = 3;              
K_p    = zeros(n_states,1);
Tau_p  = zeros(n_states,1);
Theta = zeros(n_states,1);
SSE       = zeros(n_states,1);
exitflag_all = zeros(n_states,1);

% Initial guess 

Dy_ss_P  = x_nl_noisy(end,3) - x_nl_noisy(1,3);
Kp_0      = Dy_ss_P / (deltaD + eps);
Tau_0     = (t(end)-t(1))/5;
Theta_0   = 0;
x_initial_guess = [Kp_0, Tau_0, Theta_0];

% Bounds on [Kp, Tau, Theta]
Up_Bound = [  1000,  500, 100 ];
Lo_Bound = [ -1000,  0.1,   0 ];

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

state_names = {'X (biomass)','S (substrate)','P (product)'};

for y_index = 1:n_states 
    fprintf('\n IDENTIFICATION FOR STATE %d: (%s)\n', ... 
            y_index, state_names{y_index}); 

    % call fmincon with y_index passed into optimfunc_bio
    [opt_x, fval, exitflag] = fmincon( ... 
        @(x) optimfunc_bio(x, y_index), ... 
        x_initial_guess, [], [], [], [], ...
        Lo_Bound, Up_Bound, [], options); 

    K_p(y_index)    = opt_x(1);
    Tau_p(y_index)  = opt_x(2);
    Theta(y_index) = opt_x(3);
    SSE(y_index)       = fval;
    exitflag_all(y_index) = exitflag;

    fprintf('  K_p    = %.6g\n', opt_x(1));
    fprintf('  Tau_p  = %.6g\n', opt_x(2)); 
    fprintf('  Theta = %.6g\n', opt_x(3)); 
    fprintf('  SSE       = %.6g\n', fval); 
end 

% Save all results
save('bio_foptd_identification_all.mat', ...
     'K_p','Tau_p','Theta','SSE', ...
     'x_ss','t','deltaD','exitflag_all');

disp('FOPTD identification for all three states saved to bio_foptd_identification_all.mat'); %[output:8df67949]


