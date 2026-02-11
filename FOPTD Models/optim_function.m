function fun_val = optimfunc_bio(xo, y_index)

    % Loading step response data
    data = load('bioreactor_data.mat');

    t          = data.t_vec(:);
    x_nl_noisy = data.x_nl_noisy;
    x_ss       = data.x_ss(:);
    D_0         = data.D0;
    delta_D     = data.deltaD;

    y_data = x_nl_noisy(:, y_index);
    y_ss   = x_ss(y_index);

    
    D_vec = D_0 * ones(size(t));
    D_vec(t >= 100) = D_0 - 1*delta_D;
    D_vec(t >= 200) = D_0 - 2*delta_D;
    D_vec(t >= 300) = D_0 - 3*delta_D;

    dev = D_vec - D_0;        % deviation input

    K_p    = xo(1);
    Tau   = xo(2);
    Theta = xo(3);

    if Tau <= 0
        fun_val = 1e20;
        return;
    end

    Ts = t(2) - t(1);
    N_d = max(0, round(Theta / Ts));    % samples of delay

    u_delay = zeros(size(dev));
    if N_d < length(dev)
        u_delay(N_d+1:end) = dev(1:end-N_d);
    end

    y_dev = zeros(size(t));
    for k = 1:length(t)-1
        y_dev(k+1) = y_dev(k) + (Ts/Tau)*(-y_dev(k) + K_p*u_delay(k));
    end

    y_model = y_ss + y_dev;

    % SSE 
    err = y_data - y_model;
    fun_val = sum(err.^2);

end
