function stability()
    clc; clear; close all;

    % Parameters
    Yxs = 0.4;   % g/g
    beta  = 0.2; % h^(-1)
    Pm  = 50;    % g/L
    Ki  = 22;    % g/L
    alpha = 2.2; % g/
    mum   = 0.48;  % h^(-1)
    Km    = 1.2;   % g/L
    Sf0   = 20;    % g/L
    D0    = 0.202; % h^(-1)

    % Nonlinear model
    mu = @(S,P) (mum * (1 - P / Pm) .* S) ./ (Km + S + (S.^2) / Ki);

    f  = @(x,u,Sf_par) [ ...
        -u * x(1) + mu(x(2), x(3)) * x(1);             % dX/dt
         u * (Sf_par - x(2)) - (1 / Yxs) * mu(x(2), x(3)) * x(1); % dS/dt
        -u * x(3) + (alpha * mu(x(2), x(3)) + beta) * x(1)  % dP/dt
    ];

    % Steady state 
    x0_guess = [0, 2, 0];
    optsFS = optimoptions("fsolve", ...
        "Display","off", ...
        "FunctionTolerance",1e-12, ...
        "StepTolerance",1e-12);

    x_ss = fsolve(@(x) f(x, D0, Sf0), x0_guess, optsFS);
    X_ss = x_ss(1); 
    S_ss = x_ss(2); 
    P_ss = x_ss(3);

    fprintf('Steady state at D0=%.6f, Sf=%.3f -> X=%.6f, S=%.6f, P=%.6f\n', ...
        D0, Sf0, X_ss, S_ss, P_ss);

    % Linearisation (Jacobian) 
    Aden   = Km + S_ss + (S_ss.^2 / Ki);
    dA_dS  = 1 + 2 * S_ss / Ki;
    dG_dS  = (Aden - S_ss * dA_dS) / Aden^2;          % d(S/A)/dS
    dmu_dS = mum * (1 - P_ss / Pm) * dG_dS;
    dmu_dP = -(mum / Pm) * (S_ss / Aden);
    mu_ss  = mum * (1 - P_ss / Pm) * S_ss / Aden;

    M = zeros(3,3);
    % f1
    M(1,1) = -D0 + mu_ss;
    M(1,2) = X_ss * dmu_dS;
    M(1,3) = X_ss * dmu_dP;
    % f2
    M(2,1) = -(1 / Yxs) * mu_ss;
    M(2,2) = -D0 - (1 / Yxs) * X_ss * dmu_dS;
    M(2,3) = -(1 / Yxs) * X_ss * dmu_dP;
    % f3
    M(3,1) = alpha * mu_ss + beta;
    M(3,2) = alpha * X_ss * dmu_dS;
    M(3,3) = -D0 + alpha * X_ss * dmu_dP;

    B    = [-X_ss; Sf0 - S_ss; -P_ss];
    C    = eye(3);
    Dmat = zeros(3,1);

    sys_lin  = ss(M,B,C,Dmat);
    sys_lin  = minreal(sys_lin, 1e-7);

    % Choosing controlled variable 
    fprintf('\nSelect controlled variable (CV) for stability analysis:\n');
    fprintf('  1 = Biomass X\n');
    fprintf('  2 = Substrate S\n');
    fprintf('  3 = Product P\n');
    CVchoice = input('Choice (1/2/3): ');

    switch CVchoice
        case 1
            y_index  = 1;
            cv_name  = 'X';
            cv_title = 'Biomass X';
        case 2
            y_index  = 2;
            cv_name  = 'S';
            cv_title = 'Substrate S';
        case 3
            y_index  = 3;
            cv_name  = 'P';
            cv_title = 'Product P';
        otherwise
            error('Invalid choice. Must be 1, 2, or 3.');
    end

    % SISO transfer function G(s): D to chosen CV
    G_siso = ss(sys_lin(y_index,:));
    G_siso = minreal(G_siso, 1e-7);

   % Closed-loop stability vs Kc 
    Kc_min = -100;
    Kc_max =  100;
    nK     = 4000;
    
    Kc_vec      = linspace(Kc_min, Kc_max, nK);
    maxRealPole = nan(size(Kc_vec));
    
    for i = 1:nK
        Kc    = Kc_vec(i);
        sys_cl = feedback(Kc * G_siso, 1);   % unit(1) negative feedback
        p      = pole(sys_cl);
        maxRealPole(i) = max(real(p));
    end
    
    % Stable if all poles are in LHP
    isStable = (maxRealPole < 0);
    fprintf('\n Closed-loop stability (P controller, unity negative feedback)\n');
    
    if ~any(isStable)
        fprintf('No stable Kc found in the scan range [%.2f, %.2f].\n', Kc_min, Kc_max);
    else
        % Finding stable intervals
        edges    = diff([false, isStable, false]);   % transitions
        startIdx = find(edges == 1);
        endIdx   = find(edges == -1) - 1;
    
        fprintf('Approximate stable Kc interval(s) (from numeric pole locations):\n');
        for j = 1:numel(startIdx)
            Kc_start = Kc_vec(startIdx(j));
            Kc_end   = Kc_vec(endIdx(j));
            fprintf('   Kc in [%.4f, %.4f]\n', Kc_start, Kc_end);
        end
    end

    % Root locus plot
    figure('Name', sprintf('Root locus for D -> %s', cv_name));
    rlocus(G_siso); grid on;
    title(sprintf('Root Locus of D (%s) with P controller ', cv_title));

    fprintf('\nRoot locus plotted for G(s) = D (%s).\n', cv_name);
end
