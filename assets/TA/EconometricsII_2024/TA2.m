% =========================================================================
% TA Session 2.
% Econometrics II with Michael Creel, IDEA Spring 2024
% TA: Conghan Zheng (conghan.zheng@uab.cat)
% Last update: 26 Apr 2024
% Compiled with: ver. R2024a
%
% DESCRIPTION: 
%   Some sample codes you can use in your solution.
% =========================================================================

% 0 Preliminaries ====
clear; clc; 

% - Get the current script's directory
currentscript = matlab.desktop.editor.getActive;
% - Set it as the working directory
cd(fileparts(currentscript.Filename));
% - Also save the path for later use
scriptPath = fileparts(which(matlab.desktop.editor.getActiveFilename));

% 1 Simulation ====

fprintf('\n ==== 1 SIMULATION ==== \n');

fprintf('\n-- Pseudorandom-number generators -- \n');

% Draws from Normal_ll
normal_sample_tmp1 = round(normrnd(0, 1, 1, 6), 4);
normal_sample_tmp2 = round(normrnd(0, 1, 1, 6), 4);

fprintf('\n Sample 1:\n');
disp(normal_sample_tmp1);

fprintf('\n Sample 2 from the same distribution:\n');
disp(normal_sample_tmp2);

% Seed =====

fprintf('\n-- Using Seed -- \n');

% Set the seed of the random number generator and take two samples
rng(54321); % Set the seed

% First set of samples
normal_sample_tmp3 = round(normrnd(0, 1, 1, 6), 4);
normal_sample_tmp4 = round(normrnd(0, 1, 1, 6), 4);

fprintf('\n Sample 3:\n');
disp(normal_sample_tmp3);

fprintf('\n Sample 4:\n');
disp(normal_sample_tmp4);

fprintf('\n[Seed reset to the same number ...]\n');

% Reset the seed to get the same samples again
rng(54321); % Reset the seed

% Second set of samples (should be the same as the first set)
normal_sample_tmp5 = round(normrnd(0, 1, 1, 6), 4);
normal_sample_tmp6 = round(normrnd(0, 1, 1, 6), 4);

fprintf('\n Sample 5 (same as Sample 3):\n');
disp(normal_sample_tmp5);

fprintf('\n Sample 6 (same as Sample 4):\n');
disp(normal_sample_tmp6);

% Simulation Example: OLS with chi^2 errors =====

fprintf('\n -- Simulation Example: OLS with chi^2 errors -- \n');

n1 = 15;
[b1, se1, t1, p1] = chi2ols(n1);
fprintf("\n [b, se, t, p] from a random sample of size %d = [%.4f, %.4f, %.4f, %.4f] \n", n1, b1, se1, t1, p1);

n2 = 150;
[b2, se2, t2, p2] = chi2ols(n2);
fprintf("\n [b, se, t, p] from a random sample of size %d = [%.4f, %.4f, %.4f, %.4f] \n", n2, b2, se2, t2, p2);

n3 = 1500;
[b3, se3, t3, p3] = chi2ols(n3);
fprintf("\n [b, se, t, p] from a random sample of size %d = [%.4f, %.4f, %.4f, %.4f] \n", n2, b2, se2, t2, p2);

% Grid Search =====

fprintf('\n -- Grid Search -- \n');

fprintf('\n Check the figure window. \n');

chi2ols_grid_combined();

% 2 MLE ====

fprintf('\n ==== 2 MLE ==== \n');

rng(54321);

% Specification:
% y = b0 + b1*x1 + b2*x2 + b3*x3 + u, u ~ N(0,s^2)

k = 4; % number of regressors including the constant term
theta = rand(k+1,1);
theta0 = ones(k+1,1);
Aeq = [0 1 0 0 0; 0 0 1 1 0];
beq = [1;1];
J = size(beq,1);

opts_fminunc = optimoptions(@fminunc,'Display','off');
opts_fmincon = optimoptions(@fmincon,'Display','off');

% Minimize the negative log-likelihoods
% - Note: minimized obj value estimated from fminunc/fmincon is the negative log-likelihood

n1 = 50;
data1 = random_data(n1, k, theta);
[theta_U1, ll_U1] = fminunc(@(theta) -mean(Normal_ll(theta, data1)), theta0, opts_fminunc);
ll_U1 = - ll_U1;
[theta_C1, ll_C1] = fmincon(@(theta) -mean(Normal_ll(theta, data1)), theta0,[],[],Aeq,beq,[],[],[],opts_fmincon);
ll_C1 = - ll_C1;
LR1 = 2*n1*(ll_U1 - ll_C1); 
p_LR1 = 1.0 - chi2cdf(LR1, J);
%[h1,p1] = lratiotest(ll_U1*n1,ll_C1*n1,J); % or use the built-in command

n2 = 500;
data2 = random_data(n2, k, theta);
[theta_U2, ll_U2] = fminunc(@(theta) -mean(Normal_ll(theta, data2)), theta0, opts_fminunc);
ll_U2 = - ll_U2;
[theta_C2, ll_C2] = fmincon(@(theta) -mean(Normal_ll(theta, data2)), theta0,[],[],Aeq,beq,[],[],[],opts_fmincon);
ll_C2 = - ll_C2;
LR2 = 2*n2*(ll_U2 - ll_C2);
p_LR2 = 1.0 - chi2cdf(LR2, J);

n3 = 5000;
data3 = random_data(n3, k, theta);
[theta_U3, ll_U3] = fminunc(@(theta) -mean(Normal_ll(theta, data3)), theta0, opts_fminunc);
ll_U3 = - ll_U3;
[theta_C3, ll_C3] = fmincon(@(theta) -mean(Normal_ll(theta, data3)), theta0,[],[],Aeq,beq,[],[],[],opts_fmincon);
ll_C3 = - ll_C3;
LR3 = 2*n3*(ll_U3 - ll_C3);
p_LR3 = 1.0 - chi2cdf(LR3, J);

% MLE (an alternative way): using a built-in command to write the normal density
% of the linear regression model with normal errors
% Syntax: normpdf(y,mu(y),sigma(y)), the pdf of N(mu(y),sigma(y)^2) evaluated at y

Normal_ll_alt3 = @(theta) -sum(log(normpdf(data3(:,1),  data3(:,2:end) * theta(1:end-1), theta(end))));
[theta_U3_alt, ll_U3_alt] = fminunc(Normal_ll_alt3, theta0, opts_fminunc);
ll_U3_alt = - ll_U3_alt;

fprintf('\n -- MLE of the linear regression model with Normal errors, using random data: -- \n\n ');
disp(table(theta, theta_U1, theta_U2, theta_U3, theta_U3_alt, ...
    'VariableNames',{'True Parameter',['N=', num2str(n1)],['N=', num2str(n2)],['N=', num2str(n3)],['N=', num2str(n3),', normalpdf cmd']}, ...
    'RowNames',{'beta_0','beta_1','beta_2','beta_3','sigma'}))
fprintf('\n -- The LR Tests: -- \n \n H0: (1) beta_1 = 1; (2) beta_2 + beta_3 = 1 \n\n');
disp(table([n1,n2,n3]',[ll_U1,ll_U2,ll_U3]',[ll_C1,ll_C2,ll_C3]',[LR1,LR2,LR3]',[p_LR1,p_LR2,p_LR3]',...
    'VariableNames',{'Sample Size','Log-Likelihood (U)','Log-Likelihood (C)', 'LR Stat.', 'p-Value'}));

% Local Functions ====

% Function: estimate the OLS model with Chi-sq errors
function [b, se, t, p] = chi2ols(n)
    rng(37261); % seed
    x = random('Chisquare', 1, n, 1);
    u = random('Chisquare', 1, n, 1) - 1; % demeaned chi2 error
    y = 1 + 2*x + u;
    k = 2;

    b = inv(x'*x)*x'*y;
    e = y-x*b;
    se = sqrt(diag((e'*e)/(n-k)*inv(x'*x)));

    t = b ./ se; % Element-wise division
    p = 2 * (1 - tcdf(abs(t), n - k));
end

% Function: Using Grid search method to estimate the OLS model with Chi-sq 
%           errors
function chi2ols_grid_combined()
    % Define the values for n and G for three samples
    n_values = [100, 500, 1000];
    G_values = [150, 999, 1234];

    % Initialize the figure
    figure;

    % Loop over the different (n, G) combinations
    for i = 1:length(n_values)
        n = n_values(i); % number of observations
        G = G_values(i); % the grid, a subset of the parameter space

        % Generate some random data
        rng(37261); % Set random seed
        x = random('Chisquare', 1, n, 1);
        u = random('Chisquare', 1, n, 1) - 1;
        y = 1 + 2*x + u;

        % Compute the SSR over the grid
        betas = linspace(0, 5, G);
        ssr = zeros(1, G);
        for j = 1:G
            beta = betas(j);
            e = y - 1 - beta * x;
            ssr(j) = sum(e .^ 2);
        end

        % Search for the minimum SSR, return (the index j of) the grid value
        % associated with it
        [~, ind] = min(ssr);

        % Plot the beta-SSR curve
        plot(betas, ssr, 'DisplayName', ['N = ', num2str(n), ', G = ', num2str(G), ', beta* = ', num2str(betas(ind),'%.2f')]);
        hold on;
        scatter(betas(ind), ssr(ind), 'r', 'filled', 'HandleVisibility', 'off');
    end

    % Customize the plot
    xlabel('β');
    ylabel('SSR');
    title('Combined OLS Grid Search, True Beta = 2');
    legend;
    hold off;
end

% Function: generate the random sample [y x] with N, k, beta, and sigma
function [z] = random_data(n, k, theta)

    b = theta(1:k,:);
    s = theta(k+1,:);

    x = [ones(n,1) rand(n,k-1)];
    e = s*randn(n,1);
    y = x*b + e;

	z = [y x];
end

% Function: the log-likelihood for the classical linear model
% y = x*b+e, with e~N(0,s^2)
% the parameter theta = [b' s]
function [logdensity] = Normal_ll(theta, data)

    y = data(:,1);
	x = data(:,2:end);
    [n k] = size(x);

    b = theta(1:k,:);
    s = theta(k+1,:);

	e = y - x*b;
    
    % sqrt(s^2): to make sure the value to put in the log() is positive
	logdensity = -log(sqrt(2*pi)) - log(sqrt(s^2)) - e.*e/(2*s^2);
end