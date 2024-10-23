% =========================================================================
% TA Session 1.
% Econometrics II with Michael Creel, IDEA Spring 2024
% TA: Conghan Zheng (conghan.zheng@uab.cat)
% Last update: 21 Apr 2024
% Compiled with: ver. R2024a
%
% DESCRIPTION: 
%   Tutorials on PS1 exercises: basics of numerical optimization.
%
% CONTENTS: 
% - 0 Preliminaries
% - 1 Data
% - 2 Unconstrained Minimization
% - 3 Constrained Minimization
% - 4 Jacobian and Hessian
% - * User-defined functions
% 
% INPUT:
%   - TA1.csv (Data on 1978 cars)
% =========================================================================

% 0 Preliminaries ====
clear; clc; 

% - Get the current script's directory
currentscript = matlab.desktop.editor.getActive; 
% - Set it as the working directory
cd(fileparts(currentscript.Filename));
% - Also save the path for later use
scriptPath = fileparts(which(matlab.desktop.editor.getActiveFilename));

% 1 Data ====

% Import data
% raw_data = readtable('TA1.csv');

url = 'https://vincentarelbundock.github.io/Rdatasets/csv/causaldata/auto.csv';
filename = 'auto.csv';
websave(filename, url); % Save the file temporarily
raw_data = readtable(filename);
delete(filename); % Delete the temporary file

% Some cleansing
% - Drop rows with at least one NaN observation
data = rmmissing(raw_data);

% 2 Unconstrained Minimization ====
fprintf('==== Unconstrained Minimization ==== \n');

% OLS =====
% - Specification: log(price) = b0 + b1*mileage + b2*headroom + b3*foreign + u

% - Dependent variable
data.logprice = log(data.price);
y = data.logprice;

% - A column of ones for the constant term
const = ones(size(y,1),1);

% - Regressors
x = [const, data.mpg, data.headroom, data.foreign];

% - Number of observations
n = size(y, 1);

% - Number of regressors
k = size(x, 2);

% - OLS estimates using the analytical formula
beta_analytical = inv(x'*x)*x'*y;

% - OLS estimates using unconstrained minimization, using the previous results 
% as starting values
beta0 = zeros(size(x, 2),1);
beta_numopt = fminunc(@(beta) SSR(beta,x,y), beta0, optimoptions(@fminunc,'Display', 'off'));

% - Post estimation
% -- Calculate residuals
e = y-x*beta_numopt; 
% -- Sample variance of the residuals: E(e'e) = sigma^2*(n-k)
sigma2 = (e'*e)/(n-k);
% -- VAR-COV matrix of the estimated beta: Var(beta_hat) = sigma^2*(x'x)^-1
cov_beta = sigma2*inv(x'*x); 
% -- Standard errors
se = sqrt(diag(cov_beta));
% -- t Statistic (H0: beta = 0, T = beta_est/sigma_est)
t = beta_numopt./se;
% -- p Value based on t (t-test is two-tailed)
p = 2*(1-tcdf(abs(t),n-k));

% Compare results =====
format short % Set output display format: short, 4 decimal places

fprintf('\n-- The OLS estimates from the analytical formula: -- \n');
disp(beta_analytical)

fprintf('\n-- The OLS regression table from unconstrained minimization: --\n');
disp(table(beta_numopt,se,t,p,'VariableNames',{'Estimate','SE','sStat','pValue'}, ...
    'RowNames',{'(Intercept)','mpg','headroom','foreign'}))

fprintf("\n-- The OLS regression table from the built-in 'fitlm' command: --\n");
disp(fitlm(data,'logprice~mpg+headroom+foreign'))

% 3 Constrained Minimization ====
fprintf('\n ==== Constrained Minimization ==== \n');

% OLS with one linear equality constraint =====
% log(price) = b0 + b1*mileage + b2*headroom + b3*foreign + u
%   s.t.: b1 = b2
%
% b1 = b2 <=> b1 - b2 = 0 <=> (0,1,-1,0)*(b0,b1,b2,b3)' = b1 - b2 = 0,
% Aeq = (0,1,-1), beq = 0

% Linear inequality constraint: Aeq*b = beq
Aeq = [0 1 -1 0];
beq = 0;

% Constrained minimization:
% - syntax: beta = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub)
% - In Matlab, fminunc and fmincon have a loose default tolerance, so you might 
%   want to use

%   to tighten the tolerances, in order to be able to correctly replicate the 
%   results.

% - 'SSR': defined at the end of this script
beta_ctr = fmincon(@(beta) SSR(beta,x,y),beta0,[],[],Aeq,beq,[],[],[],optimoptions(@fmincon,'MaxFunEvals', 1e6, 'TolFun',1e-10, 'Display', 'off'));

% - Post estimation
% -- Calculate residuals
e_ctr = y-x*beta_ctr; 
% -- Sample variance of the residuals
sigma2_ctr = (e_ctr'*e_ctr)/(n-k);
% -- Variance-covariance matrix of the estimates
cov_beta_ctr = sigma2_ctr*inv(x'*x); 
% -- Standard errors
se_ctr = sqrt(diag(cov_beta_ctr));
% -- t Statistic
t_ctr = beta_ctr./se_ctr;
% -- p Value
p_ctr = 2*(1-tcdf(abs(t_ctr),n-k));

% Compare results =====
format short

fprintf('\n -- The OLS regression table from numerical optimization: --\n');
disp(table(beta_ctr,se_ctr,t_ctr,p_ctr,'VariableNames',{'Estimate','SE','tStat', ...
    'pValue'}, 'RowNames',{'(Intercept)','mpg','headroom','foreign'}))

fprintf("\n -- The OLS estimates from the built-in 'lsqlin' command: --\n");
disp(lsqlin(x,y,[],[],Aeq,beq,[],[],[],optimset('Display', 'off')))

% 4 Jacobian and Hessian ====
fprintf('\n==== Jacobian and Hessian ==== \n');

% Context: Nonlinear objective function

% 4.1 Jacobian =====

fprintf('\n---- Jacobian ----\n');

% 4.1.1 Finite differences ======

Jacobian_fd0 = grad(beta0, @(beta) SSR(beta, x, y)); 
Jacobian_fdhat = grad(beta_analytical, @(beta) SSR(beta, x, y)); 

% 4.1.2 Symbolic differentiation ======

beta = sym('beta',[k 1],'real');

SSR_val = SSR(beta, x, y);

Jacobian_val = jacobian(SSR_val,beta);

Jacobian_sd0 = subs(Jacobian_val, beta, beta0);
Jacobian_sdhat = subs(Jacobian_val, beta, beta_analytical);

% 4.1.3 (Optional) Automatic differentiation (AD, more precise) ======
% ... evaluates derivatives at certain numeric values (numeric 'double' input
% and numeric 'double' output);
% ... does not construct symbolic expressions for derivatives.

% For finite differences, the solver evaluates k finite differences (see the 
% loop above in 4.1.1).

% While AD reduces the number of function evaluations the solver does, and it 
% does not need to do finite difference steps, so the derivative estimation 
% process takes fewer function evaluations and is (not always faster, but at 
% least) more accurate.

% u = u_n(u_n-1(...u_2(u_1(b))))
% du/db = (du_n/du_n-1)*...*(du_2/du_1)*(du_1(b)/db)

beta = sym('beta',[k 1],'real');
u1_val = y-x*beta;
u1 = matlabFunction(u1_val, 'vars', {beta}); % e: n*1 vector
e = sym('e',[n 1],'real');
u2_val = e'*e; % sigma^2: scalar
u2 = matlabFunction(u2_val, 'vars', {e});

u = {u1, u2};
u_val = {u1_val, u2_val};
u_arg = {beta, e};

Jacobian_ad0 = double(AD(u, u_val, u_arg, beta0));
Jacobian_adhat = double(AD(u, u_val, u_arg, beta_analytical));

% Compare results =====
format short

fprintf('\n-- Jacobian at [0,0,0,0] using different methods: --\n');
disp(table(beta0, double(Jacobian_sd0'), Jacobian_fd0, Jacobian_ad0, ...
    'VariableNames',{'b0','df/db|b=b0,SD','df/db|b=b0,FD','df/db|b=b0,AD'}, ...
    'RowNames',{'(Intercept)','mpg','headroom','foreign'}))

fprintf('\n-- Jacobian at beta_ols using different methods: --\n');
disp(table(beta_analytical, double(Jacobian_sdhat'), Jacobian_fdhat, ...
    Jacobian_adhat, 'VariableNames',{'b_ols', 'df/db|b=b_ols,SD', ...
    'df/db|b=b_ols,FD','df/db|b=b_ols,AD'}, ...
    'RowNames',{'(Intercept)','mpg','headroom','foreign'}))

% 4.2 Hessian =====
fprintf('\n---- Hessian ---- \n');

% Again, the OLS problem
y = data.logprice;
x = [const, data.mpg, data.headroom, data.foreign];
beta = sym('beta',[k 1],'real');
SSR_val = SSR(beta,x,y);

% Jacobian of the gradients
fprintf('\n-- The Jacobian of the Jacobian: -- \n');
disp(jacobian(jacobian(SSR_val,beta),beta))

% Use the built-in Hessian command
fprintf('\n-- The Hessian: -- \n');
disp(hessian(SSR_val, beta))

% APPENDIX =====

% Alternative Community toolbox: DERIVEST
% To compute numeric derivatives, you can use the package:
% https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation 
% Usage: 
% function [grad,err,finaldelta] = gradest(fun,x0)
% function [hess,err] = hessian(fun,x0)
% function [jac,err] = jacobianest(fun,x0)

% Local Functions ====
% Note: Before Ver. R2024a, local functions in scripts must be defined at the 
%       end of the file, after the last line of script code.

% - Function: Compute the SSR
function f = SSR(beta,x,y)
    % SSR of model y = beta*x + u
    f = (y-x*beta)'*(y-x*beta);
end

% Function: Compute the numerical gradient of a function at a vector
function [f] = grad(x, fct)
% INPUT:
%   - x: point at which the gradient is to be calculated, k*1
%   - fct: function to compute the gradient, scalar-valued
% OUTPUT:
%   - f: numeric gradient of the function 'fct' at point 'x', k*1

tol = 10^(-6);
k = length(x);
f = zeros(k,1);

    for i=1:k % loop until we get the first partial derivatives of 'fct' w.r.t. 
              % all the elements of x
        
        % delta: a very small change (only) in the i-th component of x, k*1
        delta=zeros(k,1);
        delta(i)=tol;

        % w: the symmetric difference quotient / symmetric derivative, scalar
        w=(fct(x+delta)-fct(x-delta))/(2*tol);
        
        % Store the results in the i-th element
        f(i)=w;
        
    end

end

% Function: Evaluate the derivatives using the automatic differentiation
%           method, forward mode (bottom up)
function [d] = AD(u, u_val, u_arg, b)
    % Inputs:
    %   u (cell array): List of functions (u_1, u_2, ..., u_n)).
    %   u_val (cell array): List of symbolic expressions of the above functions
    %   u_val (cell array): List of the input arguments of the above functions
    %   b (k*1 vector)
    % Outputs:
    %   d = derivative du/db by automatic differentiation in forward mode
    fval = b;
    for i = 1:length(u)
        % Evaluate derivative of u_n at u_n-1
        dval = subs(jacobian(u_val{i}, u_arg{i}), u_arg{i}, fval); 
        % Evaluate u_n at u_n-1
        fval = u{i}(fval);
        % Update
        if i == 1
            d = dval;
        else
            d = d'*dval';
        end
    end
end