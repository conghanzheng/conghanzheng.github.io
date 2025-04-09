% =========================================================================
% TA Session 3.
% Econometrics II with Michael Creel, IDEA Spring 2024
% TA: Conghan Zheng (conghan.zheng@uab.cat)
% Last update: 14 May 2024
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

% MLE: Bernoulli ====

fprintf('\n ==== MLE: Bernoulli ==== \n');

% Monte Carlo Simulation

n = 100;
p = 0.5;
R = 2000; % number of Monte Carlo replications
b = zeros(R, 1);
contribs = zeros(R, 1);

for r = 1:R
    y = rand(n,1) < p; % a sample of size N from distribution: Bernoulli (k=1,p=0.5)
    phat = mean(y); % estimate p
    b(r,:) = phat;
    contribs(r,:) = sqrt(n)*(phat-p); % estimate sqrt{n}*(p-p^0)
end

b0 = p;
fprintf('True asymptotic mean: %.1f\n', p);
bmc = mean(b);
fprintf('Monte Carlo mean: %.4f\n', bmc);
s0 = sqrt(p-p^2);
fprintf('True asymptotic standard error: %.1f\n', s0);
smc = sqrt(var(contribs));
fprintf('Monte Carlo standard error: %.4f\n', smc);

% GMM: Dynamic CAPM ====

fprintf('\n ==== GMM: Dynamic CAPM ==== \n');

url = 'https://raw.githubusercontent.com/mcreel/Econometrics/main/Examples/GMM/hall.m';
filename = 'hall.m';
outfilename = websave(filename, url);
hall = load(outfilename);
delete(outfilename);

n = size(hall, 1);

c = hall(3:n,1);
c1 = hall(2:n-1,1);
c2 = hall(1:n-2,1);

r = hall(3:n,2);
r1 = hall(2:n-1,2);
r2 = hall(1:n-2,2);

n = size(c,1); % update sample size

e = @(theta) theta(1,:)*r.*c.^(theta(2,:) - 1) - 1; % the Euler eq. (FOC)
inst = [ones(n,1) c1 r1 c2 r2]; % instruments
W = eye(size(inst,2)); % initial weight matrix

% moments
ms = @(theta) e(theta) .* inst; % moment contributions (N*L)
%m = @(theta) mean(ms(theta),1)'; % moment conditions (L*1)
m = @(theta) (1/n)*inst'*e(theta); % moment conditions (L*1)

theta0 = [0;0]; % starting values: discount and utility
opts_fminunc = optimoptions('fminunc', 'Display', 'off', 'TolX', 1e-15, 'TolFun', 1e-12, 'MaxFunEvals', 1e6);

% 1-step
[thetahat_1s, obj1] = fminunc(@(theta) n*m(theta)'*W*m(theta), theta0, opts_fminunc);  % note the n in there, obj function is scaled to converge to chi^2.

% standard errors (if we assume all moment contributions are independent)
Dtheta1 = r.*c.^(thetahat_1s(2,:)-1);
Dtheta2 = Dtheta1.*log(c)*thetahat_1s(1,:);
D = [mean(Dtheta1.*inst); mean(Dtheta2.*inst)];

vc = inv(D*W*D')/n;
se_1s = sqrt(diag(vc));

% standard errors (if we don't assume all moment contributions are independent)
vc = inv(D*inv(NeweyWest(ms(thetahat_1s), floor(size(e(thetahat_1s),1)^0.25)))*D')/n;
se_1s_nw = sqrt(diag(vc));

fprintf('\n-- One-Step --\n')
fprintf('theta_hat: (%.4f,%.4f) \n', thetahat_1s);
fprintf('se: (%.4f,%.4f) \n', se_1s);
fprintf('se(Newey–West): (%.4f,%.4f) \n', se_1s_nw);

% efficient weight matrix
mm = ms(thetahat_1s); % moment contributions eval at the 1s GMM est: (N,L)
omega = mm'*mm/n; % var-cov: (L,N)*(N,L) = (L,L)
% omega = cov(mm); % the same result
Wo = inv(omega);

% 2-step
[thetahat_2s, obj2] = fminunc(@(theta) n*m(theta)'*Wo*m(theta), thetahat_1s, opts_fminunc);

% standard errors (if we assume all moment contributions are independent)
Dtheta1 = r.*c.^(thetahat_2s(2,:)-1);
Dtheta2 = Dtheta1.*log(c)*thetahat_2s(1,:);
D = [mean(Dtheta1.*inst); mean(Dtheta2.*inst)];

vc = inv(D*Wo*D')/n;
se_2s = sqrt(diag(vc));

% standard errors (if we don't assume all moment contributions are independent)
vc = inv(D*inv(NeweyWest(ms(thetahat_2s), floor(size(e(thetahat_2s),1)^0.25)))*D')/n;
se_2s_nw = sqrt(diag(vc));

fprintf('\n-- Two-Step --\n');
fprintf('theta_hat: (%.4f,%.4f) \n', thetahat_2s);
fprintf('se: (%.4f,%.4f) \n', se_2s);
fprintf('se(Newey–West): (%.4f,%.4f) \n', se_2s_nw);

% iterated GMM (optional for your solution)
fprintf('\n -- Iterated GMM -- \n');

tol = 1e-6; % convergence tolerance
max_iter = 100; % maximum number of iterations

theta_old = thetahat_2s; % initialize theta_old with the two-step estimator
W_old = Wo; % initialize W_old with the efficient weight matrix from the two-step estimator

for iter = 1:max_iter
    
    % calculate weight matrix
    mm = ms(theta_old);
    omega = mm' * mm / n; 
    W_new = inv(omega);
    
    % re-estimate theta using the new weight matrix
    [theta_new, ~] = fminunc(@(theta) n * m(theta)' * W_new * m(theta), theta_old, opts_fminunc);
    
    % check for convergence
    if max(abs(theta_new - theta_old)) < tol && max(abs(W_new - W_old), [], 'all') < tol
        
        % standard errors (if we assume all moment contributions are independent)
        Dtheta1 = r.*c.^(theta_new(2,:)-1);
        Dtheta2 = Dtheta1.*log(c)*theta_new(1,:);
        D = [mean(Dtheta1.*inst); mean(Dtheta2.*inst)];
        
        vc = inv(D*W_new*D')/n;
        se_ir = sqrt(diag(vc));
        
        % standard errors (if we don't assume all moment contributions are independent)
        vc = inv(D*inv(NeweyWest(ms(theta_new), floor(size(e(theta_new),1)^0.25)))*D')/n;
        se_ir_nw = sqrt(diag(vc));

        fprintf('Convergence achieved after %d iterations.\n', iter);
        fprintf('theta_hat: (%.4f,%.4f) \n', theta_new);
        fprintf('se: (%.4f,%.4f) \n', se_ir);
        fprintf('se(Newey–West): (%.4f,%.4f) \n', se_ir_nw);
        break;
    end
    
    % update theta and W
    theta_old = theta_new;
    W_old = W_new;
   
end

if iter == max_iter
    fprintf('\nMaximum number of iterations reached without convergence.\n');
end

% Local Functions ====

% Function: the Newey-West estimator of the asymptotic variance matrix
% link: https://raw.githubusercontent.com/mcreel/Econometrics/main/Examples/GMM/NeweyWest.m
function omegahat = NeweyWest(Z,nlags)

    % INPUTS: Z - the object
    %         nlags - number of lags
    %
    % OUTPUTS: omegahat - Cov
    
    [n,k] = size(Z);
    
    % de-mean
    Z = Z - repmat(mean(Z),n,1);
    
    omegahat = Z'*Z/n; % sample variance
    if nlags > 0
       % sample autocovariances
       for i = 1:nlags
          Zlag = Z(1:n-i,:);
          ZZ = Z(i+1:n,:);
          gamma = (ZZ'*Zlag)/n;
          weight = 1 - (i/(nlags+1));
          omegahat = omegahat + weight*(gamma + gamma');
       end
    end
end