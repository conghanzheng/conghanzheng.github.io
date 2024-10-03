%==========================================================================
% Exercise 2. System GMM
% Problem Set 2, Microeconometrics Fall 2024
%
% DESCRIPTION: 
%   Performing Arellano-Bover estimation: 
%   y_it = alpha*y_it-1 + eta_i + epsilon_it 
%   using a balanced panel with four periods.
%
% CALLS: user-defined functions used in this script
%   - none 
%==========================================================================

close all; clear; clc;
rng(13);

% Setting working directory to that of the current script (no need to set it manually)
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
scriptPath = fileparts(which(matlab.desktop.editor.getActiveFilename));

% Import data. CSV file exported from Stata
data = readtable('PS1_2.csv');

% Take the dependent variable
h = data.healthy;

% Number of individuals and periods
T = 4;
N = size(h,1)/T;

% Define the vector of regressor: delta_y(i,t-1), delta_y(i,t-2), delta_y(i,t-3)
X = [];
for i=1:N
   j = 1 + (i-1)*T;
   X_ind = [h(i+1) - h(i);...
            h(i+2) - h(i+1);...
            h(i+2)];
   X = [X; X_ind];
end

% Define the vector of dependent variable: y(i,t)
Y = [];
for i=1:N
   j = 1 + (i-1)*T;
   Y_ind = [h(i+2) - h(i+1);...
            h(i+3) - h(i+2);...
            h(i+3)];
   Y = [Y; Y_ind];
end

% Define the matrix of instruments
Z = [];
for i=1:N
    j = 1 + (i-1)*T;
    Z_ind = [h(j) zeros(1,(T-2)*((T-1)/2+1)-1);...
            0 h(j) h(j+1) zeros(1,(T-2)*((T-1)/2+1)-3);...
            0 0 0 h(j+1)-h(j) h(j+2)-h(j+1)];
    Z = [Z; Z_ind];    
end

% Weighting matrix
W = [2 -1 zeros(1,(T-2)*((T-1)/2+1)-2);...
    -1 2 -1 zeros(1,(T-2)*((T-1)/2+1)-3);...
    0 -1 2 -1 zeros(1,(T-2)*((T-1)/2+1)-4);...
    0 0 -1 2 -1;
    0 0 0 -1 2];

% System-GMM estimator
alpha_sysgmm = inv(X'*Z*W*Z'*X)*X'*Z*W*Z'*Y;
disp('System-GMM estimates:')
disp(' ')
disp(alpha_sysgmm);
