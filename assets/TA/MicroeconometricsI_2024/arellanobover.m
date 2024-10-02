function [gamma_abv, e_abv] = arellanobover(y,x)
%==========================================================================
% DESCRIPTION: 
% This function estimates Arellano-Bover model for dynamic panel data with
% three lags. One can adjust it easily by including higher lags in matrix
% Z. To adjust matrix Z for different number of observations for each
% individual it needs to be defined by loop that for each row with complete
% number of lags every iteration i has first sum_{j=1}^{i-1}j zeros and
% then the required number of lags.
%
% INPUT: 
%   - y: Dependent variable
%   - x: Matrix of covariates
%
% OUTPUT:
%   - gamma_ab: coefficients for each x respectively
%   - se_ab: standard errors
%==========================================================================

% distinguish independent variables from matrix x
firm = x(:,1);
n = x(:,2);
w = x(:,3);
k = x(:,4);

% construct technical variables, which will allow to compute differenced
% variables
firm_tech=zeros(length(n),1);
n_tech=zeros(length(n),1);
k_tech=zeros(length(k),1);
w_tech=zeros(length(w),1);
y_tech=zeros(length(y),1);

for i=1:length(n)-1
    firm_tech(i+1)=firm(i);
    n_tech(i+1)=n(i);
    k_tech(i+1)=k(i);
    w_tech(i+1)=w(i);
    y_tech(i+1)=y(i);
end

%computing first diferences, replacing 0 for missing values

n_delta=n-n_tech;
n_delta(firm~=firm_tech)=0;
n_delta(n==0)=0;
n_delta(n_tech==0)=0;

k_delta=k-k_tech;
k_delta(firm~=firm_tech)=0;
k_delta(n==0)=0;
k_delta(n_tech==0)=0;

w_delta=w-w_tech;
w_delta(firm~=firm_tech)=0;
w_delta(n==0)=0;
w_delta(n_tech==0)=0;

y_delta=y-y_tech;
y_delta(firm~=firm_tech)=0;
y_delta(n==0)=0;
y_delta(n_tech==0)=0;


%defining laged differenced n
n_delta_1=zeros(length(n_delta),1);
for i=1:length(n_delta)-1
    n_delta_1(i+1)=n_delta(i);
end

%dropping observations for the lowest t
for i=1:9:length(n_delta_1)-8
    n_delta_1(i)=0;
    n_delta(i)=0;
    k_delta(i)=0;
    y_delta(i)=0;
    w_delta(i)=0;
end

% defining matrix of instruments
% first two rows are necessary since we are using zeros instead of missing
% values

Z_abv=[];
for i=1:9:length(n)-8
    z_abv=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       n(i) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 n(i) n(i+1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 n(i) n(i+1) n(i+2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 n(i) n(i+1) n(i+2) n(i+3) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 n(i) n(i+1) n(i+2) n(i+3)  n(i+4) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n(i) n(i+1) n(i+2) n(i+3) n(i+4) n(i+5) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n(i) n(i+1) n(i+2) n(i+3) n(i+4) n(i+5) n(i+6) 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n_delta(i) n_delta(i+1) n_delta(i+2) n_delta(i+3) n_delta(i+4) n_delta(i+5) n_delta(i+6) n_delta(i+7)];
   Z_abv=[Z_abv;z_abv];
end

%transforming data for arellano bover

n_abv_1=n_delta_1;
n_abv=n_delta;
k_abv=k_delta;
y_abv=y_delta;
w_abv=w_delta;

%substituting observations for the lowest t
for i=1:max(firm)
    n_abv_1(i*9)=n(i*9-1);
    n_abv(i*9)=n(i*9-1);
    k_abv(i*9)=k(i*9-1);
    y_abv(i*9)=y(i*9-1);
    w_abv(i*9)=w(i*9-1);
end

%defining optimal weighting matrix
V=[2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;...
   -1 2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

j = 1;

for i=1:33
   
    v = [zeros(1,j), -1, 2, -1, zeros(1,33-j)];
    V = [V; v];
    j = j+1;
    
end
v1 = [zeros(1,34) , -1, 2];
V = [V; v1] ;

%computing the abv estimator
X_abv=horzcat(n_abv_1, k_abv, w_abv, y_abv);

gamma_abv=(X_abv'*Z_abv*V*Z_abv'*X_abv)^(-1)*X_abv'*Z_abv*V*Z_abv'*n_abv;

%saving the matrix of residuals
for i=1:max(firm)
    e_abv(:,i)=n_abv((i-1)*9+1:i*9,1)-X_abv((i-1)*9+1:i*9,:)*gamma_abv;
end

end