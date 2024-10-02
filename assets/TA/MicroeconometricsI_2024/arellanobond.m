function [gamma_ab, se_ab, e_ab] = arellanobond(y, x)
%==========================================================================
% DESCRIPTION: 
% This function estimate Arellano-Bond model for dynamic panel data with 
% three lags. One can modify it easily by including higher lags in matrix Z.
%
% INPUT: 
%   - y: Dependent variable
%   - x: Matrix of covariates
%
% OUTPUT:
%   - gamma_ab: coefficient estimates
%   - se_ab: standard errors
%   - e_ab: matrix of residuas, necessary for Sargan test
%==========================================================================

% Distinguish independent variables from matrix x
healthy = x(:,1);
college = x(:,2);
black = x(:,3);
taxincome = x(:,4);
children03 = x(:,5);
children412 = x(:,6);
children1318 = x(:,7);
birthdum1963 = x(:,8);
birthdum1964 = x(:,9);
birthdum1965 = x(:,10);
birthdum1966 = x(:,11);
birthdum1967 = x(:,12);
birthdum1968 = x(:,13);
birthdum1969 = x(:,14);
birthdum1970 = x(:,15);
birthdum1971 = x(:,16);
birthdum1972 = x(:,17);
birthdum1973 = x(:,18);
birthdum1974 = x(:,19);
birthdum1975 = x(:,20);
birthdum1976 = x(:,21);
birthdum1977 = x(:,22);
i_agecat_27 = x(:,37);
i_agecat_32 = x(:,38);
i_agecat_37 = x(:,39);
i_agecat_42 = x(:,40);
i_agecat_47 = x(:,41);
i_agecat_52 = x(:,42);
i_agecat_57 = x(:,43);
i_agecat_62 = x(:,44);
i_ageXmarr_22 = x(:,45);
i_ageXmarr_27 = x(:,46);
i_ageXmarr_32 = x(:,47);
i_ageXmarr_37 = x(:,48);
i_ageXmarr_42 = x(:,49);
i_ageXmarr_47 = x(:,50);
i_ageXmarr_52 = x(:,51);
i_ageXmarr_57 = x(:,52);
i_ageXmarr_62 = x(:,53);

% Construct technical variables, which will allow to compute differenced
% variables
healthy_tech = zeros(length(healthy),1);
college_tech = zeros(length(college),1);
black_tech = zeros(length(black),1);
taxincome_tech = zeros(length(taxincome),1);
children03_tech = zeros(length(children03),1);
children412_tech = zeros(length(children412),1);
children1318_tech = zeros(length(children1318),1);
birthdum1963_tech = zeros(length(birthdum1963),1);
birthdum1964_tech = zeros(length(birthdum1964),1);
birthdum1965_tech = zeros(length(birthdum1965),1);
birthdum1966_tech = zeros(length(birthdum1966),1);
birthdum1967_tech = zeros(length(birthdum1967),1);
birthdum1968_tech = zeros(length(birthdum1968),1);
birthdum1969_tech = zeros(length(birthdum1969),1);
birthdum1970_tech = zeros(length(birthdum1970),1);
birthdum1971_tech = zeros(length(birthdum1971),1);
birthdum1972_tech = zeros(length(birthdum1972),1);
birthdum1973_tech = zeros(length(birthdum1973),1);
birthdum1974_tech = zeros(length(birthdum1974),1);
birthdum1975_tech = zeros(length(birthdum1975),1);
birthdum1976_tech = zeros(length(birthdum1976),1);
birthdum1977_tech = zeros(length(birthdum1977),1);
i_agecat_27_tech = zeros(length(i_agecat_27),1);
i_agecat_32_tech = zeros(length(i_agecat_32),1);
i_agecat_37_tech = zeros(length(i_agecat_37),1);
i_agecat_42_tech = zeros(length(i_agecat_42),1);
i_agecat_47_tech = zeros(length(i_agecat_47),1);
i_agecat_52_tech = zeros(length(i_agecat_52),1);
i_agecat_57_tech = zeros(length(i_agecat_57),1);
i_agecat_62_tech = zeros(length(i_agecat_62),1);
i_ageXmarr_22_tech = zeros(length(i_ageXmarr_22),1);
i_ageXmarr_27_tech = zeros(length(i_ageXmarr_27),1);
i_ageXmarr_32_tech = zeros(length(i_ageXmarr_32),1);
i_ageXmarr_37_tech = zeros(length(i_ageXmarr_37),1);
i_ageXmarr_42_tech = zeros(length(i_ageXmarr_42),1);
i_ageXmarr_47_tech = zeros(length(i_ageXmarr_47),1);
i_ageXmarr_52_tech = zeros(length(i_ageXmarr_52),1);
i_ageXmarr_57_tech = zeros(length(i_ageXmarr_57),1);
i_ageXmarr_62_tech = zeros(length(i_ageXmarr_62),1);

for i=1:length(n)-1
        
    healthy_tech(i+1)=healthy(i);
    college_tech(i+1)=college(i);
    black_tech(i+1)=black(i);
    taxincome_tech(i+1)=taxincome(i);
    children03_tech(i+1)=children03(i);
    children412_tech(i+1)=children412(i);
    children1318_tech(i+1)=children1318(i);
    birthdum1963_tech(i+1)=birthdum1963(i);
    birthdum1964_tech(i+1)=birthdum1964(i);
    birthdum1965_tech(i+1)=birthdum1965(i);
    birthdum1966_tech(i+1)=birthdum1966(i);
    birthdum1967_tech(i+1)=birthdum1967(i);
    birthdum1968_tech(i+1)=birthdum1968(i);
    birthdum1969_tech(i+1)=birthdum1969(i);
    birthdum1970_tech(i+1)=birthdum1970(i);
    birthdum1971_tech(i+1)=birthdum1971(i);
    birthdum1972_tech(i+1)=birthdum1972(i);
    birthdum1973_tech(i+1)=birthdum1973(i);
    birthdum1974_tech(i+1)=birthdum1974(i);
    birthdum1975_tech(i+1)=birthdum1975(i);
    birthdum1976_tech(i+1)=birthdum1976(i);
    birthdum1977_tech(i+1)=birthdum1977(i);
    i_agecat_27_tech(i+1)=i_agecat_27(i);
    i_agecat_32_tech(i+1)=i_agecat_32(i);
    i_agecat_37_tech(i+1)=i_agecat_37(i);
    i_agecat_42_tech(i+1)=i_agecat_42(i);
    i_agecat_47_tech(i+1)=i_agecat_47(i);
    i_agecat_52_tech(i+1)=i_agecat_52(i);
    i_agecat_57_tech(i+1)=i_agecat_57(i); 
    i_agecat_62_tech(i+1)=i_agecat_62(i);
    i_ageXmarr_22_tech(i+1)=i_ageXmarr_22(i);
    i_ageXmarr_27_tech(i+1)=i_ageXmarr_27(i);
    i_ageXmarr_32_tech(i+1)=i_ageXmarr_32(i);
    i_ageXmarr_37_tech(i+1)=i_ageXmarr_37(i);
    i_ageXmarr_42_tech(i+1)=i_ageXmarr_42(i);
    i_ageXmarr_47_tech(i+1)=i_ageXmarr_47(i);
    i_ageXmarr_52_tech(i+1)=i_ageXmarr_52(i);
    i_ageXmarr_57_tech(i+1)=i_ageXmarr_57(i);
    i_ageXmarr_62_tech(i+1)=i_ageXmarr_62(i);
end

% Compute first diferences, replacing 0 for missing values

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

% Define matrix of instruments. The first two rows are necessary since we
% are using zeros instead of missing values

Z=[];

for i=1:9:length(n)-8
    z=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       n(i) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 n(i) n(i+1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 n(i) n(i+1) n(i+2) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 n(i) n(i+1) n(i+2) n(i+3) 0 0 0 0 0 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 n(i+1) n(i+2) n(i+3) n(i+4) 0 0 0 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 n(i+2) n(i+3) n(i+4) n(i+5) 0 0 0 0;...
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n(i+3) n(i+4) n(i+5) n(i+6)];
   Z=[Z;z];
end

% Define weighting matrix
W=[2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
   -1 2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

j = 1;

for i=1:19
    w1 = [zeros(1,j), -1, 2, -1, zeros(1,19-j)];
    W = [W; w1];
    j = j+1;
end

w2 = [zeros(1,20), -1, 2];
W = [W; w2] ;

% Define lagged differenced n
n_delta_1=zeros(length(n_delta),1);

for i=1:length(n_delta)-1
    n_delta_1(i+1) = n_delta(i);
end

% Dropping observations for the lowest t
for i=1:9:length(n_delta_1)-8 
    n_delta_1(i)=0;
    n_delta(i)=0;
    k_delta(i)=0;
    y_delta(i)=0;
    w_delta(i)=0;
end

% Define matrix of independent variables
X_ab = horzcat(n_delta_1, k_delta, w_delta, y_delta);

% Compute Arellano-Bond estimator 
gamma_ab = (X_ab'*Z*W*Z'*X_ab)^(-1)*X_ab'*Z*W*Z'*n_delta;

% Save matrix of residuals
for i=1:max(firm)
    e_ab(:,i)=n_delta((i-1)*9+1:i*9,1)-X_ab((i-1)*9+1:i*9,:)*gamma_ab;
end

% Computing standard errors of Arellano-Bond model
e = y_delta - X_ab*gamma_ab;
sigma2=(e'*e)/(max(firm)-4);
cov_beta = sigma2*inv(X_ab'*X_ab); 
se_ab=sqrt(diag(cov_beta));

end