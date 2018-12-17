%% Subset Simulation

clear; close all; clc; rng(1);
tic;

%% Initial data
% Configuration parameters
mma      = 0 ;           % flag; set to 0 to use new multivariate draw, 
                         %%%% to 1 for Au/Beck MMA algorithm
tvn      = 1 ;           % flag; set to 0 for earthquake, 1 for toy example
N        = 1000;         % Total number of samples for each level
p0       = 0.1;          % Probability of each subset, chosen adaptively

% Earthquake scenario parameters
pga      = 1 ;           % peak ground acceleration
dt       = 2;            % for time vector
t_vec    = dt:dt:300;    % time vector

%%%% Performance Function: Earthquake or toy
if tvn; g = @tpf; else g = @l29ft_v8_50pc; end

%%%% Settings
if tvn
    load tvnparams % Contains rotation matrices Q,Qi=Q^{-1}
    pga=Q; % Just to simplify the code, we'll sneak in Q as peak ground acc. param
else
    load loc29params_v7 % Contains parameter means (unnormalized)
    n    = numel(C);
    C    = reshape(C,[1,n])';
    theta_bounds = [ 0.6 * C  1.4 * C ]; 
end

% Max demand, marginal PDFs of the random variables and limit-state function
if tvn
    B  = 1;       % critical threshold defining failure region of interest
else
    B  = 0.016;   % critical threshold defining failure region of interest
end

%%%% Logit transformations
logit_trans = @(x) log(x ./ (1-x));
logit_rev_trans = @(x) exp(x) ./ (1+exp(x));

%%%% Proposal and target distributions

pi_targ = @() rand();
if mma == 1 % Recall mma is flag for whether to use Au/Beck MMA or new multivariate draw with logit transform
    pi_prop = @(x,Sigma) (Sigma * normrnd(0,1) + x) ;
else
    pi_prop = @(x,Sigma)...
        logit_rev_trans(mvnrnd(logit_trans(x),Sigma));
end

%% Subset simulation
[Pf_SS,Pf,gsort,b,F_total,F_seeds,...
    theta_rec,theta_rec_u,uniques,Nf,geval] = ...
    SS(n,N,p0,B,pi_targ,pi_prop,g,pga,t_vec,mma,tvn);
fprintf('\n***SubSim Pf: %g ***\n', Pf_SS);
t=toc;
Nf

%%END