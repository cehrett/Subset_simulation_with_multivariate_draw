% Master file for running Multivariate Draw Subset Simulation
clear; close all; clc; rng(1);

%% Initial settings
%%% Configuration parameters
mma      = 0 ;           % flag; set to 0 to use new multivariate draw, 
                         %%%% to 1 for Au/Beck MMA algorithm
N        = 1000;         % Total number of samples for each level
p0       = 0.1;          % Probability of each adaptively chosen subset

%%% Set performance function, and settings for it:
g = @tpf;
n = 1000; % Dimensionality of the hyperellipse
load tpfparams
% tpfparams.mat contains two things:
% 1. rotation_matrix - included so that the hyperellipse of which we're
% finding the volume gets rotated. The semiprincipal axes of the resulting
% hyperellipse are not axis-aligned in the supporting R^n space. This
% matrix was produced using the rot.m file, with rng(1). Since the volume
% we're trying to find is actually the intersection of a hyperellipse with
% a hypercube, the volume of the resulting intersection changes depending
% on how the hyperellipse is rotated in the hypercube. Which is to say that
% a different rotation than the one provided here will yield a different
% volume.
% 2. hyperellipse_indices - this tells which of the n dimensions are the
% semiprincipal axes of the (unrotated, axis-aligned) hyperellipse.

%%% Set threshold defining region of interest: this defines the region of
% interest to be those points x such that g(x)>B. When g is the function
% tpf.m, B=1 corresponds to finding points that are within the
% hyperellipsoid with four semiprincipal axes, the jth semi-principal axis 
% of which is of length 1/(j^2*sqrt(2.314).
B = 1;

%%% Proposal and target distributions
pi_targ = @() rand(n,1); % Prior distribution
if mma == 1 % Recall mma tells whether use Au/Beck MMA or 
            % new multivar draw w/ logit transform
    pi_prop = @(x,Sigma) (Sigma * normrnd(0,1) + x) ;
else
    logit_trans = @(x) log(x ./ (1-x));
    logit_rev_trans = @(x) exp(x) ./ (1+exp(x));
    pi_prop = @(x,Sigma) logit_rev_trans(mvnrnd(logit_trans(x),Sigma));
end

%% Perform Subset simulation
[Pf_SS,Pf,gsort,b,F_total,F_seeds,...
    theta_rec,theta_rec_u,uniques,Nf,geval] = ...
    SS(n,N,p0,B,pi_targ,pi_prop,g,gsettings,mma);
fprintf('\n***SubSim Pf: %g ***\n', Pf_SS);

