% Master file for running Multivariate Draw Subset Simulation
clear; close all; clc;

%% Initial settings
%%% Configuration parameters
N        = 1000;         % Total number of samples for each level
p0       = 1/10;          % Probability of each adaptively chosen subset

%%% Set performance function, and settings for it:
g = @tpf;
% B is the threshold defining region of interest: it defines the region of
% interest to be those points x such that g(x)>B. When g is the function
% tpf.m, B=1 corresponds to finding points that are within the
% hyperellipsoid with four semiprincipal axes, the jth semi-principal axis 
% of which is of length 1/(j^2).
B = 1;
n = 100; % Dimensionality of the hyperellipse
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

%% Perform Subset simulation
% Pf_SS is the estimated probability of the failure region. When g is the
% function tpf.m and B is 1 with the rotation defined in tpfparams, the
% true value here is approximately  1.5e-3. The estimate Pf_SS will vary
% around this value.

m = 10 ; 
% We'll run each of the two versions of SS m times, recording each time the
% estimating probability Pf_SS of the failure region and the time taken to
% run the algorithm. After the two algorithms have each run m times, we can
% compare their estimates and the time taken to achieve those estimates.

%%% First use the ``stock'' version of SS
mma      = 1 ;           % flag; set to 0 to use new multivariate draw, 
                         %%%% to 1 for Au/Beck MMA algorithm
Pf_SS_stock = zeros(m,1);
times_stock = zeros(m,1);
for ii=1:m
    tic;
    fprintf('\n\n==================');
    fprintf('\nSTOCK SS LOOP %g/%g',ii,m);
    fprintf('\n==================\n\n');
    [Pf_SS,Pf,gsort,b,F_total,F_seeds,...
        theta_rec,theta_rec_u,uniques,Nf,geval] = ...
        SS(n,N,p0,B,g,gsettings,mma);
    fprintf('\n***SubSim Pf: %g ***\n', Pf_SS);
    Pf_SS_stock(ii) = Pf_SS;
    times_stock(ii) = toc;
end

%%% Now use the new proposed multivariate draw version of SS
mma      = 0 ;           % flag; set to 0 to use new multivariate draw, 
                         %%%% to 1 for Au/Beck MMA algorithm
Pf_SS_mvd = zeros(m,1);
times_mvd = zeros(m,1);
for ii=1:m
    tic;
    fprintf('\n\n==============================');
    fprintf('\nMULTIVARIATE DRAW SS LOOP %g/%g\n\n',ii,m);
    fprintf('\n==============================\n\n');
    [Pf_SS,Pf,gsort,b,F_total,F_seeds,...
        theta_rec,theta_rec_u,uniques,Nf,geval] = ...
        SS(n,N,p0,B,g,gsettings,mma);
    fprintf('\n***SubSim Pf: %g ***\n', Pf_SS);
    Pf_SS_mvd(ii) = Pf_SS;
    times_mvd(ii) = toc;
end

%%% Now, compare the results
fprintf('\n=====================');
fprintf('\nCOMPARISON OF RESULTS');
fprintf('\n=====================\n');
fprintf...
    ('\n\nMean estimate found via stock SS: %g\n',mean(Pf_SS_stock));
fprintf('Std. dev. of estimate found via stock SS: %g\n',std(Pf_SS_stock));
fprintf('Mean time in seconds to complete stock draw SS: %g\n',...
    mean(times_stock));
fprintf...
    ('\nMean estimate found via multivariate draw SS: %g\n',...
    mean(Pf_SS_mvd)); 
fprintf('Std. dev. of estimate found via multivariate draw SS: %g\n',...
    std(Pf_SS_mvd));
fprintf('Mean time in seconds to complete multivariate draw SS: %g\n',...
    mean(times_mvd));


