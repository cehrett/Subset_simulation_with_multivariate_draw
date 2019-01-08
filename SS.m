function [Pf_SS,Pf,gsort,b,F_total,F_seeds,...
    theta_rec,theta_rec_u,uniques,Nf,geval] = ...
    SS(n,N,p0,B,pi_targ,pi_prop,g,gsettings,mma)
% Subset Simulation 

% Check inputs

if (nargin ~= 9)
   error('Incorrect number of parameters in function "SS"');
end
tic;

% Define logit
logit = @(x) log(x ./ (1-x));

%% Initialization of variables and storage
j      = 1; % initial conditional level
Ns     = 1/p0; % number of samples simulated form each Markov chain
if mod(N,Ns) ~= 0
    N = N + Ns-mod(N,Ns) ; % Necessary to make Nc an integer
end
Nc     = N*p0; % number of markov chains
est_it = 20; % estimated # of iteratns (serves as ad hoc max # of iteratns)
theta_j  = zeros(n,N);
geval    = zeros(1,N);
gsort    = zeros(1,N);
Nf       = zeros(est_it,1);
Pf       = zeros(est_it,1);
b        = zeros(est_it,1);
F_total  = cell(est_it,1);
F_seeds  = cell(est_it,1);
theta_l  = zeros(n,Ns);
theta_jk = cell(Nc,1);
theta_rec = cell(1,1);
theta_rec_u = cell(1,1);
uniques  = zeros(est_it,1);
Sigma  = 1e-1; % This is changed below before being used if mma==0
sig_c1 = 1e-1;  % These constants scale the variance and the identity 
sig_c2 = 1e-6; % addition to it to improve convergence; 1e0,1e-8; 1e1,1e-8
               % work well for hyperellipsoid, earthquake respectively

%%%% SS procedure
%%%% initial MCS stage
fprintf('Evaluating performance function:\t');
msg=0;
for i = 1:N
   theta_j(:,i) = pi_targ();
   geval(j,i)   = g(gsettings,theta_j(:,i)');
   if geval(j,i) >= B
      Nf(j) = Nf(j)+1;
   end
   
   % Update user via console
   fprintf(repmat('\b',1,msg));
   msg = fprintf('\n**Evaluating sample no. %g/%g',i,N);
   
end
% Update console
fprintf('\nInitial MCS complete! Nf(1)=%g\n',Nf(1));

%%%% SS stage
last = 0; % This indicator will turn to 1 when we're on final loop
while last ~= 1
   if Nf(j)/N >= p0
      last = 1    % need one more run to compute the curve
   end   
   j = j+1;   % next level
   
   % Sort values and define new intermediate failure region threshold
   [gsort(j-1,:),idx] = sort(geval(j-1,:)); % ascending
   F_total{j-1} = theta_j(:,idx); % store ordered level samples
   F_seeds{j-1} = theta_j(:,idx(N-Nc+1:N));% store ordered level seeds
   % Recall mma tells whether use Au/Beck MMA algorthm or new multivar draw
   if mma == 0 
       % Get new cov matrix for proposal dist:
       Sigma = sig_c1*(cov(logit(F_seeds{j-1})') + sig_c2 * eye(n)); 
   end
   b(j-1) = (gsort(j-1,N-Nc)+gsort(j-1,N-Nc+1))/2; % Get the p0 quantile
   fprintf('\n-Current threshold level %g is: %g \n', j-1, b(j-1));
   
   % Sampling process using MCMC
   gvals = zeros(N,1); % Will hold g values for this round of theta_j draws
   for k = 1:Nc % For each MCMC chain used this round:
      msg = fprintf('***Sampling from seed no. %g/%g', k, Nc);
      theta_l(:,1) = F_seeds{j-1}(:,k); % Set the seed for MCMC
      gvals((k-1)*Ns+1) = g(gsettings,theta_l(:,1)'); % Get seed's g-value
      for l = 1:Ns-1 % For each step in the MCMC chain:
         % Make a draw from the MCMC chain (using MAS function):
         [theta_l(:,l+1),gval] = ...
             MAS(theta_l(:,l),pi_prop,Sigma,g,...             
             b(j-1),gsettings,mma,gvals((k-1)*Ns+l));
         gvals((k-1)*Ns+1+l) = gval; % Record g-value of new draw
      end
      theta_jk{k} = theta_l; % This chain is recorded in theta_jk as cell k
      fprintf(repmat('\b',1,msg)); % Remove old console output
   end
   
   % Recordkeeping
   theta_j = cell2mat(theta_jk'); % Transform all chains drawn to matrix
   theta_rec{j-1}=theta_j; % Record matrix of this round's draws in a cell
   theta_rec_u{j-1}=unique(theta_j','rows'); % Record unique samples only
   % This records how many unique samples we got in each step:
   uniques(j-1) = size(theta_rec_u{j-1},1); 
   fprintf('Unique samples: %g\n',uniques(j-1));
   
   % Count new samples in failure region of interest
   fprintf('\nEvaluating performance function:\t');
   msg=0;
   for i = 1:N
      %geval(j,i) = g(pga,t_vec,theta_j(:,i)');
      geval(j,i) = gvals(i);
      if geval(j,i) >= B
         Nf(j) = Nf(j)+1;
      end
   end
   fprintf('\nOK! Nf(%g)=%g\n',j,Nf(j));
   
   if j > est_it % Change as desired to stop loop early; o/w set to est_it
       last=1;
   end
end

% delete unnecesary data
if j < est_it
  F_total(j:end) = [];
  F_seeds(j:end) = [];
  Nf(j:end)      = [];
  Pf(j:end)      = [];
  b(j:end)       = [];
end

% %% pf calculation
% Pf(1)        = p0;    
% Pf_line(1,:) = linspace(p0,1,1000);
% b_line(1,:)  = gsort(1,round((N-N*Pf_line(1,:))+1));
% for i = 2:j-1
%    Pf(i)        = Pf(i-1)*p0;
%    Pf_line(i,:) = Pf_line(i-1,:)*p0;
%    b_line(i,:)  = gsort(i,round((N-N*Pf_line(1,:))+1));
% end
% Pf_line = sort(Pf_line(:),'descend');  
% pf_ss   = Pf_line;
% b_line  = sort(b_line(:));
% b_ss    = b_line;
% Pf_SS   = p0^(j-2)*(Nf(j-1)/N);

Pf_SS   = 0; % This will be our prob estimate of F
% % OLD VERSION:
% for ii = 1:(j-1)
%    Pf_SS = Pf_SS + p0^(ii-1) * (Nf(ii)/N);
% end
% Pf_SS   = Pf_SS / (j-1);
% New version:
Pf_SS = p0^(j-2) * Nf(j-1)/N ;

% gvals = sort(gsort(:));
% b_line(1,:) = linspace(min(gvals),max(gvals),1000);
% Pf_line = zeros(1,size(b_line,2));
% for ii = 1:size(b_line,2)
%     pfv=0;
%     n_levels = min(sum(b<=b_line(ii))+1,j-1);
%     for jj = 1:n_levels
%         pfv = pfv + sum(gsort(jj,:)>= b_line(ii))/N * p0^(jj-1);
%     end
%     Pf_line(ii) = pfv/n_levels;
% end

%b_ss=b_line;
%pf_ss=Pf_line;

toc;

return;

end
%%END