function [theta_ip1, gval] = MAS(theta_i,pi_prop,Sigma,...
    g,b,gsettings,mma,gval_old)

if (nargin ~= 8)
   error('Incorrect number of parameters in function "MAS"');
end

%% Initialize variables
n      = length(theta_i);   % Number of uncertain parameters
xi     = zeros(1,n);        % \xi  
hat_xi = zeros(1,n);        % \tilde{\xi}  

%% 1. Sampling from proposal PDF with media the current state
if mma==1 % mma tells whether use Au/Beck MMA algorthm or new multivar draw
    accept = 0; % Change accept below if any value changes in theta
    for k = 1:n % We loop through each element of theta
       hat_xi(k) = pi_prop(theta_i(k),Sigma); % Propose new theta_i(k)
       % Get likelihood ratio: must use MH since we using truncated normals
       lhr_num = (normcdf(1-theta_i(k))-normcdf(-theta_i(k)) ) * (hat_xi(k)<1) * (hat_xi(k)>0);
       lhr_den = normcdf(1-hat_xi(k)) - normcdf(-hat_xi(k));
       r = lhr_num/lhr_den;
       
       if rand <= min(1,r) % Accept new draw with probability r
          xi(k) = hat_xi(k);    % Accept the candidate
          accept = 1;
       else
          xi(k) = theta_i(k);   % Reject the candidate and use the same state
       end
    end
else % If using new multivar draw instead of Au/Beck MMA algorithm:
    hat_xi = pi_prop(theta_i,Sigma); % Get proposal for new theta_i draw
    % Get likelihood ratio on log scale:
    lhr_num = sum( log(hat_xi .* (1-hat_xi)));
    lhr_den = sum( log(theta_i .* (1-theta_i)));
    r = (lhr_num-lhr_den); % Here's the likelihood ratio
    if log(rand) <= r % We accept the new draw with probability r
        xi = hat_xi;
        accept=1;
    else
        xi = theta_i;
        accept=0;
    end
end

%% 2. Check the location of xi by system analysis
% The support for the distribution we're drawing from is limited to those
% points with g-value above the current threshold b. So if our newly
% accepted draw has g-value below b, then it falls outside of the support
% of the distribution we're trying to draw from, and we should reject it
% after all. So here we test its g-value, and reject the new draw if it
% does not rise above the current threshold b.
if accept == 1 % We only need to check the g-value if we accepted new draw
    gval = g(gsettings,xi);
    if gval >= b
        theta_ip1 = xi;          % Accept the candidate in failure region
        new = 1;
    else
        theta_ip1=theta_i;
        gval = gval_old;
        new=0;
    end
else
    gval = gval_old;
    theta_ip1 = theta_i;
    new = 0;
end

return;
%%END