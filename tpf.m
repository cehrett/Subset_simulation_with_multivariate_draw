% Toy performance function
% Just distance from the origin scaled by a hyperellipsoid.
% The jth semi-principal axis is of length 1/(j^2*sqrt(c)), c given below.

function [gval] = tpf(gsettings,x) % settings here includes the rotation mat
% Unpack settings
Q = gsettings.rotation_matrix ;
idx = gsettings.hyperellipse_indices ; 

y = 2 * x - 1; % Convert from [0,1] to [-1,1]
unrotx = Q * y'; % Q "unrotates" so the hyperellipse is axis-aligned
% The axes of the hyperellipse are given by idx
pt = unrotx(idx);
n=length(idx); % Dimensionality of the ellipse
c=1;
sdist = 0; % We'll add to this in below loop to get point's scaled norm
for ii = 1:n
    %denoms[ii]=ii^2
    sdist = sdist + pt(ii)^2 * ii^4 * c;
end
gval=1/sqrt(sdist);
return;

