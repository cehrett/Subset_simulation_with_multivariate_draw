% Toy performance function
% Just distance from the origin scaled by a hyperellipsoid.
% The jth semi-principal axis is of length 1/(j^2*sqrt(2.314)).

function [gval] = tpf(gsettings,x) % settings here includes the rotation mat
% Unpack settings
Q = gsettings.rotation_matrix ;
idx = gsettings.hyperellipse_indices ; 

y = 2 * x - 1; % Convert from [0,1] to [-1,1]
unrotx = Q * y'; % Q "unrotates" so the hyperellipse is axis-aligned
% The axes of the hyperellipse are given by idx
pt = [unrotx(idx(1)); unrotx(idx(2)); unrotx(idx(3)); unrotx(idx(4))];
n=4; % Dimensionality of the ellipse
c=2.314; % Normalizing constant
sdist = 0; % We'll add to this in below loop to get point's scaled norm
for ii = 1:n
    %denoms[ii]=ii^2
    sdist = sdist + pt(ii)^2 * ii^4 * c;
end
gval=1/sqrt(sdist);
return;

