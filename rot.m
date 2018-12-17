% Produce rotation matrix.
% Details: input is dimension n of vector space for rotation. 
% A random permutation of 1:n is generated and stored as indx. 
% Then for i=1,...,n-1, a rotation of pi/4 rad between indx(i) 
% and indx(i+1) is performed. 
% The resulting rotation matrix Q is output.
% Carl Ehrett 2017

function [Q] = rot(n)

rindx = randperm(n); 

Q = eye(n);
stot=sqrt(2)/2;
msg=0;

for ii = 1 : (n-1)
    rii = rindx(ii);
    rjj = rindx(ii+1);
    q = speye(n);
    q(rii,rii) = stot;
    q(rjj,rjj) = stot;
    q(rjj,rii) = stot;
    q(rii,rjj) = -stot;
    Q = q * Q;
    fprintf(repmat('\b',1,msg));
    msg=fprintf('Current column: %g/%g\n', ii,n);
end

end

% Saving the below as backup 
% function [Q] = rot(n)
% 
% rindx = randperm(n); 
% 
% Q = eye(n);
% stot=sqrt(2)/2;
% msg=0;
% 
% for ii = 1 : (n-1)
%     rii = rindx(ii);
%     for jj = (ii+1):n
%         rjj = rindx(jj);
%         q = speye(n);
%         q(rii,rii) = stot;
%         q(rjj,rjj) = stot;
%         q(rjj,rii) = stot;
%         q(rii,rjj) = -stot;
%         Q = q * Q;
%     end
%     fprintf(repmat('\b',1,msg));
%     msg=fprintf('Current column: %g/%g\n', ii,n);
% end
% 
% end