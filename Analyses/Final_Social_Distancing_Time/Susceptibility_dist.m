% @ copyright
% Authors:
%   Ricardo Aguas
%   Rodrigo M Corder
%   Jessica G King
%   Guilherme Goncalves
%   Marcelo U Ferreira
%   M Gabriela M Gomes
%
% This work is protected under the @Attribution-NonCommercial 4.0 International intellectual property license.
% You are free to:
%   Share - copy and redistribute the material in any medium or format
%   Adapt - remix, transform, and build upon the material Under the following terms:
%   Attribution - You must give appropriate credit to the authors, and indicate if any changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
%   NonCommercial - You may not use the material for commercial purposes.
%   ShareAlike - If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

global xvec qvec var N

% Distribution of susceptibility given a defined variance:
n = 30;
nstep = 0.01;
xvec = (nstep:nstep:n);
N = n/nstep;
theta = max(var,1e-18);
k = 1/theta;
qvec = zeros(N,1);

qvec(1) = gamcdf(nstep,k,theta);
for j = 2:N
    qvec(j) = gamcdf(j*nstep,k,theta) - gamcdf((j-1)*nstep,k,theta);
end
qvec(N) = 1 - sum(qvec(1:N-1));

F = @(xN)... 
    [((xvec(1:N-1)-xvec(1:N-1)*qvec(1:N-1)-xN*qvec(N)).^2)*qvec(1:N-1)+(xN-xvec(1:N-1)*qvec(1:N-1)-xN*qvec(N))^2*qvec(N)-var;... 
    xvec(1:N-1)*qvec(1:N-1)+xN*qvec(N)-1];
warning('off')
options = optimset('Display','off');
[xN,eval] = fsolve(F,xvec(N),options);
xvec(N) = xN;
 