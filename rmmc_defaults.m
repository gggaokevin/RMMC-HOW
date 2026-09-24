function [lambda, rho1, rho2] = rmmc_defaults(m, n)
% RMMC_DEFAULTS  reported RMMC settings, switched by the DATA FAMILY.
%
% V6 (observed-support IQR, eta = (2.5,3.5)): ONE penalty pair serves BOTH data
% families, (rho1,rho2) = (9.05,1.45).  With v = 1.74331 (5v = 8.7165) it satisfies
% the sufficient conditions eq. (9) of the paper at either lambda:
%
%   lambda = 5e-3 : binding condition rho1 > 5v = 8.7165          (margin +3.83%)
%   lambda = 0.1  : binding condition rho1 > (sqrt(41a^2+20)-a)/2
%                   = 9.0176 with a = 1+rho2+2*lambda*tau_L       (margin +0.36%)
%
% The pair is the smallest one on the 0.05 grid that works for both families: on the
% ORL side the smallest feasible rho1 is 9.05 and the admissible window is
% rho2 in [1.44,1.46]; the unconstrained minimum 8.9945 sits at rho2 ~ 1.441.
%
%   synthetic S1-S3 : 100 x 200 slices  ->  lambda = 5e-3, (rho1,rho2) = (9.05,1.45)
%   ORL       S4    : 112 x  92 slices  ->  lambda = 0.1,  (rho1,rho2) = (9.05,1.45)
%
% The switch keys on the slice geometry (m,n) only, so the K-frame subsets used by
% the selection studies resolve to the same family.  ``RMMC.m`` calls this when the
% caller does not pass lambda / rho1 / rho2 explicitly; explicit arguments always
% win.
%
% Unknown geometry: warn and fall back to lambda = 0.1 (the more conservative of the
% two), so that a forgotten argument is loud rather than silent.

rho1 = 9.05; rho2 = 1.45;
if m == 100 && n == 200
    lambda = 5e-3;
elseif m == 112 && n == 92
    lambda = 0.1;
else
    warning('RMMC:unknownGeometry', ...
        ['slice size %dx%d matches no reported data family; falling back to ' ...
         'lambda=0.1 with (rho1,rho2)=(9.05,1.45) -- pass the parameters explicitly'], m, n);
    lambda = 0.1;
end
end
