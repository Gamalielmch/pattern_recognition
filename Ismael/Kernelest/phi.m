function normprob = phi(data) 
% PHI, Standard Gaussian (normal) c.d.f.
% Input:
%     data  - column vector, or matrix of data
% Output:
%     normprob  - for each entry in data, returns P(Z <= data)
%                   where Z has the Standard (mean 0, variance 1)
%                   Normal Distribution
%
%    For c.d.f. of X ~ N(my,sig^2), use phi((data - mu) / sig) 
%
%    Inverse is phiinv.m
%
%    Copyright (c) J. S. Marron 1997


%  Do main calculation
%
normprob = 0.5 * (1 + erf(data / sqrt(2))) ;


%  Truncate anything below 0, or above 1 (caused by round off errors)
%
flag = normprob < 0 ;
nflag = sum(sum(flag)) ;
if nflag > .5 ;    %  then need to fix negative values
  normprob(flag) = zeros(nflag,1) ;
end ;
%
flag = normprob > 1 ;
nflag = sum(sum(flag)) ;
if nflag > .5 ;    %  then need to fix values > 1
  normprob(flag) = ones(nflag,1) ;
end ;

