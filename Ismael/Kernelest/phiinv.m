function normq = phiinv(prob) 
% PHIINV, Inverse of Standard Gaussian (normal) c.d.f.
% Input:
%     prob  - matrix of probabilities 
%                  (assumed between 0 and 1)
% Output:
%     normq  - corresponding Standard Normal (mean 0, variance 1) Quantiles
%
%    For quantiles of X ~ N(my,sig^2), use sig * phiinv(prob) + mu
%
%    Inverse is phi.m
%
%    Copyright (c) J. S. Marron 1997


%  Do main calculation
%
normq = sqrt(2) * erfinv((2 * prob) - 1) ;

