function hrot = bwrot(data,irot) 
% BWROT, Silverman's Rule of Thumb BandWidths, for 1-d kde
%     from Silverman's density estimation book
%     Variations on the Simple normal reference:  bwsnr.m
%     Assumes the kernel is standard Gaussian
% Inputs:
%     data - column vector of data
%     irot - Silverman's version:
%               1 - uses Silverman's eqn. 3.28 + the eqn 3.30 adjustment
%               2 (or unspecified) - uses eqn. 3.31
% Output:
%     hrot - simple "Rule Of Thumb bandwidth"
%
% Used by: gpkde
%
% Assumes path can find personal functions:
%    iqr
%    cquant

%    Copyright (c) J. S. Marron 1996

if nargin == 1 ;   %  Then use default: ROT2
  iirot = 2 ;
else ;             %  Use what was read in
  iirot = irot ;
end ;

n = length(data) ;
dsd = std(data) ;
diqr = iqr(data) ;


a = min([dsd; (diqr / 1.34)]) ;
          %  Using Silverman's notation

if iirot == 1 ;
  hrot = 1.06 * a * n^(-1/5) ;
else ;
  hrot = .9 * a * n^(-1/5) ;
end ;

