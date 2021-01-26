function hsnr = bwsnr(data) 
% BWSNR, Simple Normal Reference BandWidth, for 1-d kernel density estimation
%     This is the MISE asy optimal bandwidth, if the data are 
%     normally distributed, the sample variance as variance
%     Assumes the kernel is standard Gaussian
% Inputs:
%     data - column vector of data
% Output:
%     hsnr - simple "Normal reference bandwidth"
%
% Used by: gpkde
%

%    Copyright (c) J. S. Marron 1996

n = length(data) ;
dsd = std(data) ;

  rk = 1 / (2 * sqrt(pi)) ;
          %  Integral of the square of the Gaussian kernel
  rfpp = 3 / (sqrt(pi) * 8) ;
          %  Integral of the square of the second derivative
          %      of the Gaussian density
c = (rk / rfpp)^(1/5) ;

hsnr = c * dsd * n^(-1/5) ;

