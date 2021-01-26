function [kde,xgrid,mker] = gpkde(data,vh,vxgrid,imptyp,eptflag) 
% GPKDE, General Purpose Kernel Density Estimate (1-d, Gaussian Kernel)
%     Does 1-d kernel density estimation, using binned (default) or 
%     direct (either matrix, or loops for bigger data sets), 
%     implementations, with the bandwidth either user specified 
%     (can be vector), or data driven (SJPI, Normal Reference, 
%     Silverman's ROT2, Oversmoothed).
%   Can use first 1, 2, 3, 4 or 5 arguments.
% Inputs:
%     data   - either n x 1 column vector of 1-d data
%                or vector of bincts, when imptyp = -1
%     vh     - vector of bandwidths, or specifies data driven:
%                  0 (or not specified)  -  Sheather Jones Plug In
%                  -1  -  Simple Normal Reference
%                  -2  -  Silverman's Rule Of Thumb 2 
%                            (20% Smaller than min of sd and IQR)
%                  -3  -  Oversmoothed
%                      Note: <0 only works for imptype = 0
%                  >0  -  Use input number (numbers if vector)
%                      Note: this MUST be >0 for imptyp >= 1
%                                 (the direct implementations)
%     vxgrid - vector of parameters for, or values of, grid to evaluate at:
%                  0 (or not specified)  -  use endpts of data and 401 bins
%                  [le; lr]  -  le is left end, re is right, 401 bins
%                         (get error message and no return if le > lr)
%                  [le; lr; nb] - le left, re right, and nb bins
%                  xgrid  -  Use input values 
%                            Note:  need to have more than 3 entries,
%                                 and only works when imptyp = 1 or 2
%    imptyp  - flag indicating implementation type:
%                 -1  -  binned version, and "data" is assumed to be
%                                   bincounts of prebinned data
%                  0 (or not specified)  -  linear binned version
%                                   and bin data here
%                  1  -  Direct matrix implementation
%                  2  -  Slow looped implementation (only useful when
%                            1 creates matrices that are too large)
%    eptflag - endpoint truncation flag (only has effect when imptyp = 0):
%                  0 (or not specified)  -  move data outside range to
%                                   nearest endpoint
%                  1  -  truncate data outside range
% Output:
%     (none)  -  Draws a graph of the result (in the current axes)
%     kde     -  col vector of heights of kernel kernel density estimate,
%                    unless vh is a vector (then have matrix, with
%                    corresponding cols as density estimates)
%     xgrid   -  col vector grid of points at which estimate(s) are 
%                    evaluted (useful for plotting), unless grid is input,
%                    can also get this from linspace(le,re,nb)'  
%     mker    -  matrix (vector) of kernel functions, evaluated at xgrid,
%                    which can be plotted to show "effective window 
%                    sizes" (currently scaled to have mass 0.05, may
%                    need some vertical rescaling)
%
% Used by:   gpfam.m 
%
% Assumes path can find personal functions:
%    vec2mat.m
%    gplbinr.m
%    bwsjpib.m
%    bwos.m
%    bwrot.m
%    bwsnr.m

%    Copyright (c) J. S. Marron 1996-1997



%  Set parameters and defaults according to number of input arguments
%
if nargin == 1 ;    %  only 1 argument input, use default bandwidth
  ivh = 0 ;         %  use default SJPI
else ;              %  bandwidth was specified, use that
  ivh = vh ;
end ;

if nargin <= 2 ;    %  at most 2 arguments input, use default xgrid
  ivxgrid = 0 ;
else ;              %  xgrid was specified, use that
  ivxgrid = vxgrid ;
end ;

if nargin <= 3 ;    %  at most 3 arguments input, use default implementation
  iimptyp = 0 ;
else ;              %  implementation was specified, use that
  iimptyp = imptyp ;
end ;

if nargin <= 4 ;    %  at most 4 arguments input, use default endpt trunc
  ieptflag = 0 ;    %  Default
else ;
  ieptflag = eptflag ;    %  Endpt trunc was specified, so use it
end ;



%  Calculate kde
%
if iimptyp > 0 ;    %  Then do direct implementation

  if min(ivh) > 0 ;    %  Then have valid bandwidths, so proceed

    n = length(data) ;

    if length(ivxgrid) > 3 ;  %  Then use input grid
      xgrid = ivxgrid ;
      nbin = length(xgrid) ;
    else ;                    %  Need to generate a grid
      nbin = 401 ;         %  Default
      lend = min(data) ;   %  Default
      rend = max(data) ;   %  Default
      if length(ivxgrid) >= 2 ;      %  use input endpoints
        lend = ivxgrid(1) ;
        rend = ivxgrid(2) ;
      end ;
      if length(ivxgrid) == 3 ;      %  use number of grid points
        nbin = ivxgrid(3) ;
      end ;

      if lend > rend ;    %  Then bad range has been input
        disp('!!!   Error in gpkde: invalid range chosen  !!!') ;
        xgrid = [] ;
      else ;
        xgrid = linspace(lend,rend,nbin)' ;
      end ;
    end ;


    %  Loop through bandwidths
    kde = [] ;
    for ih = 1:length(ivh) ;
      h = ivh(ih) ;

      if iimptyp ~= 2 ;  %  Then do direct matrix implementation
        kdeh = vec2mat((data ./ h),nbin) - vec2mat((xgrid' ./ h),n) ;
          %  efficient way to divide all dif's by h
          %  variable name "kde" is used to avoid creating too many biggies
        kdeh = exp(-(kdeh .^2) / 2) ;
          %  exponential part of Gaussian density
        kdeh = sum(kdeh)' ;
          %  sum part of kde, and make result a column vector
        kdeh = kdeh / (n * h * sqrt(2 * pi)) ;
          %  normalize, and mult by Gaussain density constant
        kde = [kde kdeh] ;
      else ;   %  Do slower looped implementation
        kdeh = [] ;
        for ixg = 1:nbin ;    %  Loop through grid points
          kdehx = (data - xgrid(ixg)) / h ;
          kdehx = sum(exp(-(kdehx .^2) / 2)) ;
          kdeh = [kdeh; kdehx] ;
        end ;
        kdeh = kdeh / (n * h * sqrt(2 * pi)) ;
        kde = [kde kdeh] ;
      end ;
    end ;

  else ;    %  Have invalid bandwidths

    disp('!!!   Error in gpkde: A bandwidth is invalid   !!!') ;
    disp('    (Note: cannot use data driven, with direct impl''s)') ;

  end ;

else ;     %  Then do binned implementation

  if iimptyp == -1 ;   %  Then data have already been binned

    if (length(ivxgrid) == 1) | (length(ivxgrid) > 3) ;
                         %  Then can't proceed because don't have bin ends
      disp('!!!   Error: gpkde needs to know the endpoints   !!!') ;
      disp('!!!            to use this implementation        !!!') ;
      bincts = [] ;
    else ;
      n = sum(data) ;
      bincts = data ;

      nbin = 401 ;
      lend = ivxgrid(1) ;
      rend = ivxgrid(2) ;
      if length(ivxgrid) == 3 ;          %  then use number of grid points
        nbin = ivxgrid(3) ;
      end ;

      if nbin ~= length(bincts) ;    %  Then something is wrong
        disp('!!!   Warning: gpkde was told the wrong number of bins   !!!') ;
        disp('!!!            will just use the number of counts.       !!!') ;
        nbin = rows(bincts) ;
      end ;
    end ;

  else ;               %  Then need to bin data

    if length(ivxgrid) > 3 ;  %  Then need to warn of change to default
      disp('!!!   Warning: gpkde was given an xgrid, and also   !!!') ;
      disp('!!!       asked to bin; will bin and ignore xgrid   !!!') ;
    end ;

    %  Specify grid parameters
    nbin = 401 ;         %  Default
    lend = min(data) ;   %  Default
    rend = max(data) ;   %  Default
    if (length(ivxgrid) == 2) | (length(ivxgrid) == 3) ;
                                     %  then use input endpoints
      lend = ivxgrid(1) ;
      rend = ivxgrid(2) ;
    end ;
    if length(ivxgrid) == 3 ;          %  then use number of grid points
      nbin = ivxgrid(3) ;
    end ;

    if lend > rend ;    %  Then bad range has been input
      disp('!!!   Error in gpkde: invalid range chosen  !!!') ;
      bincts = [] ;
    else ;
      bincts = gplbinr(data,[lend,rend,nbin],ieptflag) ;
    end ;

    %  Can do data-based bandwidth selection here, if specified
    if ivh == -1 ;        %  Then use Simple Normal Reference
      ivh = bwsnr(data) ;
    elseif ivh == -2 ;    %  Then use Silverman's Rule Of Thumb 2 
                          %  (~10% Smaller than min of sd and IQR)
      ivh = bwrot(data) ;
    elseif ivh == -3 ;    %  Then use Terrell's Oversmoother
      ivh = bwos(data) ;
    elseif min(ivh) <= 0 ;     %  Then be sure to use default SJPI 
                          %    (in case an unsupported value was input)
      ivh = 0 ;
    end ;

  end ;
  n = sum(bincts) ;
          %  put this here in case of truncations during binning

  %  Get bandwidth (if still not yet specified)
  if ivh == 0 ;    %  Then use SJPI bandwidth
      ivh = bwsjpib(bincts,[lend; rend],0,-1) ;
  end ;


  %  Loop through bandwidths
  kde = [] ;
  for ih = 1:length(ivh) ;
    h = ivh(ih) ;

    %  Create vector of kernel values, at equally spaced grid
    delta = (rend - lend) / (nbin - 1) ;    %  binwidth
    k = nbin - 1 ;    %  index of last nonzero entry of kernel vector
    arg = linspace(0,k * delta / h,k + 1)' ;
    kvec = exp(-(arg.^2) / 2) / sqrt(2 * pi) ;
    kvec = [flipud(kvec(2:k+1)); kvec] ;

    %  Do actual kernel density estimation
    kdeh = conv(bincts,kvec) ;
    kdeh = kdeh(k+1:k+nbin) / (n * h) ;

    if h < 3 * delta ;    %  Then need to normalize
                             %  to make numerical integral roughly 1
      kdeh = kdeh / (sum(kdeh) * delta) ;
    end ;

    kde = [kde kdeh] ;
  end ;

  xgrid = linspace(lend,rend,nbin)' ;

end ;


%  Create matrix of kernels, if this is needed
%
if nargout == 3 ;
  cent = mean([lend; rend]) ;
          %  centerpoint of evaluation grid
  if length(ivh) > 1 ;
    mih = vec2mat(1 ./ ivh',nbin) ;
    mker = vec2mat((xgrid - cent),length(ivh)) .* mih;
          %  argument of gaussian kernel
  else ;
    mih = 1 / ivh ;
    mker = (xgrid - cent) .* mih;
          %  argument of gaussian kernel
  end ;
  mker = exp(-mker.^2 / 2) .* mih / sqrt(2 * pi) ;
          %  Gaussian kernels with mass 1
  mker = 0.05 * mker ;
          %  Make masses = 0.05
end ;



%  Make plots if no numerical output requested
%
if nargout == 0 ;  %  Then no numerical output, but make a plot
  plot(xgrid,kde) ;
end ;

