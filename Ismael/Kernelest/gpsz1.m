function [mapout,xgrid] = gpsz1(data,vxgp,vhgp,eptflag,alpha,simflag,llflag) 
% GPSZ1, General Purpose Significant derivative Zero crossings
%     Creates gray level map (function of location and bandwidth),
%     which is:
%         "dark"   at points where deriv sig > 0
%         "light"  at points where deriv sig < 0
%         "gray"   at points where deriv roughly 0
%         "darker gray" where "effective sample size" < 5
%   Can use first 1, 2, 3, 4, 5 or 6 arguments.
% Inputs:
%     data   - either n x 1 column vector of density estimation data
%                  or n x 2 matrix of regression data:
%                            X's in first column,  Y's in second
%     vxgp   - vector of x grid parameters:
%                  0 (or not specified)  -  use endpts of data and 401 bins
%                  [le; lr; nb] - le left, re right, and nb bins
%     vhgp   - vector of h (bandwidth) grid parameters:
%                  0 (or not specified)  -  use (2*binwidth) to (range),
%                                                and 21 h's
%                  [hmin; hmax; nh]  -  use hmin to hmax and nh h's.
%    eptflag - endpoint truncation flag (only has effect when imptyp = 0):
%                  0 (or not specified)  -  move data outside range to
%                                   nearest endpoint
%                  1  -  truncate data outside range
%    alpha   -  Usual "level of significance".  I.e. C.I.s have coverage
%                  probability 1 - alpha.  (0.05 when not specified)
%    simflag -  Confidence Interval type (simultaneous vs. ptwise)
%                  0  -  Use Pointwise C.I.'s
%                  1 (or not specified)  -  Use Simultaneous C.I.'s
%    llflag  -  Type of Regression Estimator (LL  vs. NW)
%                  0  -  Use Nadaraya Watson est. (in regression problems)
%                  1 (or not specified)  -  Use Local Linear est.
%                            (has no effect for density est)
%
% Output:
%     (none)  -  Draws a gray level map (in the current axes)
%     mapout  -  output of gray level map
%     xgrid   -  col vector grid of points at which estimate(s) are 
%                    evaluted (useful for plotting), unless grid is input,
%                    can also get this from linspace(le,re,nb)'  
%
% Assumes path can find personal functions:
%    vec2mat.m
%    gplbinr.m
%    phiinv.m
%

%    Copyright (c) J. S. Marron 1997


%  Set parameters and defaults according to number of input arguments
%
if nargin == 1 ;    %  only 1 argument input, use default vxgp
  ivxgp = 0 ;
else ;              %  xgrid was specified, use that
  ivxgp = vxgp ;
end ;

if nargin <= 2 ;    %  at most 2 arguments input, use default vhgp
  ivhgp = 0 ;
else ;              %  hgrid parameters were specified, use them
  ivhgp = vhgp ;
end ;

if nargin <= 3 ;    %  at most 3 arguments input, use default endpt trunc
  ieptflag = 0 ;    %  Default
else ;
  ieptflag = eptflag ;    %  Endpt trunc was specified, so use it
end ;

if nargin <= 4 ;    %  at most 4 arguments input, use default alpha
  ialpha = 0.05 ;    %  Default
else ;
  ialpha = alpha ;    %  Endpt trunc was specified, so use it
end ;

if nargin <= 5 ;    %  at most 5 arguments input, use default simul CI's
  isimflag = 1 ;    %  Default
else ;
  isimflag = simflag ;    %  Endpt trunc was specified, so use it
end ;

if nargin <= 6 ;    %  at most 6 arguments input, use default LL in reg
  illflag = 1 ;    %  Default
else ;
  illflag = llflag ;    %  Endpt trunc was specified, so use it
end ;


%  detect whether density or regression data
%
if size(data,2) == 1 ;    %  Then is density estimation
  xdat = data(:,1) ;
  idatyp = 1 ;
else ;                    %  Then assume regression ;
  xdat = data(:,1) ;
  ydat = data(:,2) ;
  idatyp = 2 ;
end ;

%  Set x grid stuff
%
n = length(xdat) ;
if ivxgp == 0 ;   %  then use standard default x grid
  ivxgp = [min(xdat),max(xdat),401] ;
end ;
left = ivxgp(1) ;
right = ivxgp(2) ;
ngrid = ivxgp(3) ;
xgrid = linspace(left,right,ngrid)' ;
          %  row vector to evaluate smooths at

%  Set h grid stuff
%
range = right - left ;
binw = range / (ngrid - 1) ;
if ivhgp == 0 ;   %  then use standard default h grid
  ivhgp = [2 * binw,range,21] ;
end ;
hmin = ivhgp(1) ;
hmax = ivhgp(2) ;
nh = ivhgp(3) ;
vh = logspace(log10(hmin),log10(hmax),nh) ;



%  Bin the data
%
if idatyp == 1 ;        %  Treating as density estimation
  bincts = gplbinr(xdat,ivxgp,ieptflag) ;
elseif idatyp == 2 ;    %  Treating as regression
  bincts = gplbinr([xdat ydat],ivxgp,ieptflag) ;
  bincts2 = gplbinr([xdat, ydat.^2],ivxgp) ;
end ;
n = sum(bincts(:,1)) ;
          %  put this here in case of truncations during binning



%  Construct Surfaces
%
mdsurf = [] ;
          %  Derivative surface
mesurf = [] ;
          %  Effective sample size surface
mvsurf = [] ;
          %  Estimated Variance of Derivative
vgq = [] ;
          %  Vector of Gaussian Quantiles (for simultaneous CI's)

%  Create grid for kernel values
delta = (right - left) / (ngrid - 1) ;    %  binwidth
k = ngrid - 1 ;    %  index of last nonzero entry of kernel vector

%  Loop through bandwidths
for ih = 1:nh ;
  h = vh(ih) ;

  %  Set common values
  arg = linspace(0,k * delta / h,k + 1)' ;
  kvec = exp(-(arg.^2) / 2) ;
  kvec = [flipud(kvec(2:k+1)); kvec] ;
        %  construct symmetric kernel


  %  Vector of Effective sample sizes 
  %          (notation "s0" is consistent with below)
  ve = conv(bincts(:,1),kvec) ;
  ve = ve(k+1:k+ngrid) ;
          %  denominator of NW est.
          %    (same as sum for kde)


  if idatyp == 1 ;       %  Doing density estimation
          %  main lines from gpkde.m, via nmsur5.m

    kvecd = -arg .* exp(-(arg.^2) / 2) / sqrt(2 * pi) ;
    kvecd = [-flipud(kvecd(2:k+1)); kvecd] ;
        %  construct symmetric kernel

    vd = conv(bincts,kvecd) ;
    vv = conv(bincts,kvecd.^2) ;
    vd = vd(k+1:k+ngrid) / (n * h^2) ;
    vv = vv(k+1:k+ngrid) / (n * h^4) ;
    vv = vv - vd.^2 ;
    vv = vv / n ;
          %  did this outsidea loop in nmsur5.m

  else ;    %    Doing regression
          %  using modification of lines from gpnpr.m

    if illflag == 0 ;    %  then do Nadaraya-Watson type regression
          %  main lines from gpnpr.m, via szeg2.m

      %  start with derivative kernel
      kvecd = -arg .* exp(-(arg.^2) / 2) ;
      kvecd = [-flipud(kvecd(2:k+1)); kvecd] ;
        %  construct symmetric kernel

      s1 = conv(bincts(:,1),kvecd) ;
      s1 = s1(k+1:k+ngrid) ;
      t0 = conv(bincts(:,2),kvec) ;
      t0 = t0(k+1:k+ngrid) ;
      t1 = conv(bincts(:,2),kvecd) ;
      t1 = t1(k+1:k+ngrid) ;
      sig2 = conv(bincts2(:,2),kvec) ;
      sig2 = sig2(k+1:k+ngrid) ;

      flag = (ve < 1) ;
          %  locations with effective sample size < 1
      ve(flag) = ones(sum(flag),1) ;
          %  replace small sample sizes by 1 to avoid 0 divide

      mhat = t0 ./ ve ;
      vd = (ve .* t1 - s1 .* t0) ./ ve.^2 ;
          %  derivative estimate

      sig2 = sig2 ./ ve - mhat.^ 2 ;
      v0 = conv(sig2,kvec.^2) ;
      v0 = v0(k+1:k+ngrid) ;
      v1 = conv(sig2,kvec .* kvecd) ;
      v1 = v1(k+1:k+ngrid) ;
      v2 = conv(sig2,kvecd.^2) ;
      v2 = v2(k+1:k+ngrid) ;
      vv = (ve.^2 .* v2 - 2 * ve .* s1 .* v1 + s1.^2 .* v0) ./ ve.^4 ;
          %  conditional variance of deriv. est.


    else ;              %  then do Local Linear regression
          %  main lines from gpnpr.m, via szeg4.m

      flag = (ve < 1) ;
          %  locations with effective sample size < 1
      ve(flag) = ones(sum(flag),1) ;
          %  replace small sample sizes by 1 to avoid 0 divide
          %  no problem below, since gets grayed out
      s1 = conv(bincts(:,1) .* xgrid , kvec) ;
      s1 = s1(k+1:k+ngrid) ;
      s2 = conv(bincts(:,1) .* xgrid.^2 , kvec) ;
      s2 = s2(k+1:k+ngrid) ;
      t0 = conv(bincts(:,2),kvec) ;
      t0 = t0(k+1:k+ngrid) ;
          %  numerator of NW est.

      xbar = conv(bincts(:,1) .* xgrid , kvec) ;
      xbar = xbar(k+1:k+ngrid) ;
          %  Weighted sum of X_i
      xbar = xbar ./ ve ;
          %  Weighted avg of X_i

      t1 = conv(bincts(:,2) .* xgrid , kvec) ;
      t1 = t1(k+1:k+ngrid) ;

      numerd = t1 - t0 .* xbar ;
          %  sum(Y_i * (X_i - xbar)) * W      (weighted cov(Y,X))
      denomd = s2 - 2 * s1 .* xbar + ve .* xbar.^2 ;
          %  sum((X_i - xbar)^2) * W      (weighted var(X))

      flag2 = denomd < (10^(-10) * mean(denomd)) ;
          %  for local linear, need to also flag locations where this
          %  is small, because could have many observaitons piled up
          %  at one location
      ve(flag2) = ones(sum(flag2),1) ;
          %  also reset these, because could have more than 5 piled up
          %  at a single point, but nothing else in window
      flag = flag | flag2 ;
          %  logical "or", which flags both types of locations
          %  to avoid denominator problems

      denomd(flag) = ones(sum(flag),1) ;
          %  this avoids zero divide problems, OK, since grayed out later

      mhat = t0 ./ ve ;
      vd = numerd ./ denomd ;
          %  linear term from local linear fit (which est's deriv).
          %       (sometimes called betahat)

      sig2 = conv(bincts2(:,2),kvec) ;
      sig2 = sig2(k+1:k+ngrid) ;
      sig2 = sig2 ./ ve - mhat.^ 2 ;

      flag2 = sig2 < (10^(-10) * mean(sig2)) ;
      ve(flag2) = ones(sum(flag2),1) ;
          %  also reset these
      flag = flag | flag2 ;
          %  logical "or", which flags both types of locations
          %  to avoid denominator problems
      sig2(flag) = ones(sum(flag),1) ;
          %  this avoids zero divide problems, OK, since grayed out later

      rho2 = vd.^2 .* (denomd ./ (sig2 .* ve)) ;
      sig2res = (1 - rho2) .* sig2 ;
          %  get the residual variance from local linear reg.

      u0 = conv(bincts(:,1) .* sig2res , kvec.^2) ;
      u0 = u0(k+1:k+ngrid) ;
      u1 = conv(bincts(:,1) .* sig2res .* xgrid , kvec.^2) ;
      u1 = u1(k+1:k+ngrid) ;
      u2 = conv(bincts(:,1) .* sig2res .* xgrid.^2 , kvec.^2) ;
      u2 = u2(k+1:k+ngrid) ;
      vv = u2 - 2 * u1 .* xbar + u0 .* xbar.^2 ;
      vv = vv ./ denomd.^2 ;
          %  vector of variances of slope est. (from local linear)
    end ;

  end ;


  %  Get Gaussian quantile, for CI's
  flag = (ve >= 5) ;
          %  locations where effective sample size >= 5
  if sum(flag) > 0 ;
    if isimflag == 0 ;         %  do pt'wise CI's
      gquant = -phiinv((1 - ialpha) / 2) ;
    else ;                     %  do simultaneous CI's
      nxbar = mean(ve(flag)) ;
          %  Crude average effective sample size
      numind = n / nxbar ;
          %  Effective number of independent groups
      beta = (1 - ialpha)^(1/numind) ;
      gquant = -phiinv((1 - beta) / 2) ;
    end ;
  else ;
    gquant = inf ;
  end ;


  mdsurf = [mdsurf vd] ;
  mesurf = [mesurf ve] ;
  mvsurf = [mvsurf vv] ;
  vgq = [vgq gquant] ;

end ;


%  Construct scale space CI surfaces
%
if length(vgq) > 1 ;    %  then have full matrices
  mloci = mdsurf - vec2mat(vgq,ngrid) .* sqrt(mvsurf) ;
          %  Lower confidence (simul.) surface for derivative
  mhici = mdsurf + vec2mat(vgq,ngrid) .* sqrt(mvsurf) ;
          %  Upper confidence (simul.) surface for derivative
else ;       %  have only vectors (since only one h)
  mloci = mdsurf - vgq * sqrt(mvsurf) ;
          %  Lower confidence (simul.) surface for derivative
  mhici = mdsurf + vgq * sqrt(mvsurf) ;
          %  Upper confidence (simul.) surface for derivative
end ;



%  Construct "gray level map", really assignment
%    of pixels to integers, with idea:
%          1 (very dark)    - Deriv. Sig. > 0 
%          2 (darker gray)  - Eff. SS < 5
%          3 (lighter gray) - Eff. SS >= 5, but CI contains 0
%          4 (very light)   - Deriv. Sig. < 0 

mapout = 3 * ones(size(mloci)) ;
            %  default is middle lighter gray

flag = (mloci > 0) ;
            %  matrix of ones where lo ci above 0
ssflag = sum(sum(flag)) ;
if ssflag > 0 ;
  mapout(flag) = ones(ssflag,1) ;
            %  put dark grey where significantly positive
end ;


flag = (mhici < 0) ;
            %  matrix of ones where hi ci below 0
ssflag = sum(sum(flag)) ;
if ssflag > 0 ;
  mapout(flag) = 4 * ones(ssflag,1) ;
            %  put light grey where significantly negative
end ;

flag = (mesurf <= 5) ;
            %  matrix of ones where effective np <= 5 ;
ssflag = sum(sum(flag)) ;
if ssflag > 0 ;

  mapout(flag) = 2 * ones(ssflag,1) ;
            %  put middle darker grey where effective sample size < 5
end ;


%  Transpose for graphics purposes
mapout = mapout' ;         

%  Make small h come out at the bottom
%mapout = flipud(mapout) ;
%  (not needed when "specifying axes backwards", as below)

%  Make plots if no numerical output requested
%
if nargout == 0 ;  %  Then no numerical output, but make a plot
                   %  on the current axes

  %  Set up gray level color map
  comap = [.2, .2, .2; ...
           .35, .35, .35; ...
           .5, .5, .5; ...
           .8, .8, .8] ;

  image([left,right],[log10(hmin),log10(hmax)],mapout) ;
    set(gca,'YDir','normal') ;
    colormap(comap) ;  
    xlabel('x') ;
    ylabel('log10(h)') ;


end ;

