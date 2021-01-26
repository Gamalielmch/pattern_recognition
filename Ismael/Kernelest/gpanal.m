function makeplot = gpanal(data,icolor,vrange,nhp,psfname,idatap) 
% GPANAL, General Purpose data ANALysis
%     Makes a single page with family analysis, and SiZer 1
%     for a quick first look at new data sets. 
%   Can use first 1, 2, 3, 4 or 5 arguments.
% Inputs:
%     data    - either n x 1 column vector of density estimation data
%                  or n x 2 matrix of regression data:
%                            X's in first column,  Y's in second
%     icolor  - 0  use gray levels
%               1 (or not specified) do red, blue purple
%     vrange  - 0 (or not specified) use the min and max of the data
%               2 vector:  use minx as vrange(1) and maxx as vrange(2)
%               3 vector:  use minx as vrange(1) and maxx as vrange(2)
%                               and number of grid points as vrange(3)
%     nhp     - 0 use default value of 11
%               otherwise, value is number in family, and rows in SiZer
%     psfname - string with prefix for name of postscript file
%                    (will automatically add ".ps")
%                    (will be a color postscript file when icolor = 1)
%     idatap  - 0  don't add data to family plot
%               1  (or not specified) add data to family plot, as
%                    jitter plot (density estimation) or scatterplot
%
% Outputs:
%     Only graphics, in current Figure, unless a postscript file is created
%
% Assumes path can find personal functions:
%    bwsjpib.m
%    bwrswb.m
%    gpkde.m
%    gpnpr.m
%    gpsz1.m

%    Copyright (c) J. S. Marron 1997, 1998




%  Set parameters and defaults according to number of input arguments
%
if nargin == 1 ;    %  only 1 argument input, use default icolor
  iicolor = 1 ;
else ;              %  icolor was specified, use that
  iicolor = icolor ;
end ;

if nargin <= 2 ;    %  at most 2 arguments input, use default range
  ivrange = 0 ;
else ;              %  vrange was specified, use that
  ivrange = vrange ;
end ;

if nargin <= 3 ;    %  at most 3 arguments input, use default nhp
  inhp = 11 ;
else ;              %  vrange was specified, use that
  if nhp == 0 ;
    inhp = 11 ;
  else ;
    inhp = nhp ;
  end ;
end ;

if nargin <= 4 ;    %  at most 4 arguments input, create no postscript file
  ipsfname = [] ;    %  Default
else ;
  ipsfname = psfname ;    %  Create Postscript file, with given prefix
end ;

if nargin <= 5 ;    %  at most 5 arguments input, so plot data
  iidatap = 1 ;    %  Default
else ;
  iidatap = idatap ;    %  Plot or not, depending on input value
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


disp('    Working on family') ;
subplot(2,1,1) ;
  if length(ivrange) == 1 ;
    mind = min(xdat) ;
    maxd = max(xdat) ;
    ngrid = 401 ;
  else ;
    mind = vrange(1) ;
    maxd = vrange(2) ;
    if length(ivrange) == 3 ;
      ngrid = vrange(3) ;
    else ;
      ngrid = 401 ;
    end ;
  end ;
  centd = mean([mind;maxd]) ;

  %  Set h grid stuff, as in SiZer1
  range = maxd - mind ;
  binw = range / (ngrid - 1) ;
  hmin = 2 * binw ;
  hmax = range ;
  vh = logspace(log10(hmin),log10(hmax),inhp) ;

  if idatyp == 1 ;    %  doing density estimation
    hcent = bwsjpib(data) ;
  else ;              %  doing regression
    hcent = bwrswb(data) ;
  end ;
  if hcent < min(vh) ;    %  if data based h below range
    disp('!!!   Warning from gpanal: databased h below this range   !!!') ;
    hcflag = 0 ;
  elseif hcent > max(vh) ;    %  if databased h above this range
    disp('!!!   Warning from gpanal: databased h above this range   !!!') ;
    hcflag = 0 ;
  else ;     %  if data based in range, then get its nearest index
    [temp, hcind] = min(abs(log10(hcent) - log10(vh))) ;
    hcind = inhp + 1 - hcind ;
                        %  since handles appear in reverse order in plot
    hcflag = 1 ; 
  end ;
  if idatyp == 1 ;    %  doing density estimation
    [msmoo, xgrid] = gpkde(data,vh,[mind; maxd; ngrid]) ;
    bottom = 0 ;
  else ;              %  doing regression
    [msmoo, xgrid] = gpnpr(data,vh,[mind; maxd; ngrid]) ;
    bottom = min(min(msmoo)) ;
  end ;
  top = max(max(msmoo)) ;
  range = top - bottom ;
  bottom = bottom - 0.05 * range ;
  top = top + 0.05 * range ;

           %  Next lines from gpfam.m
  plot(xgrid,msmoo,'c') ;
          %  Plot most curves in cyan
    vcurvh = get(gca,'Children') ;
          %  Vector of handles for curves
    if hcflag == 1 ;
      set(vcurvh(hcind),'LineWidth',2) ;
          %  use fatter line width for data based choice
      set(vcurvh(hcind),'Color',[1 0 0]) ;
          %  use red color for one in middle
    end ;
    vax = axis ;

    vax([1,2]) = [mind,maxd] ;
    axis(vax) ;
    title(['Family Plot, from gpanal, ' date]) ;

    if iidatap == 1 ;    %  then overlay a plot of the raw data
      hold on ;
        if idatyp == 1 ;    %  doing density estimation
          yrand = vax(3) + (0.8 + 0.1 * rand(size(data))) ...
                                                 * (vax(4) - vax(3)) ;
               %  y coords for jittering
          plot(data,yrand,'g.') ;
        else ;              %  doing regression
          plot(data(:,1),data(:,2),'g.') ;
        end ;
      hold off ;
    end ;



disp('    Working on SiZer') ;
subplot(2,1,2) ;
  gpsz1(data,[mind,maxd,ngrid]) ;
    title(['SiZer 1 Plot, from gpanal, ' date]) ;
          %  next add marking lines
    hold on ;
      if hcflag == 1 ;
        plot([mind; maxd], ones(2,1)*log10(hcent),'-k') ;
      end ;
      plot(centd + 2*vh, log10(vh),':k') ;
      plot(centd - 2*vh, log10(vh),':k') ;
    hold off ;


if iicolor ~= 0 ;     %  Then go for nice colors in sizer
  %  Set up colorful color map
  cocomap = [0,    0,   1; ...
            .35, .35, .35; ...
            .5,    0,  .5;
             1,    0,   0] ; ...
  colormap(cocomap) ;
end ;


if length(ipsfname) > 0 ;     %  then create postscript file with given name
    orient tall ;

  if iicolor ~= 0 ;     %  then make color postscript
    eval(['print -dpsc \matlab\steve\ps\' ipsfname '.ps']) ;
  else ;                %  then make black and white
    eval(['print -dps \matlab\steve\ps\' ipsfname '.ps']) ;
  end ;

end ;

