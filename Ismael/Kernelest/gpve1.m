function [ve,vind] = gpve1(xgrid,v1,v2,iscale,isymm,inorm,iloop) 
% GPVE1, General Purpose Visual Error Measure
%     As developed in Marron and Tsybakov, JASA 1995.
%     These are intended to correspond more to "visual
%     impression of error", than the usual integral norms.
%     For each point on [xgrid,v1] curve, finds "closest
%     point on [xgrid,v2] curve", then summarizes these
%     distances, according to norm choice (with optional 
%     symmetrization).
%   Can use first 3, 4, 5, 6 or 7 arguments.
%     (just using first 3, gives M&T's recommended assymteric
%      VE criteria)
% Inputs:
%     xgrid  - column vector of (common) x-coords of two curves
%                  (norms most sensible, when this is equally spaced)
%     v1     - column vector of y-coords of first curve
%                   (fhat in M&T's recommended visual criterion)
%     v2     - column vector of y-coords of second curve
%                   (f in M&T's recommended visual criterion)
%     iscale - 0 for using distance basedon raw x and y values
%              1 for using distance based on x and y mapped to
%                      [0,1], assuming x min = 0 as in KDE
%              2 for using distance based on x and y mapped to
%                      [0,1], no restrictions, as in regression
%                (2 is the default, when unspecified)
%     isymm  - 0 for a one sided summary
%                (0 is the default, giving M&T's recommended VE
%              1 for a symmetrized summary
%     inorm  - 0 for all 3 norms below (in order 1,2,3), as a col. vector
%              1 for the L1 summary of the distances
%              2 for the L2 summary of the distances
%                (2 is the default, giving both recommendations in M&T)
%              3 for sup norm, which gives Hausdorff type distance
%     iloop  - 0 for using a direct matrix algorithm, which creates
%                    a large nxn square matrix, that could be too large
%                (0 is the default)
%              1 for using a slow looped implementation (could be useful
%                    when input vectors are quite long)
% Outputs:
%     ve     - visual error, according to input parameters
%               (for inorm = 0, col. vector, in order L1, L2, sup)
%     vind   - vector of indices of closest points on v2 curve
%                      to each point on v1 curve
%                  for nice picture, use arrows starting at [xgrid,v1]
%                  and ending at [xgrid(vind),v2(vind)] 
%                  (Note: for isymm = 1, get second set of arrows, only)

%    Copyright (c) J. S. Marron 1997




%  Set parameters and defaults according to number of input arguments
if nargin == 3 ;    %  only 3 arguments input
  iiscale = 2 ;      
          %  Use default: regression type scaling
else ;              %  more than 3 arguments input
  iiscale = iscale ;
          %  Use input scaling type
end ;
if nargin <= 4 ;    %  3 or 4 arguments input
  iisymm = 0 ;
          %  Use default: no symmetrization
else ;              %  more than 4 arguments input
  iisymm = isymm ;
          %  Use input symmetrization type
end ;
if nargin <= 5 ;    %  3, 4 or 5 arguments input
  iinorm = 2 ;
          %  Use default: L2 norm
else ;              %  more than 5 arguments input
  iinorm = inorm ;
          %  Use input norm type
end ;
if nargin <= 6 ;    %  3, 4, 5 or 6 arguments input
  iiloop = 0 ;
          %  Use default: no looping, direct matrix version
else ;              %  more than 6 arguments input
  iiloop = iloop ;
          %  Use input looping option
end ;



%  Do rescaling if needed
%
%  First versions for rescaling
x = xgrid ;
m1 = v1 ;
m2 = v2 ;

%  Rescale x's if needed
if iiscale ~= 0 ;    %  then need to rescale x's to [0,1]
  denom = max(xgrid) - min(xgrid) ;
  if denom ~= 0 ;    %  then have more than x, OK to rescale
    x = (xgrid - min(xgrid)) ./ denom ;
  else ; 
    disp('!!!   Error from gpve1.m, need at least one x to rescale   !!!') ;
  end ;
end ;

%  Rescale v's if needed
if iiscale == 1 ;     %  then do density estimation type rescaling
  denom = max([v1; v2]) ;
  if denom ~= 0 ;     %  then have several y's, OK to rescale
    m1 = v1 ./ denom ;
    m2 = v2 ./ denom ;
  else ;
    disp('!!!   Error from gpve1.m, need different y''s to rescale   !!!') ;
  end ;
elseif iiscale == 2 ;     %  then do regression estimation type rescaling
  mi = min([v1; v2]) ;
  denom = max([v1; v2]) - mi ;
  if denom ~= 0 ;     %  then have several y's, OK to rescale
    m1 = (v1 - mi) ./ denom ;
    m2 = (v2 - mi) ./ denom ;
  else ;
    disp('!!!   Error from gpve1.m, need different y''s to rescale   !!!') ;
  end ;
end ;



%  get vector with indices of nearest point on second curve
%        (to each point on first curve)
%
n = length(xgrid) ;
if iiloop == 1 ;    %  do slow looped version
  vind = [] ;
  for i = 1:n ;
    dx = (x(i) - x).^2 ;
    dy = (m1(i) - m2).^2 ;
    d = dx + dy ;
    [temp,ind] = min(d) ;    %  for each v1, find min'r in v2 direction
    vind = [vind; ind] ;
  end ;
else ;    %  do direct matrix version
  mdx = (vec2mat(x',n) - vec2mat(x,n)).^2 ;
  mdy = (vec2mat(m1',n) - vec2mat(m2,n)).^2 ;
  md = mdx + mdy ;   %  matrix of pairwise distances
  [temp,vind] = min(md) ;    %  for each v1, find min'r in v2 direction
  vind = vind' ;     %  make this a column vector
end ;



%  Calculate norms
%
mdx = (x - x(vind)).^2 ;
mdy = (m1 - m2(vind)).^2 ;
md = mdx + mdy ;

del = abs(x(1) - x(2)) ;
if iinorm == 1 ;    %  then calculate L1 norm
  ve = sum(sqrt(md)) * del ;
elseif iinorm == 2 ;    %  then calculate L2 norm
  ve = sqrt(sum(md) * del) ;
elseif iinorm == 3 ;    %  then calculate sup norm
  ve = max(sqrt(md)) ;
elseif iinorm == 0 ;    %  then calculate all norms
  ve = [sum(sqrt(md)) * del; ...
        sqrt(sum(md) * del); ...
        max(sqrt(md))] ;
end ;  



%  Symmetrize if needed
%
if iisymm == 1 ;    %  then redo with "mirror image", and then combine

  %  get vector with indices of nearest point on second curve
  %        (to each point on first curve)
  %
  if iiloop ~= 1 ;    %  do slow looped version
    vind = [] ;
    for i = 1:n ;
      dx = (x(i) - x).^2 ;
      dy = (m2(i) - m1).^2 ;
           %  note: these are reversed from above
      d = dx + dy ;
      [temp,ind] = min(d) ;    %  for each v1, find min'r in v2 direction
      vind = [vind; ind] ;
    end ;
  else ;    %  do direct matrix version
    mdx = (vec2mat(x',n) - vec2mat(x,n)).^2 ;
    mdy = (vec2mat(m2',n) - vec2mat(m1,n)).^2 ;
           %  note: these are reversed from above
    md = mdx + mdy ;   %  matrix of pairwise distances
    [temp,vind] = min(md) ;    %  for each v1, find min'r in v2 direction
    vind = vind' ;     %  make this a column vector
  end ;

  %  Calculate norms
  %
  mdx = (x - x(vind)).^2 ;
  mdy = (m2 - m1(vind)).^2 ;
           %  note: these are reversed from above
  md = mdx + mdy ;

  if iinorm == 1 ;    %  then calculate symmetrized L1 norm
    ve = ve + sum(sqrt(md)) * del ;
  elseif iinorm == 2 ;    %  then calculate symmetrized L2 norm
    ve = sqrt(ve^2 + sum(md) * del) ;
  elseif iinorm == 3 ;    %  then calculate symmetrized sup norm
    ve = max([ve; sqrt(md)]) ;
  elseif iinorm == 0 ;    %  then calculate symmetrized version of all norms
    ve = [(ve(1) + sum(sqrt(md)) * del); ...
          sqrt(ve(2)^2 + sum(md) * del); ...
          max([ve(3); sqrt(md)])] ;
  end ;  

end ;

