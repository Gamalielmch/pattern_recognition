function iqrange = iqr(data,isort) 
% IQR, Inter Quartile Range
% Inputs:
%     data  - column vector, or matrix of data
%     isort - flag indicating need to sort:
%                0  ===>  Don't need to sort
%                        !!!  DATA ASSUMED TO BE IN INCREASING ORDER  !!!
%                1  ===>  Data unsorted, so first do a sort
%                            (default, when isort not specified)
% Output:
%     iqrange  - interquartile range(s)
%                   scalar when data is a column vector
%                   row vector when data is a matrix
%
% Assumes path can find personal functions:
%    vec2mat
%    cquant

%    Copyright (c) J. S. Marron 1996


%  First decide whether or not to sort, based on number of inputs
if nargin == 1 ;
  iisort = 1 ;
          %  default is to sort, when isort unpsecified    
else ;
  iisort = isort ;
end ;

%  Do sort if needed
if iisort ~= 0 ;    %  then do a sort
  sdata = sort(data) ;
else ;
  sdata = data ;
end ;

%  now get IQR
iqrange = cquant(sdata,.75,0) - cquant(sdata,.25,0) ;
          %  last 0 turns of need to sort

