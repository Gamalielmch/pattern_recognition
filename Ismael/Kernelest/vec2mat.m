function mat = vec2mat(vec,num) 
% VEC2MAT, Turns vectors into matrices (which are duplicates)
% Inputs:
%     vec  - row or column vector
%     num  - number of duplicates of vec to make
% Output:
%     mat  - matrix of duplicates of vec,
%               when vec is an nc x 1 column vector, mat is nc x num
%               when vec is a  1 x nr column vector, mat is num x nr
% Note:
%     To do this for BOTH a row and column vector, it is more 
%     convenient to use "meshgrid"

%    Copyright (c) J. S. Marron 1996, 1997

[nr,nc] = size(vec) ;

if nc == 1 & nr == 1 ;     %  start with scalar, give warning
                           %  and return empty matrix
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!! Caution from vec2mat: not a vector, returning empty matrix !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  mat = [] ;
elseif nc == 1 ;     %  then treat as column vector
  mat = vec * ones(1,num) ;
elseif nr == 1 ; %  then treat as row vector
  mat = ones(num,1) * vec ;
else ;            %  then not a vector, give warning
                  %  and return empty matrix
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!! Caution from vec2mat: not a vector, returning empty matrix !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  mat = [] ;
end ;

