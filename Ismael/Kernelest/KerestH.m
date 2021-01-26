function [fctres,xeval] = KerestH(xv,d,neval,xeval)

%
%	[fctres,xeval] = KerestH(xv,d,neval,xeval)
%	-------------------------------------------
%
%	evaluation de la ddp par somme de fonctions symetriques Hilbertienne
%	=> Principe de Masry ( Kernel estimation )
%	
%	. xv    : realisations du processus (vecteur aleatoire)
%	. d     : Definit le volume d'une balle unitaire d = 1 (def. 1)
%	. neval : nombre de points pour l'evaluation (def. n)
%	. xeval : domaine d'evaluation (def. [3min-2moy,3max-2moy])
%
%	. fctres : amplitude du modele de la loi
%	. xeval  : domaine de variation de la variable
%
	
if (nargin < 2), d = 1;  end,
if (nargin < 3), neval = length(xv); end,
if (nargin < 4), 
xeval = linspace(min(xv),max(xv),neval); end,


N = length(xv);

 for j = 1:neval,

  xj = xeval(j);
  if j+1 > neval, xi = xeval(j-1);
   else, xi = xeval(j+1); end,

  fctres(j) = sum(1./(abs(xv - xj).^(2*d) + abs(xv - xi).^(2*d)));
    
 end,


fctres = (fctres*4/(0.5*pi*N*(N-1)*log(N))).^(1/2);
