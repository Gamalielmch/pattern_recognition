function [fctres,xeval] = Kerestco(xv,Cn,neval,xeval)

%
%	[fctres,xeval] = Kerestco(xv,Cn,neval,xeval)
%	-------------------------------------------
%
%	evaluation de la ddp par somme de fonctions symetriques cosinus
%	=> Principe de Masry ( Kernel estimation )
%
%	. xv    : realisations du processus (vecteur aleatoire)
%	. Cn    : Puissance du cosinus Cn = 1 (def. 1)
%	. neval : nombre de points pour l'evaluation (def. n)
%	. xeval : domaine d'evaluation (def. [3min-2moy,3max-2moy])
%
%	. fctres : amplitude du modele de la loi
%	. xeval  : domaine de variation de la variable
%

if (nargin < 2), Cn = 1;  end,
if (nargin < 3), neval = length(xv); end,
if (nargin < 4),
xeval = linspace(min(xv),max(xv),neval); end,

hn = 2*sqrt(pi)/sqrt(Cn);

den = 2;

N = length(xv);

for i=1:neval,
  xi = xeval(i);
  fctres(i) = sum(((1 + cos(xv - xi))/den).^Cn);
end,

fctres = fctres/(hn*N);
