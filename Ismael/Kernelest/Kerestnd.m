function [fctres,xeval] = Kerestnd(xv,hn,neval,xeval)

%
%	[fctres,xeval] = Kerestnd(xv,hn,neval,xeval)
%	-------------------------------------------
%
%	evaluation de la ddp par somme de gaussiennes normalisees
%	=> Principe de Masry ( Kernel estimation )
%
%	. xv    : realisations du processus (vecteur aleatoire)
%	. hn    : bandwidth hn = 0.1 (def. 0.1)
%	. neval : nombre de points pour l'evaluation (def. n)
%	. xeval : domaine d'evaluation (def. [3min-2moy,3max-2moy])
%
%	. fctres : amplitude du modele de la loi
%	. xeval  : domaine de variation de la variable
%

if (nargin < 2), hn = 0.1;  end,
if (nargin < 3), neval = length(xv); end,
if (nargin < 4), 
xeval = linspace(min(xv),max(xv),neval); end,

den = 2*hn^2;
N = length(xv);

for i=1:neval,
  xi = xeval(i);
  fctres(i) = sum(exp(-(xv - xi).^2/den).*exp(-(xv + xi).^2/den));
end,

fctres = fctres/(sqrt(2*pi)*hn*N);
