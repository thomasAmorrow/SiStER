function [cohes]=SiStER_get_cohesion(im,ep,MAT)
% [cohes]=SiStER_get_cohesion(im,ep,MAT)
% compute cohesion on markers based on ep
%G.Ito 8/2016

Cmax=[MAT(im).Cmax];
Cmin=[MAT(im).Cmin];
epscrit=[MAT(im).ecrit];

%if im==2
%if (GEOM(2).fzstrong(bb)~=1)
% FZ weakness T.Morrow JAN 3 2017 - reseeds ep every time step - is this needed?
%	ep(im>1 & ym<GEOM(2).top+GEOM(2).FZseed(2) & xm<=GEOM(2).fz(bb)*max(xm)+GEOM(2).FZseed(1) ...
%			& xm>GEOM(2).fz(bb)*max(xm)-GEOM(2).FZseed(1))=MAT(2).ecrit;%1/(1/MAT(3).pre_diff+1/MAT(3).pre_disc);
%end

% get cohesion
cohes=max(Cmax+(Cmin-Cmax).*ep./epscrit,Cmin);

return
