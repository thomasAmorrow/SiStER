% SiStER_FZisoset

% sets initial isostatic offset for FZ models
%       T Morrow 09 May 2017

ZZ=0:0.001:(ysize-GEOM(2).top);

TT=zeros(size(ZZ))';

for kk=1:length(GEOM(2).age);
TT(:,kk)=PARAMS.T0+(BCtherm.bot(2)-PARAMS.T0).* ...
    erf((ZZ)./(2*sqrt(PARAMS.kref/(MAT(2).rho0*PARAMS.cpref).*(GEOM(2).age(kk)))));
end

RHO=[MAT(2).rho0].*(1-[MAT(2).alpha].*(TT-PARAMS.T0));

for  kk=1:length(GEOM(2).age);
    DEFICIT(kk)=sum(RHO(:,1))-sum(RHO(:,kk));
end

H_add=0.001*DEFICIT./min(RHO(:,1));

GEOM(2).step=H_add;
