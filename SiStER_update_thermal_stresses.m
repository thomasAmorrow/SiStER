% SiStER_Update_thermal_stresses
%
% 	adds Sxx due to thermal stresses at temperatures below 300C (Turcotte, 1974)
%
%
% T.Morrow 04 Jul 2017

sxxmTOLD=sxxmT;

% calculate thermal stresses with new temperature solution, remove old 
%sxxmT=MAT(2).alpha*MAT(2).G.*(600/T0-Tm)-sxxmTOLD;
% why remove old stress? 
dsxxmT=MAT(2).alpha*MAT(2).G.*(Tm_old-Tm);

% remove any sxxT > 600 deg
dsxxmT=(dsxxmT).*(Tm<600/T0);
dsxxmT=dsxxmT.*(im>0.5);

% add to stresses
sxxm=sxxm+dsxxmT;
