function [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
% [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
% sets the advection time step
% J.-A. Olive, November 2014

dt_m=PARAMS.fracCFL*0.5*min(min(dx),min(dy))./max(max(max(abs(vx))),max(max(abs(vy))))

%if dt_m>25e3*356.25*24*3600;
%	dt_m=25e3*356.25*24*3600;
%end
