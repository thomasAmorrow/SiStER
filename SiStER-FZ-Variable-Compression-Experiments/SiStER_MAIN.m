% SiStER_MAIN.m
%
% Simple Stokes solver with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
% jaolive <at> ldeo.columbia.edu
% March 2011 - Sept 2016

%close all

% INITIALIZATION

% Input File: loads parameter values, model geometry, boundary conditions
if exist('running_from_SiStER_RUN','var')==0
    clear 
    InpFil = 'SiStER_Input_File_fracture_zone_ND';%input('Input file ? ','s');
end
run(InpFil)

tic
% construct grid and initialize marker / node arrays
SiStER_Initialize


%% Initialize weak zones for FZ
if GEOM(2).fzSW==1
%	[etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
	for bb=1:length(GEOM(2).fz)-1
			if (GEOM(2).fzstrong(bb)~=1)
				WKm(xm<GEOM(2).fz(bb)*max(max(X))+GEOM(2).FZseed(1) & ...
				xm>GEOM(2).fz(bb)*max(max(X))-GEOM(2).FZseed(1) & ...
				ym<GEOM(2).top+GEOM(2).FZseed(2) & ...
				ym>GEOM(2).top-(GEOM(2).step(bb+1)) & ...
				im>0)=1;
			end
	end
end

% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0;
t=0;

while time < stoptime
%for t=1:Nt % time loop
    
    % update time
    time=time+dt_m*t0
	t=t+1;
    
    
    % Here we prepare nodal arrays to feed the solver 
    % and get a (vx,vy,P) solution

	%if t<=5
	%	epsIIm(WKm>0)=GEOM(2).weakEPSIIm;
	%	PARAMS.conv_crit_ResL2=1e-15;
	%else
	%	PARAMS.conv_crit_ResL2=1e-10;
	%end
	
    %if (time<GEOM(2).weakTime*3.13e13 && t>5)
		%ep(:)=0;
        %ep(WKm>0)=MAT(2).ecrit;
		%epsIIm(WKm>0)=1e-12;
	%elseif t<=5
		%epsIIm(WKm>0)=1e-12;
		%ep(WKm>0)=MAT(2).ecrit;
		%epsIIm(:)=0;
		%ep(:)=0;
    %end
    
    if time<GEOM(2).weakTime*t0
		
		if t<10		
			epsIIm(WKm>0)=1e-12;
		end
		ep(WKm>0)=2*MAT(2).ecrit;
		%epsIIm(:)=0;
		%ep(:)=0;
    end

	% UPDATE BC
	BC.left(3)=(10*3.13e13-time)/(10*3.13e13)*Uin
	BC.right(3)=-(10*3.13e13-time)/(10*3.13e13)*Uin
        
    SiStER_material_props_on_nodes

    %%% SOLVE STOKES HERE 
    SiStER_flow_solve
    
    % GET STRAIN RATE FROM CURRENT SOLUTION
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
	%epsIIm(GEOM(2).WK_zones>0)=MAT(2).ecrit;
    
    % BUILD UP ELASTIC STRESSES IF ELASTICITY IS ACTIVATED
    if (PARAMS.YNElast==1) 
        SiStER_update_marker_stresses;
    end
    
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    
	% OLD CONDITION
	%if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt % SAVING SELECTED OUTPUT

	% TIME CONDITION
	if (mod(time/356.25/24/3600,dt_out)<=max_step || t==1 || time > stoptime-(max_step*356.25*24*3600))
        disp('SAVING SELECTED VARIABLES') 
        filename=num2str(t);
        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
        save(filename,'X','Y','vx','vy','p','time','xm','ym','etam','rhom','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y')
    end
    
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS,max_step/t0);
    
    % BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1) 
        SiStER_update_ep;
    end

    % ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end
    
    % EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if (PARAMS.Tsolve==1)
        SiStER_thermal_update;
    end

	% ADD THERMAL STRESSES TO ACCOUNT FOR COOLING
	if (PARAMS.YNTstresses==1)
		SiStER_update_thermal_stresses;
		disp(' ')
		disp('THERMAL STRESSES ON')
		disp(' ')
	end

    
    % MARKER UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_update_markers;
    % advect markers in current flow field
    % remove markers if necessary
    % add markers if necessary
    SiStER_update_topography_markers
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    

end

disp('FIN')
toc
    
