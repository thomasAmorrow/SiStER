% SiStER_Input_File


% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=1000; % max number of time iterations turned off currently
stoptime=10*3600*24*365.25*1e6 % myr
dt_out=100e3; % output every dt_out years
max_step=25e3; % maximum time step (assuming CFL criterion not met)
% OLD - output files every "dt_out" iterations

H0=100e3;
eta0=1e20;
rho0=3300;
V0=1e-7;
SS0=eta0*V0/H0;
t0=H0/V0;
T0=1300;
StNum=rho0*9.8*H0^2/eta0/V0;
R=8.314; %need this earlier
Cp0=1000;
K0=3;

% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=600e3/H0;
ysize=100e3/H0;
% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y
GRID.dx(1)=2000/H0;
GRID.x(1)=280e3/H0;
GRID.dx(2)=1000/H0;
GRID.x(2)=320e3/H0;
GRID.dx(3)=2000/H0;
GRID.x(3)=600e3/H0;
%GRID.dx(4)=1000;
%GRID.x(4)=400e3;
%GRID.dx(5)=4000;
%GRID.x(5)=600e3;
GRID.dy(1)=500/H0;
GRID.y(1)=60e3/H0;
GRID.dy(2)=2000/H0;
GRID.y(2)=100e3/H0;
%GRID.dy(3)=2000;
%GRID.y(3)=150e3;
%GRID.dy(4)=1000;
%GRID.y(4)=100e3;
%GRID.dy(5)=4000;
%GRID.y(5)=150e3;


% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=6; % number of markers in the smallest quadrant
Mquad_crit=4; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphase=3; % number of phases

% phase 1
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=20e3/H0;

% phase 2
GEOM(2).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=20e3/H0;
GEOM(2).bot=80e3/H0;

% phase 3
GEOM(3).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(3).top=80e3/H0;
GEOM(3).bot=100e3/H0;

% phase 4
GEOM(4).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(4).top=150e3/H0;
GEOM(4).bot=150e3/H0;

% FZ GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% triggers temp structure and bathymetric step for n oceanic lithosphere
% blocks
GEOM(2).fzSW=1;                                 % fracture zone switch
GEOM(2).fz=[0.4833 0.5167 1.0];           % fraction of domain for each lithosphere block (n entries)
GEOM(2).age=3.154e13*[10 2 0.01]./t0;  % age of each block [Ma]*convert
GEOM(2).step=2*3.1e-5*3300*1280/2300*sqrt(3*abs(GEOM(2).age-max(GEOM(2).age))/(3300*1000)/pi)/H0;             % bathymetric step size
GEOM(2).FZseed=1e3*[2 20]./H0;   %weak seed dimensions [width depth] km

GEOM(2).fzstrong=[0 0 0];

GEOM(2).weakEPSIIm=1e-11;
GEOM(2).weakTime=0.01*3.13e13/t0;        %ma
% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1
MAT(1).phase=1;
% density parameters
MAT(1).rho0=0.01/rho0;     
MAT(1).alpha=0;
% elasticity 
MAT(1).G=1e18/SS0;
% diffusion creep parameters
MAT(1).pre_diff=.5/1e18*eta0;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=.5/1e18*eta0;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).Cmax=44e6/SS0;
MAT(1).Cmin=0.01e6/SS0;
MAT(1).ecrit=0.1;


% phase 2
MAT(2).phase=2;
% density parameters
%MAT(2).rho0=1e3;
%MAT(2).alpha=0;
MAT(2).rho0=3300/rho0;
MAT(2).alpha=3.1e-5*T0;
% elasticity 
MAT(2).G=30e9/SS0;
% diffusion creep parameters
%MAT(2).pre_diff=.5/1e18;
%MAT(2).Ediff=0;
%MAT(2).ndiff=1;
% dislocation creep parameters
%MAT(2).pre_disc=.5/1e18;
%MAT(2).Edisc=0;
%MAT(2).ndisc=1;
% diffusion creep parameters
MAT(2).pre_diff=0.5/1e40;
MAT(2).Ediff=0;
MAT(2).ndiff=1;
% dislocation creep parameters
MAT(2).pre_disc=1e-3*eta0;
MAT(2).Edisc=540000/T0;
MAT(2).ndisc=1;
% plasticity
MAT(2).mu=0.001;
MAT(2).Cmax=44e6/SS0;
MAT(2).Cmin=0.01e6/SS0;
MAT(2).ecrit=0.5;


% phase 3
MAT(3).phase=3;
% density parameters
MAT(3).rho0=3330/rho0;
MAT(3).alpha=3.1e-5*T0;
% elasticity 
MAT(3).G=30e9/SS0;
% diffusion creep parameters
MAT(3).pre_diff=0.5/1e40;
MAT(3).Ediff=0;
MAT(3).ndiff=1;
% dislocation creep parameters
MAT(3).pre_disc=1e-3*eta0;
MAT(3).Edisc=540000/T0;
MAT(3).ndisc=1;
% plasticity
MAT(3).mu=0.01;
MAT(3).Cmax=44e6/SS0;
MAT(3).Cmin=4e6/SS0;
MAT(3).ecrit=0.5;

% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=1; % elasticity on (1) or off (0)
PARAMS.YNPlas=1; % plasticity on (1) or off (0)
PARAMS.epsII_from_stress=1; % get strain rate from stresses (1, default) or from velocity field (0)
PARAMS.tau_heal=1e80/t0; % healing time for plasticity (s)
PARAMS.gx=0; % gravity along x
PARAMS.gy=StNum; % gravity along y
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=1e28/eta0; % maximum viscosity
PARAMS.etamin=1e16/eta0; % minimum viscosity
PARAMS.Tsolve=1; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.YNTstresses=0; % yes (1) or no (0) add thermal stresses below 300 C
PARAMS.a0=0;
PARAMS.a1=0;
PARAMS.a2=0;
PARAMS.a3=1000/(30e3)^3;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=2/T0;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=1; 
PARAMS.kref=1/(H0*V0*rho0*Cp0);
PARAMS.cpref=1;


% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1e4; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-9; % diffusivity of topography (m^2/s)

% Solver iterations
PARAMS.Npicard_min=4; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=50; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-10;
PARAMS.pitswitch=30; % number of Picard iterations at which the solver switches to quasi-Newton
% convergence is assumed if a fraction of the domain smaller 
% than PARAMS.conv_crit2 * domain size has more than relative
% change of PARAMS.conv_crit1 in strain rate and velocity field


% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=0; % pressure in the top-left corner of the domain (anchor point)


% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant, 1=triangular profile, set up for top and bottom only)
% 3/ value of normal velocity (or max. value along the profile, if 2/ is 1)
Uin=0.5e-5;%-6.634e-6;

BC.top=[1 0 0];%1.0563e-10];
BC.bot=[1 3 0];%4e-5]; %-1.0563e-10];
BC.left=[1 0 Uin]; %-3.1688e-10];
BC.right=[1 0 -Uin]; %3.1688e-10];
PARAMS.BalanceStickyLayer=1; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% Mixed BC; Only matters if BC(1,2)=2
BC.left_w  = [1 0 0]; % Water
BC.left_l  = [1 0 0.5*Uin]; % Lithospheere
BC.left_m  = [1 0 0]; % Mantle

BC.right_w = [1 0 0]; % Water
BC.right_l = [1 0 -0.5*Uin]; % Lithosphere
BC.right_m = [1 0 0]; % Mantle

BC.sliding_lith  = [1 10]; % [a b]; if a=1, applies lith. BC to b elements
                          % above and below lithosphere. This is to allow
                          % lithosphere sliders to move without sliders in
                          % the water and mantle 
BC.mech_lith = 0;         % Use viscosity to define lithosphere vs layer 2 
BC.min_lith_visc = 1e20; % Viscosities above this are considered to be
                           % lithosphere for mixed BC purposes
                           
% BC.sticky_bot_left = 20e3;
% BC.sticky_bot_right = 20e3;
% BC.lith_bot_left = 250e3;
% BC.lith_bot_right = 250e3;
% 
% BC.p_lith_bot=8.200e9;
% BC.p_lith_top=0;

BC.INJ.on=0;

% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 0=Neumann)
% 2/ value
BCtherm.top=[1 2/T0];
BCtherm.bot=[1 1300/T0];
BCtherm.left=[2 0];
BCtherm.right=[2 0];
