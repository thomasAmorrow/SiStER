% SiStER Initialize

% modified for FZ MorrowTA 25 Oct 2016

% construct staggered grids
[X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID);

% initialize marker arrays and positions
[xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);

% locate markers with respect to grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);

% balance isostasy
SiStER_FZisoset

% assign marker phases
[im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);
% initialize marker plastic strain (to zero) and strain rate (to one)
ep=zeros(size(xm));
epNH=ep;
epsIIm=ones(size(xm));
WKm=zeros(size(xm));

% FZ weak seeding ----------------------------------------- % 

if GEOM(2).fzSW==1
    for bb=1:length(GEOM(2).fz)-1
		if (GEOM(2).fzstrong(bb)~=1)
	%im(im>1 & ym>=GEOM(2).top-10e3 & ym<GEOM(2).top+GEOM(2).FZseed(2) & xm<=GEOM(2).fz(bb)*max(xm)+GEOM(2).FZseed(1) ...
        %    & xm>GEOM(2).fz(bb)*max(xm)-GEOM(2).FZseed(1))=3;
			%ep(im>1 & ym<GEOM(2).top+GEOM(2).FZseed(2) & xm<=GEOM(2).fz(bb)*max(xm)+GEOM(2).FZseed(1) ...
			%& xm>GEOM(2).fz(bb)*max(xm)-GEOM(2).FZseed(1))=40*MAT(2).ecrit;%1/(1/MAT(3).pre_diff+1/MAT(3).pre_disc);
		end
	%MAT(2).ecrit;
			%1/(1/MAT(3).pre_diff+1/MAT(3).pre_disc);        
	%ep(im==2 & ym>=GEOM(2).top-5000 & ym<GEOM(2).top+GEOM(2).FZseed(2) & xm<=0.5*max(xm)+GEOM(2).FZseed(1) ...
        %    & xm>0.5*max(xm)-GEOM(2).FZseed(1))=MAT(2).ecrit;
        
%GEOM(2).fz*max(xm)+GEOM(2).FZseed(bb) & xm>GEOM(2).fz*max(xm)+GEOM(2).FZseed(bb))=MAT(2).ecrit;
    end
end

% --------------------------------------------------------- %

% initialize marker stresses
sxxm=zeros(size(xm));
sxym=sxxm;

% initialize marker index (a unique number to identify and track each marker)
idm=1:length(xm);

% initialize temperature structure on nodes
%T=PARAMS.a0+PARAMS.a1*Y+PARAMS.a2*Y.^2+PARAMS.a3*Y.^3;
%T=T+PARAMS.amp*sin(2*pi*X/PARAMS.lam);

% FZ geometry modification -------------------------------- %
SiStER_FZisoset

T=0;
for kk=1:length(GEOM(2).fz)
    if kk==1
        T=T+PARAMS.T0+(BCtherm.bot(2)-PARAMS.T0)*erf((Y-GEOM(2).top-GEOM(2).step(kk))/(2*sqrt(PARAMS.kref/(MAT(2).rho0*PARAMS.cpref)*(GEOM(2).age(kk))))).*(X<=GEOM(2).fz(kk)*max(max(X)));
    else
        T=T+PARAMS.T0+(BCtherm.bot(2)-PARAMS.T0)*erf((Y-GEOM(2).top-GEOM(2).step(kk))/(2*sqrt(PARAMS.kref/(MAT(2).rho0*PARAMS.cpref)*(GEOM(2).age(kk))))).*((X<=GEOM(2).fz(kk)*max(max(X)) & X>GEOM(2).fz(kk-1)*max(max(X))));
    end
end

% Matches profile to boundary condition
T(T>PARAMS.T0)=T(T>PARAMS.T0)-PARAMS.T0;
T(end,:)=T(end,:)+2*PARAMS.T0;
% --------------------------------------------------------- %

if PARAMS.ynTreset==1 % reset T=T0 in top layer
    T(T<PARAMS.T0)=PARAMS.T0;
end
% pass initial nodal T to markers
[Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm0=Tm;

% initialize nodal strain rate and other useful arrays
EXX=zeros(size(X));
EXY=zeros(size(X));
vx=zeros(size(X));
vy=zeros(size(X));
v=vx;
p=zeros(size(EXX));%1e12*ones(size(EXX));  %initialize to be high so plasticity doesnt activate at t=1, pit=1;
etan_new=zeros(Ny,Nx);
%-------------------------------------------------------------------------
% initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
%-------------------------------------------------------------------------
if (exist('dt_m','var')==0);
    dt_m=1e2/t0;
end

% initialize marker chain to track base of layer 1 (sticky layer)
%Ntopo=PARAMS.Ntopo_markers;
%topo_x=linspace(0,xsize,Ntopo);
%topo_y=GEOM(1).bot*ones(size(topo_x));
%topo_marker_spacing=mean(diff(topo_x)); % initial mean spacing of topography markers

% % FZ TOPO MARKER SETUP ----------------------------------- %
Ntopo=1e4;
topo_x=linspace(0,xsize,Ntopo);
topo_y=GEOM(1).bot*ones(size(topo_x));
for bb=2:length(GEOM(2).fz);
        topo_y(topo_x<=GEOM(2).fz(bb)*xsize & topo_x>GEOM(2).fz(bb-1)*xsize)=topo_y(topo_x<=GEOM(2).fz(bb)*xsize & topo_x>GEOM(2).fz(bb-1)*xsize)-GEOM(2).step(bb);
end
topo_marker_spacing=mean(diff(topo_x)); % initial mean spacing of topographic markers
% --------------------------------------------------------- %

