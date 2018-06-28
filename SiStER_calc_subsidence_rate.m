% Calculates the subsidence velocity for time t at model boundaries
% Used to impose a lithospheric subsidence rate based on Turcotte and Schubert 
%  equations

L_Age = GEOM(2).age(1);
R_Age = GEOM(2).age(end);

C = MAT(2).rho0*MAT(2).alpha*(BCtherm.bot(2)-BCtherm.top(2))/ ...
	(MAT(2).rho0-MAT(1).rho0)*sqrt(PARAMS.kref/(MAT(2).rho0*PARAMS.cpref)/pi)*0.5;

L_rate = C*1/(sqrt(time+L_Age));
R_rate = C*1/(sqrt(time+R_Age));

	%BC.left(4)=L_rate;
	%BC.right(4)=R_rate;

%BC.bot(3)=min(L_rate,R_rate);%(L_rate+R_rate)/2; %temp approach MorrowTA 25 Oct 2016
BC.top(3)=(L_rate+R_rate)/2;%min(L_rate,R_rate); %temp approach MorrowTA 25 Oct 2016%(L_rate+R_rate)/2;
			

%BC.top_prof=(R_rate-L_rate)/2*sin(linspace(-pi/2,pi/2,Nx))+(R_rate-L_rate)/2+L_rate; % attempted sine outflow MorrowTA 25 Oct 2016
%BC.bot_prof=(R_rate-L_rate)/2*sin(linspace(-pi/2,pi/2,Nx))+(R_rate-L_rate)/2+L_rate; % attempted sine outflow MorrowTA 25 Oct 2016
