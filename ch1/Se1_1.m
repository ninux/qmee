%% 1D FDTD simulation of psi


%% cleanup environment

clear all;
close all;

%% general definitions and constants

h 			= 6.626e-34;		% [J*s]	Plank's constant
hbar		= h/(2*pi);			% [J*s]	reduced Plank's constant
m0			= 9.1e-31;			% [kg] 	free space mas of an electron

meff_Si		= 1.08;				% effective mass of Si
meff_Ge		= 0.067;			% effective mass of Ge
meff_GaAs 	= 0.55;				% effective mass of GaAs

meff		= meff_Si;			% [?] 	effective mass
melec		= meff*m0;			% [?]	mass of an electron
ecoul		= 1.6e-19;			% [C]	charge of an electron
epsz		= 8.85e-9;			% [A*s/V/m]	dielectric of free space
eV2J		= 1.6e-19;			% [-]	energy conversion factor (eV to J)
J2eV		= 1/eV2J;			% [-]	energy conversion factor (J to eV)

%% simulation parameters

NN 			= 400;				% [-]	number of points in the problem space
del_x		= .1e-9;			% [?]	cell size
dt			= 2e-17;			% [s]	time steps
ra			= (0.5*hbar/melec)*(dt/del_x^2);	% <-- must be < 0.15
DX			= del_x*1e9;		% [_]	cell size in nm
XX			= (DX:DX:DX*NN);	% [-]	length in nm for plotting

%% specification of the potential

V			= zeros(1,NN);	% init potential vector

% barrier
for n = NN/2:NN/2+50
%	V(n) = .15*eV2J;
end

% semiconductor condution band
for n = 1 : NN/2
%	V(n) = .1*eV2J;
end

for n = NN/2+1 : NN
%	V(n) = .2*eV2J;
end

% electric field
%for n = 1 : NN
%	V(n) = -(0.2/400)*(n-400)*eV2J;
%end

%% initialize a sine wave in a gaussian envelope

lambda		= -10;			% pulse wavelength
sigma		= 50;			% pulse width
nc			= 200;			% starting position
prl			= zeros(1,NN);	% real part of the state variable
pim			= zeros(1,NN);	% imaginary part of the state variable
ptot		= 0.;			% ?

for n=2:NN-1
	prl(n)	= exp(-1.*((n-nc)/sigma)^2)*cos(2*pi*(n-nc)/lambda);
	pim(n)	= exp(-1.*((n-nc)/sigma)^2)*sin(2*pi*(n-nc)/lambda);
	ptot	= ptot + prl(n)^2 + pim(n)^2;
end

pnorm = sqrt(ptot);			% normalization constant

% normalize and check
ptot = 0.;

for n=1:NN
	prl(n)	= prl(n)/pnorm;
	pim(n)	= pim(n)/pnorm;
	ptot	= ptot + prl(n)^2 + pim(n)^2;
end

fprintf("normalization = %f\n", ptot);	% should have the value 1

T			= 0;
n_step		= 1;

while n_step > 0
	n_step = input("How many time steps -->");

	% this is the core FDTD program
	for m=1:n_step
		T = T+1;
	
		for n=2:NN-1
			prl(n) = prl(n) - ra*(pim(n-1) -2*pim(n) + pim(n+1)) ...
				+ (dt/hbar)*V(n)*pim(n);
		end
		
		for n=2:NN-1
			pim(n) = pim(n) + ra*(prl(n-1) -2*prl(n) + prl(n+1)) ...
				- (dt/hbar)*V(n)*prl(n);
		end
	end
	
	% calculate the expected values
	PE	= 0.;	% potential energy
	
	for n=1:NN
		psi(n) = prl(n) + i*pim(n);
		PE = PE + psi(n)*psi(n)'*V(n);
	end

	fprintf("normalization = %f\n", psi*psi');

	PE = PE*J2eV;	% potential energy

	ke = 0. + j*0.;
	
	for n=2:NN-1
		lap_p = psi(n+1) - 2*psi(n) + psi(n-1);
		ke = ke + lap_p*psi(n)';
	end

	KE = -J2eV*((hbar/del_x)^2/(2*melec))*real(ke);	% kinetic energy

	% plot results
	subplot(2,1,1);
	
	plot(XX, prl, "k");
	hold on;
	plot(XX, pim, "-.k");
	plot(XX, J2eV*V, "--k");
	hold off;
	
	axis([1 DX*NN -.2 .3] );
	
	TT = text(5, .15, sprintf("%7.0f fs", T*dt*1e15));
	set(TT, "fontsize", 12);

	TT = text(5, -.15, sprintf("KE = %5.3f eV", KE));
	set(TT, "fontsize", 12);

	TT = text(25, -.15, sprintf("PE = %5.3f eV", PE));
	set(TT, "fontsize", 12);

	TT = text(25, .13, sprintf("E_t_o_t = %5.3f eV", KE+PE));
	set(TT, "fontsize", 12);

	xlabel("nm");
	set(gca, "fontsize", 12);
	title("Se1-1 - 1D FDTD Simulation");
	grid on;

	T

	saveas(gcf, "Se1_1.png");

end
