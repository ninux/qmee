##############################################################################
##
## Copyright (C) 2023, Ervin Mazlagic
##
## This file is part of the personal solutions manual for "Quantum Mechanics
## for Electrical Engineers" by Dennis M. Sullivan (ISBN 978-0-470-87409-7),
## see "github.com/ninux/qmee".
##
##############################################################################

## Chapter 1 - Exercise 2-5
## 1D FDTD simulation of psi

# cleanup environment
clear all;
close all;

% measure time for complete script to run
tic();

# general definitions and constants
h           = 6.626e-34;    # [J*s]     Plank's constant
hbar        = h/(2*pi);     # [J*s]     reduced Plank's constant
m0          = 9.1e-31;      # [kg]      free space mass of electron (rest mass)

meff_Si     = 1.08;         # [-]       effective mass factor of Si
meff_Ge     = 0.067;        # [-]       effective mass factor of Ge
meff_GaAs   = 0.55;         # [-]       effective mass factor of GaAs

meff        = meff_Si;      # [m]       effective mass
melec       = meff*m0;      # [?]       mass of an electron
ecoul       = 1.6e-19;      # [C]       charge of an electron
epsz        = 8.85e-9;      # [A*s/V/m] dielectric of free space
eV2J        = 1.6e-19;      # [-]       energy conversion factor (eV to J)
J2eV        = 1/eV2J;       # [-]       energy conversion factor (J to eV)

# simulation parameters
NN          = 400;          # [-]       number of points in the problem space
del_x       = .1e-9;        # [?]       cell size
dt          = 2e-17;        # [s]       time steps
DX          = del_x*1e9;    # [_]       cell size in nm
XX          = (DX:DX:DX*NN);# [-]       length in nm for plotting
ra          = (0.5*hbar/melec)*(dt/del_x^2);    # must be < 0.15 for stability

# check stability
if (ra < 0.15)
    # stable
elseif
    # unstable
    printf("WARNING: Simulation not stable, check parameters ...");
endif

# specification of the potential
# original code with loops
#{
V = zeros(1,NN);  # init potential vector

## potential definitions

# constant barrier band in center
for n = NN/2:NN/2+50
    V(n) = .15*eV2J;
end

# semiconductor condution band
for n = 1 : NN/2
    V(n) = .1*eV2J;
end

for n = NN/2+1 : NN
    V(n) = .2*eV2J;
end

# electric field
for n = 1 : NN
    V(n) = -(0.2/400)*(n-400)*eV2J;
end
#}
# new code with potential selection and without loops
V = zeros(1,NN);    # initialize potential vector with 0

# potential options
#{
    A = zero potential
    B = constant barrier band in center
    C = semiconductor conduction band
    D = one-sided finite wall
    E = linear electric field
#}
potentialOption = "A";

if (potentialOption == "A")
    # zero potential
elseif (potentialOption == "B")
    # constant barrier band in center
    bw = 50;                    # barrier width
    vb = 0.15*eV2J;             # barrier potential
    V(NN/2:NN/2+bw) = vb;
elseif (potentialOption == "C")
    # semicondctor conduction band
    vb1 = 0.1*eV2J;             # 1st potential
    vb2 = 0.2*eV2J;             # 2nd potential
    V(1:NN/2) = vb1;
    V(NN/2+1:NN) = vb2;
elseif (potentialOption == "D")
    # one-sided finite wall
    vb = 0.2*eV2J;              # potential
    V(NN/2:NN) = vb;
elseif (potentialOption == "E")
    # linear electric field
    vbm = (0.2/NN)*eV2J;        # max. potential
    V(1:NN) = ([NN:-1:1]')*vbm;
elseif
    printf("ERROR: invalid potential option, aborting script ...\n");
    return
endif

#initialize a sine wave in a gaussian envelope
lambda  = 10;           # pulse wavelength
sigma   = 50;           # pulse width
nc      = 200;          # starting position
prl     = zeros(1,NN);  # real part of the state variable
pim     = zeros(1,NN);  # imaginary part of the state variable
ptot    = 0.;           # total energy

# original code with loops
#{
for n=2:NN-1
    prl(n)  = exp(-1.*((n-nc)/sigma)^2)*cos(2*pi*(n-nc)/lambda);
    pim(n)  = exp(-1.*((n-nc)/sigma)^2)*sin(2*pi*(n-nc)/lambda);
    ptot    = ptot + prl(n)^2 + pim(n)^2;
end
#}
# new code without loops
prl(2:NN-1) = exp(-1.*(([2:NN-1]-nc)./sigma).^2) ...
            .* cos((2*pi).*([2:NN-1]-nc)./lambda);
pim(2:NN-1) = exp(-1.*(([2:NN-1]-nc)./sigma).^2) ...
            .* sin((2*pi).*([2:NN-1]-nc)./lambda);
ptot = sum(prl.^2) + sum(pim.^2);

pnorm = sqrt(ptot); # normalization constant

# normalize and check
ptot = 0.;

# original code with loops
#{
for n=1:NN
    prl(n)  = prl(n)/pnorm;
    pim(n)  = pim(n)/pnorm;
    ptot    = ptot + prl(n)^2 + pim(n)^2;
end
#}
# new code without loops
prl(1:NN) = prl(1:NN)./pnorm;
pim(1:NN) = pim(1:NN)./pnorm;
ptot = sum(prl(1:NN).^2) + sum(pim(1:NN).^2);

# check normalization for total energy
    norm_limit = (1e-3)*[-1, 1] + 1;
    if (ptot > norm_limit(1) && ptot < norm_limit(2))
        # normalization ok
        fprintf("INFO: Normalization for energy ok (n = %f, error = %2.3f %%).\n",
                ptot, (ptot-1)*1e2);
    else
        # normalization not ok
        fprintf("WARNING: Normalization for energy not ok (n = %f, error = %2.3f %%).\n",
                ptot, (ptot-1)*1e2);
    endif

# parameters for the FDTD simulation
T       = 0;    # absolute time in steps
n_step  = 1;    # steps per loop
count   = 0;    # number of loops
nos     = 200;  # number of timestepts between plots (equidistant)
nop     = 3;    # number of plots to be made

# perform simulation
while nop > count

    # this is the core FDTD program
    for m=1:n_step
        T = T+1;

        #{
        for n=2:NN-1
            prl(n) = prl(n) ...
                   - ra*(pim(n-1) -2*pim(n) + pim(n+1)) ...
                   + (dt/hbar)*V(n)*pim(n);
        end

        for n=2:NN-1
            pim(n) = pim(n) ...
                   + ra*(prl(n-1) -2*prl(n) + prl(n+1)) ...
                   - (dt/hbar)*V(n)*prl(n);
        end
        #}
        # new code without loops
        rng = [2:NN-1];

        prl(rng) = prl(rng) ...
                 - ra.*(pim(rng-1) - 2.*pim(rng) + pim(rng+1)) ...
                 + (dt/hbar).*V(rng).*pim(rng);

        pim(rng) = pim(rng) ...
                 + ra.*(prl(rng-1) - 2.*prl(rng) + prl(rng+1)) ...
                 - (dt/hbar).*V(rng).*prl(rng);
    end

    # calculate the expected values
    PE = 0.;    # initialize potential energy

    # original code with loops
    #{
    for n=1:NN
        psi(n) = prl(n) + i*pim(n);
        PE     = PE + psi(n)*psi(n)'*V(n);
    end
    #}
    # new code without loops
    psi(1:NN) = prl(1:NN) + i.*pim(1:NN);
    PE = sum(psi(1:NN).*conj(psi(1:NN)).*V(1:NN));

    # check normalization
    check_norm = psi*psi';
    norm_limit = (1e-3)*[-1, 1] + 1;
    if (check_norm > norm_limit(1) && check_norm < norm_limit(2))
        # normalization ok
        fprintf("INFO: Normalization for psi ok (n = %f, error = %2.3f %%).\n",
                check_norm, (check_norm-1)*1e2);
    else
        # normalization not ok
        fprintf("WARNING: Normalization for psi not ok (n = %f, error = %2.3f %%).\n",
                check_norm, (check_norm-1)*1e2);
    endif

    PE = PE*J2eV;    # potential energy

    ke = 0. + j*0.;

    # original code with loops
    #{
    for n=2:NN-1
        lap_p = psi(n+1) - 2*psi(n) + psi(n-1);
        ke = ke + lap_p*psi(n)';
    end
    #}
    # new code without loops
    rng = [2:NN-1];
    ke = sum((psi(rng+1)-2.*psi(rng)+psi(rng-1)).*conj(psi(rng)));

    KE = -J2eV*((hbar/del_x)^2/(2*melec))*real(ke);    # kinetic energy

    # plot results
    figure(count+1);
    plot(XX, prl, "k");
    hold on;
    plot(XX, pim, "-.k");
    plot(XX, J2eV*V, "--k");
    hold off;

    axis([1 DX*NN -.2 .2] );

    # insert plot annotations
    text( 5,  .15, sprintf("%7.0f fs",           T*dt*1e15)); # show time in plot
    text( 5, -.15, sprintf("KE = %5.3f eV",      KE));        # show KE in plot
    text(30, -.15, sprintf("PE = %5.3f eV",      PE));        # show PE in plot
    text(30,  .15, sprintf("E_t_o_t = %5.3f eV", KE+PE));     # show Etot in plot

    xlabel("nm");

    set(gca, "fontsize", 10);
    if (count == 0)
        title("Exercise 2-5a - 1D FDTD Simulation");
    else
        title(" ");
    endif
    grid on;

    fig = gcf();
    set(fig, "PaperUnits", "centimeters");
    x_width = 21;
    y_width = 7;
    set(fig, "PaperPosition", [0 0 x_width y_width]);
    print(fig, strcat("ex_2_5_a_", num2str(count), ".png"), "-dpng");

    count++;
    n_step = nos;
end

# show time
elapsed_time = toc();
printf("INFO: Total time for simulation is %f s\n", elapsed_time);
