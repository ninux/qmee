##############################################################################
##
## Copyright (C) 2023, Ervin Mazlagic
##
## This file is part of the personal solutions manual for "Quantum Mechanics
## for Electrical Engineers" by Dennis M. Sullivan (ISBN 978-0-470-87409-7),
## see "github.com/ninux/qmee".
##
##############################################################################

## Chapter 2 - Exercise 2-4
## Stationary States

# cleanup environment
clear all;
close all;

# general definitions and constants
h           = 6.626e-34;    # [J*s]     Plank's constant
hbar        = h/(2*pi);     # [J*s]     reduced Plank's constant
m0          = 9.1e-31;      # [kg]      free space mass of electron (rest mass)
A           = 10e-10;       # [m]       Angstrom

meff_vac    = 1.0;          # [-]       effective mass factor for vacuum

meff        = meff_vac;     # [-]       effective mass
melec       = meff*m0;      # [kg]      mass of an electron

eV2J        = 1.6e-19;      # [-]       energy conversion factor (eV to J)
J2eV        = 1/eV2J;       # [-]       energy conversion factor (J to eV)

# eigenvalues for ground state
n   = 1;
a1  = 10;
a2  = 100;
e1  = (hbar^2 * pi^2 *n^2) / (2 * melec * (a1*A)^2) * J2eV;
e2  = (hbar^2 * pi^2 *n^2) / (2 * melec * (a2*A)^2) * J2eV;

# show result
printf("KE(n=1) = %f eV\n", e1);
