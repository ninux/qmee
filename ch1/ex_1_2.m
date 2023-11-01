##############################################################################
##
## Copyright (C) 2023, Ervin Mazlagic
##
## This file is part of the personal solutions manual for "Quantum Mechanics
## for Electrical Engineers" by Dennis M. Sullivan (ISBN 978-0-470-87409-7),
## see "github.com/ninux/qmee".
##
##############################################################################

## Chapter 1 - Exercise 1-2

# cleanup environment
clear all;
close all;

# definitions
phi_Ti  = 4.33;         # [eV]    work function of Titanium
c0      = 3e8;          # [m/s]   speed of light in vacuum
h       = 6.625e-34;    # [Js]    Planck's constant

ev2j    = 1.6e-19;      # [Joule]

# calculations
E       = phi_Ti*ev2j;  # [J]     required energy
lambda  = (c0*h)/(E);   # [m]     wavelength

# show results
printf("The required wavelength is %3.1f nm\n", lambda*1e9);
