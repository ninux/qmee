##############################################################################
##
## Copyright (C) 2023, Ervin Mazlagic
##
## This file is part of the personal solutions manual for "Quantum Mechanics
## for Electrical Engineers" by Dennis M. Sullivan (ISBN 978-0-470-87409-7),
## see "github.com/ninux/qmee".
##
##############################################################################

## Chapter 1 - Exercise 1-3

# cleanup environment
clear all;
close all;

# definitions
me        = 9.11e-31;     # [kg]    electron mass (in vacuum)
h         = 6.625e-34;    # [Js]    Planck's constant
lambda_1  = 10e-9;        # [m]     initial wavelength
V         = 0.02;         # [eV]    potential

eV2J      = 1.6e-19;      # [Joule] electron-Volt to Joule factor

# calculations
KE1       = (1/2/me)*(h/lambda_1)^2;
KE2       = KE1 + V*eV2J;
lambda_2  = h/sqrt(2*KE2*me);

# show results
printf("The wavelength is #3.1f nm\n", lambda_2*1e9);
