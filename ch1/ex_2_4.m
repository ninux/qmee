##############################################################################
##
## Copyright (C) 2023, Ervin Mazlagic
##
## This file is part of the personal solutions manual for "Quantum Mechanics
## for Electrical Engineers" by Dennis M. Sullivan (ISBN 978-0-470-87409-7),
## see "github.com/ninux/qmee".
##
##############################################################################

## Chapter 1 - Exercise 2-4

# cleanup environment
clear all;
close all;

# definitions
dx = 5e-9;  # [m] interval

k3 = 1;     # relative magnitude in interval 3
k4 = 1/2;   # relative magnitude in interval 4
k5 = 1/10;  # relative magnitude in interval 5

# calculations
q = 1/(dx*(k3^2 + 2*(k4^2) + 2*(k5^2)));  # probability density factor
                                          # q = A^2

p3 = k3^2*q*dx; # probability for interval 3
p4 = k4^2*q*dx; # probability for interval 4 and 2
p5 = k5^2*q*dx; # probability for interval 5 and 1

p = p3 + 2*p4 + 2*p5;;

# show results
printf("Inteval probabilities:  \
    \n\tp3 = %2.3f\
    \n\tp4 = %2.3f\
    \n\tp5 = %2.3f\n",
    p3, p4, p5);
printf("Total probability:      \
    \n\tp  = %2.3f\n", p);
