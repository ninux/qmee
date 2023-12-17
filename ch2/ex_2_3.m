##############################################################################
##
## Copyright (C) 2023, Ervin Mazlagic
##
## This file is part of the personal solutions manual for "Quantum Mechanics
## for Electrical Engineers" by Dennis M. Sullivan (ISBN 978-0-470-87409-7),
## see "github.com/ninux/qmee".
##
##############################################################################

## Chapter 2 - Exercise 2-3
## Stationary States

# cleanup environment
clear all;
close all;

A       = 10e-10;
a_start = 0;
a_end   = 100;

x = [a_start:1:a_end].*A;

lambda = ((x(end)-x(1))*2)./[1:1:4];

f = 1./lambda;

psi = [...
    sin(2*pi*f(1)*x);
    sin(2*pi*f(2)*x);
    sin(2*pi*f(3)*x);
    sin(2*pi*f(4)*x)
];

pd = [...
    abs(psi(1,:)).^2;
    abs(psi(2,:)).^2;
    abs(psi(3,:)).^2;
    abs(psi(4,:)).^2
];

figure(1);
subplot(4,2,1);
    plot(x/A, psi(1,:), "-k");
    title("\\psi(x) for n=1");
    xlabel("Distance (A)");
    ylabel("\\psi(x)");
    grid on;
    ylim([-1, 1]);
subplot(4,2,3);
    plot(x/A, psi(2,:), "-k");
    title("\\psi(x) for n=2");
    xlabel("Distance (A)");
    ylabel("\\psi(x)");
    grid on;
    ylim([-1, 1]);
subplot(4,2,5);
    plot(x/A, psi(3,:), "-k");
    title("\\psi(x) for n=3");
    xlabel("Distance (A)");
    ylabel("\\psi(x)");
    grid on;
    ylim([-1, 1]);
subplot(4,2,7);
    plot(x/A, psi(4,:), "-k");
    title("\\psi(x) for n=4");
    xlabel("Distance (A)");
    ylabel("\\psi(x)");
    grid on;
    ylim([-1, 1]);

subplot(4,2,2);
    plot(x/A, pd(1,:), "-k");
    title("Probability Density |\\psi(x)|^2 for n=1");
    xlabel("Distance (A)");
    ylabel("|\\psi(x)|^2");
    grid on;
    ylim([-1, 1]);
subplot(4,2,4);
    plot(x/A, pd(2,:), "-k");
    title("Probability Density |\\psi(x)|^2 for n=2");
    xlabel("Distance (A)");
    ylabel("|\\psi(x)|^2");
    grid on;
    ylim([-1, 1]);
subplot(4,2,6);
    plot(x/A, pd(3,:), "-k");
    title("Probability Density |\\psi(x)|^2 for n=3");
    xlabel("Distance (A)");
    ylabel("|\\psi(x)|^2");
    grid on;
    ylim([-1, 1]);
subplot(4,2,8);
    plot(x/A, pd(4,:), "-k");
    title("Probability Density |\\psi(x)|^2 for n=3");
    xlabel("Distance (A)");
    ylabel("|\\psi(x)|^2");
    grid on;
    ylim([-1, 1]);

fig = gcf();
    set(fig, "PaperUnits", "centimeters");
    x_width = 21;
    y_width = 21;
    set(fig, "PaperPosition", [0 0 x_width y_width]);
    print(fig, "ex_2_3.png", "-dpng");
