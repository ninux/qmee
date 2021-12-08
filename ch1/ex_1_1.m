% exercise 1-1
printf("Exericise 1-1\n");

% cleanup environment
clear all;
close all;

% definitions
phi_Ti  = 4.33;       % [eV]    work function of Titanium
c0      = 3e8;        % [m/s]   speed of light in vacuum
h       = 6.625e-34;  % [Js]    Planck's constant

ev2j    = 1.6e-19;    % [Joule]

% calculations
E = phi_Ti*ev2j;      % [J]     required energy
lambda  = (c0*h)/(E); % [m]     wavelength

% show results
printf("The required wavelength is %3.1f nm\n", lambda*1e9);