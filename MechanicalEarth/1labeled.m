clear all % clears all variables
close all % closes all figures
clc

% Vary magma viscosity
p   = linspace(1,3,100);  % sets of an array to vary magma viscosity expoentially, between 1 and 3 with incriements of .02
eta = 10.^p;   % creates an array varying viscosities of powers of p, that are linearly incremented by .02 in log_10 space
  
  
% Vary dike width
  w   = linspace(1,100,100); % creates an array of varying dike widths from 1 to 100 with step size 1
  [ETA,W] = meshgrid(eta,w); % creates a 2d array of the combinations of width and viscosity, where viscosity is the first entry and width is the second for each pair.
  
  
% Make plots
Figure
  contourf(log10(ETA),W,log10(Re)) % plots the assosiated reynolds number using contours with each pair (viscosity, width)
xlabel('log viscosity (Pa s)')   % labels the plot

