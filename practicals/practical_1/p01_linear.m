% Enrico Bussetti, 210304

% Advanced Catalytic Reactor Design
% Practical 1

% Solution of a dispersion PFR with first order kinetics using the 
% finite-difference differentiation scheme (with variable step-size)
% and the Matlab built-in solver for linear systems

close all
clear variables

% -------------------------------------------------------------------------
% Data (arbitrary units)
% -------------------------------------------------------------------------

L  = 1;     % Length
v  = 1;     % Velocity
Pe = 10;    % Peclét material number (L*v/Di)

Cin = 10;   % Concentration of the feed (not at the inlet)
k   = 0.1;  % Kinetic constant 

% Vector of grid points
y  = 0:0.01:1;
% y = y.^2;
Np = length(y);

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

A = zeros(Np);
b = zeros(Np, 1);

A(1, 1) =   1/Pe/(y(2) - y(1)) + 1;
A(1, 2) = - 1/Pe/(y(2) - y(1));

for i = 2:Np-1
    
    A(i, i-1) =   2/Pe/((y(i+1) - y(i))^2 + (y(i-1) - y(i))^2);
    A(i, i)   = - 4/Pe/((y(i+1) - y(i))^2 + (y(i-1) - y(i))^2) - 1/Pe/(y(i+1) - y(i)) - L/v*k;
    A(i, i+1) =   2/Pe/((y(i+1) - y(i))^2 + (y(i-1) - y(i))^2) + 1/Pe/(y(i+1) - y(i));
    
end

A(Np, Np-1) = -1;
A(Np, Np)   =  1;

b(1) = Cin;

C = A\b;

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

% Analytical solution


% -------------------------------------------------------------------------
% Graphical-Post-Processing
% -------------------------------------------------------------------------

figure
plot(y, C, '-o')
title('C vs y')
xlabel('y [-]')
ylabel('C [au ]')
