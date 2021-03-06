% Enrico Bussetti, 210304

% Advanced Catalytic Reactor Design
% Practical 1

% Solution of a dispersion PFR with first order kinetics using the 
% finite-difference differentiation scheme (with variable step-size)
% and a solver for linear systems based the Gauss-Sidler method with 
% a Successive-Over-Relaxation convergence accelerator

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
Np = length(y);

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

A = zeros(Np);
b = zeros(Np, 1);

A(1, 1) = 1 + 1/Pe/(y(2) - y(1));
A(1, 2) = - 1/Pe/(y(2) - y(1));

for i = 2:Np-1
    
    A(i, i-1) = 2/Pe/((y(i+1) - y(i))^2 + (y(i-1) - y(i))^2);
    A(i, i)   = - L/v*k - 4/Pe/((y(i+1) - y(i))^2 + (y(i-1) - y(i))^2) - 1/Pe/(y(i+1) - y(i));
    A(i, i+1) = 1/Pe/(y(i+1) - y(i)) + 2/Pe/((y(i+1) - y(i))^2 + (y(i-1) - y(i))^2);
    
end

A(Np, Np-1) = -1;
A(Np, Np)   =  1;

b(1) = Cin;

% First-Guess-Solution
FGS = linspace(Cin, 1, Np)';

% Find a solution using the Jacobi method with vector algebra
tic
C = gauss_seidler(A, b, FGS, 1e-5, 250e3, 1.9, 0);
toc

% Find a solution using the Jacobi method with for loops
tic
C = gauss_seidler(A, b, FGS, 1e-5, 250e3, 1.9, 0, 1);
toc

% Find solution with back-slash
tic
C_b = A\b;
toc

% -------------------------------------------------------------------------
% Graphical-Post-Processing
% -------------------------------------------------------------------------

figure
plot(y, C, y, C_b)
title('C vs y')
xlabel('y [-]')
ylabel('C [au ]')
