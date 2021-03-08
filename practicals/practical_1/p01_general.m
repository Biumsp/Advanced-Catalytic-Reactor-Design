% Enrico Bussetti, 210304

% Advanced Catalytic Reactor Design
% Practical 1

% Solution of a dispersion PFR with linear kinetic mechanisms and no 
% change of moles. 
% Finite-Difference differentiation scheme (with variable step-size).
% CDS for internal points.
% Performance comparison between different linear systems solvers.

close all
clear variables

% -------------------------------------------------------------------------
% Data [SI units]
% -------------------------------------------------------------------------

% =========================================================================
% Geometry and flows

d  = 2.54e-2;           % Diameter [m]
V  = 1.5e-3;            % Volume [m^3]
L  = V*4/pi/d^2;        % Length [m]

Q = 0.01;               % Volumetric flow rate [m^3/s]
v  = 1;                 % Velocity [m/s]

% =========================================================================
% Properties

 % Dispersion coefficients [m^2/s]
Di = [1e-5, 1e-5, 1e-5]'; 
Pe = L*v./Di;              % Peclét material number [-]

% Concentration of the feed (not at the inlet) [mol/m^3]
Cin = [10, 0, 0];          

% =========================================================================
% Mechanism 

% Stoichiometric Coefficients matrix 
SC = [-1,  1, 0;
       0, -1, 1]';       
   
% Reaction Orders matrix
RO = [1, 0, 0;
      0, 1, 0]';         
  
% Kinetic constants  
k  = [0.1, 0.1];        

% =========================================================================
% Discretization

y = 0:0.01:1;
y = y.^2;

d = [0 , y(2:end)-y(1:end-1)];
g = @(n, m) d(n)^2 + d(n)*d(m);

Np = length(y);
[Ns, Nr] = size(SC);

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

% =========================================================================
% Definition of the system

A = zeros(Np*Ns);
b = zeros(Np*Ns, 1);

r = zeros(1, Nr);

for kk = 1:Np
    for uu = 1:Ns
        
        ii = (kk - 1)*Ns + uu;
        
        % Boundary points
        if kk == 1
            
            A(ii, ii)    =   1/Pe(uu)/d(kk+1) + 1;
            A(ii, ii+Ns) = - 1/Pe(uu)/d(kk+1);
            
            b(ii) = Cin(uu);
 
            continue
        elseif kk == Np
            
            A(ii, ii-Ns) =   1;
            A(ii, ii)   =  -1;
            
            break
        end
        disp(kk)
        % Central points
        A(ii, ii-Ns) =  d(kk+1)/g(kk, kk+1) + 2/Pe(uu)/g(kk, kk+1);
        A(ii, ii)    =  1/d(kk+1) - g(kk+1, kk)/g(kk, kk+1)/d(kk+1)...
                        - 2/Pe(uu)/g(kk, kk+1)*(1 + d(kk)/d(kk+1));
        A(ii, ii+Ns) =  2*d(kk)/Pe(uu)/d(kk+1) - 1/d(kk+1) + d(kk)/g(kk, kk+1);
        
        for mm = 1:Ns
            A(ii, (kk - 1)*Ns + mm) =  L/v*SC(uu, :)*(k'.*RO(mm, :)');
        end
        
    end
end

% =========================================================================
% Backslash

t_in = cputime;
C = A\b;
t_out = cputime;

fprintf('Method: Backslash, time = %f\n', t_out-t_in);

% =========================================================================
% Gauss-Seidler (vector algebra)

FGS = zeros(Np*Ns);
for ii = 1:Np
    for jj = 1:Ns
        kk = (ii-1)*Ns + jj;
        FGS(kk) = Cin(jj);
    end
end


t_in = cputime;
C = gauss_seidler(A, b, FGS, 1e-5, 250e3, 1.8, 0, 1);
t_out = cputime;

fprintf('Method: Gauss-Seidler (vectors), time = %f\n', t_out-t_in);

% =========================================================================
% Gauss-Seidler (for loops)

t_in = cputime;
C = gauss_seidler(A, b, FGS, 1e-5, 250e3, 1.8, 0, 1);
t_out = cputime;

fprintf('Method: Gauss-Seidler (for loops), time = %f\n', t_out-t_in);

% =========================================================================
% Biconjugate-Gradient-Stabilized

t_in = cputime;
C = bicgstabl(A, b);
t_out = cputime;

fprintf('Method: Bicgstabl, time = %f\n', t_out-t_in);

% -------------------------------------------------------------------------
% Graphical-Post-Processing
% -------------------------------------------------------------------------

for ii = 1:Np
    for jj = 1:Ns
        
        kk = (ii -1)*Ns + jj;
        CC(ii, jj) = C(kk);
        
    end
end

figure
plot(y, CC, '-o')
title('C vs y')
xlabel('y [-]')
ylabel('C [au ]')
