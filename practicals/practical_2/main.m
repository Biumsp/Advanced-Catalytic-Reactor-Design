% Enrico Bussetti
% Advanced Catalytic Reactor Design 2020/2021
% Practical 2 

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

% =========================================================================
% Species and Reactions

% Species = N2 O2 o-Xyl PA H2O CO2

S = struct();

% Stoichiometric Coefficients
S.SC = [0 -3 -1 +1 +3 0;
        0 -10.5 -1 0 5 8;
        0 -7.5 0 -1 2 8]';

[S.Ns, S.Nr] = size(S.SC);
  
% Reaction enthalpies
S.Dhr = [-1285409 -4564000 -3278591];     % [kJ/kmol]

% Molecular Weights
S.MW = [28 32 106.16 148.12 18 44]';      % [kg/kmol]

% =========================================================================
% Reactor

R = struct();

% To Be Designed
R.d = 0.0254;           % Diameter [m] 
R.Across = pi*R.d^2/4;  % Coss Section Area [m^2]
R.L = 3;                % Length [m]
R.U = 0.096*3600;       % External Cooling Coefficient [kJ/m2/s/K]

% Fraction of tube with different catalyst density
frac_cat = 0.5;

% Dilution of the catalyst [-]
dilution = 0.5;

% =========================================================================
% Catalyst 

C = struct();

C.d   = 0.005;            % Particle diameter [m]
C.rho = 2100;             % Catalyst density [kg/m3]

% Void fraction of the reactor
DiameterRatio = R.d/C.d;
R.eps = 0.363+0.35*(exp(-0.39*DiameterRatio)); 

% =========================================================================
% Operating Conditions and Parameters

O = struct();

O.Tsalt = 335+273.15;         % Molten salts T [K]
O.Tin   = O.Tsalt;            % Feed temperature [K]
O.G     = 4900;               % Specific mass flow rate [kg/m2/h]

% Inlet molar fractions
oXylToAirRatio = 0.015;
n_in  = [0.79 0.21 oXylToAirRatio 0 0 0]'; 
O.x_in = n_in/sum(n_in);     % [-]

O.Cp  = 0.992;               % Specific Heat [kJ/kg/K]
O.mi  = 2.95e-5;             % Viscosity [Pa*s]
O.Pin = 1.3;                 % Inlet P [bar]

% Inlet mass fractions
O.omega_in = O.x_in.*S.MW/sum(O.x_in.*S.MW);  % [-]

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

% First part 
C.rho = C.rho*dilution;
IC = [O.omega_in', O.Tin, O.Pin];

[z,y] = ode15s(@PFR_pseudohomogeneous, 0:0.001:R.L*frac_cat, IC, [], ...
                S, R, C, O);

% Second part
C.rho = C.rho/dilution;

[z1,y1] = ode15s(@PFR_pseudohomogeneous, R.L*frac_cat:0.001:R.L, y(end,:), [], ...
                S, R, C, O);
            
z = [z', z1'];
y = [y', y1']';
            
omega = y(:, 1:S.Ns);
T = y(:, end-1);
P = y(:, end);

% -------------------------------------------------------------------------
% Results 
% -------------------------------------------------------------------------

NP=size(y(:,1)); %number of points
n = omega;

for j=1:NP
        n(j,:) = O.G*R.Across*omega(j,:)./S.MW';    
end
x = n./sum(n, 2);

% Conversion of ortho-xylene
Xoxyl = (n(1, 3) - n(:, 3))/n(1, 3);

% Selectivity of ortho-xylene into phtalic anhydride (PA)
sel_PA = n(:,4)./(n(1,3)-n(:,3) + 1d-12);
sel_PA(1) = sel_PA(2);

% Yield of PA
yield = sel_PA.*Xoxyl;

% Productivity of PA per tube
prod_per_tube = yield*n(1,3)*S.MW(4);      % [kg_PA/h/tube]

prod = prod_per_tube(end)*24;                   % [kg_PA/d/tube]

fprintf("Productivity per tube: %f [kg/d]\n", prod)
fprintf("Number of tubes: %f \n", 8000e3/(prod*365))

% -------------------------------------------------------------------------
% Graphical Post-Processing
% -------------------------------------------------------------------------

figure(1)
plot(z, Xoxyl, z, yield)
legend('Conversion', 'Yield')

figure(2)
yyaxis right
plot(z, T - 273.15)
yyaxis left
plot(z, sel_PA)
ylim([min(sel_PA), max(sel_PA)])
legend('T [°C]', 'Selectivity [-]')

figure(3)
plot(z, P)
legend('P [bar]')


    
    
