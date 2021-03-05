% -y'+(1/Pe)*y'' + nu(i)*k*y(A) = 0
% y0 = y(0) - (1/Pe9 + y'(0)
% y'(0) = 0

% A, B

% r = k*cA
% -------------------------------------------------------
clc
clearvars -except Np

NEQ = 2;

% Reactor 
d = 0.0254; %m
V = 0.0015; % m3
L = V/(pi*d^2/4); %m

Q = 0.01; % m3/s

v = Q/(pi*d^2/4); %m/s

D(1) = 1e-5; % m2/s
D(2) = 1e-5; % m2/s

k = 20;

% A => B
nu = [-1, 1];

Pe = L*v./D;

% Discretization of the domain
delta_z_star = 1/(Np - 1);
z_star = 0:delta_z_star:1;

% ------------------------------------
y0(1) =1; 
y0(2) =0;

% ------------------------------------

for j = 1:NEQ
    for i = 1:Np
        
        y_guess((i - 1)*NEQ + j) = y0(j);
        
    end
end



y = fsolve('pfr_with_dispersion', y_guess, [],...
           L, v, Pe, k, Np, NEQ, y0, delta_z_star, nu);


for j = 1:NEQ
    for i = 1:Np
        if j == 1
            y_a(i) = y((i-1)*NEQ + j);
        else
            y_b(i) = y((i-1)*NEQ + j);
        end
    end
end

plot(z_star, y_a,  z_star, y_b)

hold on 
if Np < 200
    Np = Np*2;
    p01_inclass
end