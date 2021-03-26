function yy = PFR_pseudohomogeneous(z, var, S, R, C, O)

    x = var(1:6)./S.MW/sum(var(1:6)./S.MW);
    T = var(end-1);
    P = var(end);
    
    yy = zeros(length(x) + 2, 1);
    
    Ctot = P*1e5/8.314/T;               % [mol/m^3]
    
    MW_mix = sum(x.*S.MW);
    rho = Ctot*MW_mix/1000;             % [kg/m^3]
    
    v0 = O.G/rho/3600;                  % [m/s]
    
    %===============================
    % Kinetics
    
    K(1) = exp(19.837- 13636/T);
    K(2) = exp(18.970- 14394/T);
    K(3) = exp(20.860- 15803/T);
    
    % The exp reaction rate needs P in bar    
    RR(1) = K(1)*P*x(3)*P*x(2);
    RR(2) = K(2)*P*x(3)*P*x(2);
    RR(3) = K(3)*P*x(4)*P*x(2);
    
    r = zeros(6,1);
    
    for ii = 1:6
        r(ii) = r(ii) + sum(S.SC(ii,:).*RR);
    end
    
    yy(1:S.Ns) = r*(1-R.eps)*C.rho.*S.MW/O.G;
    
    yy(end-1) = (-sum(S.Dhr.*RR)*(1-R.eps)*C.rho + ...
                R.U*(4/R.d)*(O.Tsalt - T))/O.G/O.Cp;
    
    yy(end) = -(150*(1-R.eps)^2/(R.eps^3)*O.mi*v0/C.d^2 + ...
              + 7/4*rho*v0^2/C.d*(1-R.eps)/(R.eps^3))/1e5;
    a = 2;
end