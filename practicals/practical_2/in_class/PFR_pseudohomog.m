function PFR_pseudohomog=PFR_pseudohomog(t,y)
%
%--ACRD (Prof. Matteo Maestri) - 2020/2021 - March 19, 2021  
%
global Stoichiometry deltaH TubeDiameter ParticleDiameter ReactorLength...
    CatalystDensity ExternalCoolingCoefficient VoidFraction MoltenSaltsTemperature...
    FeedTemperature G SpecificHeatP Viscosity FeedPressure mw

%molar_fraction:

for i=1:6
    mol(i)=y(i)/mw(i); %kmol (assumed: bc=>1 kg)
end
MolarFraction=mol./(sum(mol));

Temperature=y(end-1);
Pressure_bar=y(end); %bar
Pressure_Pa=Pressure_bar*1e5; %Pa

DensityMolarGas=1e-3*Pressure_Pa/8.314/Temperature; %kmol/m3

mw_av=0;
for i=1:6
    mw_av=mw_av+MolarFraction(i)*mw(i);
end

DensityMassGas=DensityMolarGas*mw_av; %kg/m3

SuperficialVelocity=G/DensityMassGas/3600; %m/s
%--
%KINETICS:
ReactionK(1)=exp(19.837-13636/Temperature);
ReactionK(2)=exp(18.970-14394/Temperature);
ReactionK(3)=exp(20.860-15803/Temperature);
ReactionRate(1)=ReactionK(1)*Pressure_bar*MolarFraction(3)*...
    Pressure_bar*MolarFraction(2);
ReactionRate(2)=ReactionK(2)*Pressure_bar*MolarFraction(3)*...
    Pressure_bar*MolarFraction(2);
ReactionRate(3)=ReactionK(3)*Pressure_bar*MolarFraction(4)*...
    Pressure_bar*MolarFraction(2);
%--------
for i=1:6
    NetRateProduction(i)=0;
end

for j=1:3
    for i=1:6
        NetRateProduction(i)=NetRateProduction(i)+...
            Stoichiometry(j,i)*ReactionRate(j);
    end
end

for i=1:6
    PFR_pseudohomog(i)=NetRateProduction(i)*(1-VoidFraction)...
        *CatalystDensity*mw(i)/G;
end

%--Energy equation:
RateH=0;
for j=1:3
    RateH=RateH+deltaH(j)*ReactionRate(j); %kJ/kg_cat/h
end

PFR_pseudohomog(7)=(-RateH*(1-VoidFraction)*CatalystDensity+...
    ExternalCoolingCoefficient*(4/TubeDiameter)*...
    (MoltenSaltsTemperature-Temperature))/G/SpecificHeatP;

%--Ergun Equation:
PFR_pseudohomog(8)=-(150*(1-VoidFraction)^2/(VoidFraction^3)...
    *Viscosity*SuperficialVelocity/ParticleDiameter^2+...
    +7/4*DensityMassGas*SuperficialVelocity^2/ParticleDiameter*...
    (1-VoidFraction)/(VoidFraction^3))/1e5;

PFR_pseudohomog=PFR_pseudohomog';
main
a = 2;

end


