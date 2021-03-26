clear all
%
%--ACRD (Prof. Matteo Maestri) - 2020/2021 - March 19, 2021  
%
%=> N2 O2 o-xyl PA H2O CO2 Temperature Pressure
%
% example: y(2)=omega(O2)
%

global Stoichiometry deltaH TubeDiameter ParticleDiameter ReactorLength...
    CatalystDensity ExternalCoolingCoefficient VoidFraction MoltenSaltsTemperature...
    FeedTemperature G SpecificHeatP Viscosity FeedPressure mw

%Reactions:

Stoichiometry=[0 -3 -1 +1 +3 0
               0 -10.5 -1 0 5 8
               0 -7.5 0 -1 2 8];
           
deltaH=[-1285409 -4564000 -3278591]; %kJ/mol

%Species:
mw=[28 32 106.16 148.12 18 44]; % kg/kmol

%-REACTOR PARAMETERS:
TubeDiameter=0.0254; %m - TO BE DESIGNED
ParticleDiameter=0.005; %m 
ReactorLength=3; %m - TO BE DESIGNED 
CatalystDensity=2100; %kg/m3
ExternalCoolingCoefficient=0.096*3600; %kJ/m2/h/K
DiameterRatio=TubeDiameter/ParticleDiameter;
VoidFraction=0.363+0.35*(exp(-0.39*DiameterRatio)); 

%-OPERATING CONDITIONS:
MoltenSaltsTemperature=335+273.15; %K
FeedTemperature=MoltenSaltsTemperature; 
G=4900; %kg/m2/h
oXylToAirRatio=0.013;
n_in=[0.79 0.21 oXylToAirRatio 0 0 0]; 
MolarFractionIN=n_in./sum(n_in);
SpecificHeatP=0.992; %kJ/kg/K
Viscosity=2.95e-5; %Pa*s
FeedPressure=1.3; %bar

%==
basis_calc=1; %kmol YOUR CHOICE!
for i=1:6
    mass(i)=basis_calc*MolarFractionIN(i)*mw(i); %kg
end
MassFractionIN=mass./sum(mass);

%--SYSTEM INTEGRATION:
[t,y]=ode15s('PFR_pseudohomog',...
              [0:0.001:ReactorLength],...
              [MassFractionIN FeedTemperature FeedPressure]);
          
%POST PROCESSING:

NP=size(y(:,1)); %number of points

Across=pi*TubeDiameter^2/4;

for j=1:NP
    for i=1:6
        n(j,i)=G*Across*y(j,i)/mw(i); %kmol/h
    end
    MolarFraction(j,:)=n(j,:)./sum(n(j,:));
end

for j=1:NP
    conv_oxyl(j)=(n(1,3)-n(j,3))/n(1,3);
end

for j=1:NP
    sel_PA(j)=n(j,4)/(n(1,3)-n(j,3)+1d-12);
end

for j=1:NP
    yield(j)=sel_PA(j)*conv_oxyl(j);
end

for j=1:NP
    prod_per_tube(j)=yield(j)*n(1,3)*mw(4); %kg_PA/h/tube
end




    
    
