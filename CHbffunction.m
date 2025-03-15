clear all
clc
close all

tic;
% Read initial Chemical and Reaction Data from .csv files:

ChemData = readtable('CHchemdata.csv');
% Columns: 1:Element,2:Name,3:Density,4:MMass,5:Mvol,6:stdHForm,7:Cpn
% Units: Density:[kg/m^3],MMass:[kg/mol],Mvol:[m^3/mol],stdHForm:[J/mol],Cpn:[J/mol·K]
% Rows: 1:H2,2:O2,3:H20,4:C,5:CO,6:CO2,7:Fe2O3,8:Fe3O4,9:FeO,10:Fe,11:N2,12:CaCO3,13:CaO,14:SiO2,15:CaSiO3,16:Air

ReactionData = readtable('CHreactdata.csv');
% Columns: 1:Reaction,2:Description,3:f_factor,4:a_energy
% Units: f_factor:[s^-1],a_energy:[J]

% Set initial parameters:

tspan = 50; % Time span for simulation
CokeOreRatio = 1; % Coke : Ore molar ratios to investigate
HearthDiameter = 20; % Diameter of furnace [m]
HearthHeight = 10; % Height of furnace [m]
WallThickness = 2; % Thickness of BF walls [m]
ThermalConductivity = 30; % W/mK. Assuming Corundrum (Al2O3).
HotBlastRate = 60; % Volumetric flow rate of air [m^3/s]
HotBlastTemp = 1200; % Hot Blast Air Temperature [°C]
OxygenEnrichment = 0.5; % Assumes 20.95% of O2 in dry air. [%]
HydrogenInflow = 0.5; %%% Inflow of hydrogen at the bottom of the blast furnace
Temp_ext = 40; % Temperature outside Blast Furnace, assumed constant. [°C]
Temp_initial = 40; % Initial Blast Furnace temperature at start of model [°C]
SilicaImpurity = 2; % Silica content in Fe2O3[%]
nZones = 10; % Number of zones to model. Must be at least 3.
Res = 1; % Size of each time step for  modelling [s]. Lower is more accurate but requires more power.

% Pre-define arrays for efficiency:

[RateCoef,Rate,Hr,SpeciesVolume] = deal(zeros([tspan/Res height(ReactionData) nZones]));
[nTotal,CpnAv,VolUsed,VolUsedFraction,VolAvailable,VolAvailableFraction,HrNet] = deal(zeros([tspan/Res nZones]));
[dTdt,dhWall,dhProducts,dhBurdenDescentNet,dhGasAscentNet,dhHotBlast,HNet] = deal(zeros([tspan/Res nZones]));
[FeedRate,dnOut] = deal(zeros([tspan/Res 16]));
[dnBurdenDescent,dnGasAscent] = deal(zeros([tspan/Res 16 nZones]));
[VolUsedNew,VolAvailableBurden] = deal(zeros([(tspan/Res) 1]));
n = zeros([(tspan/Res)+1 16 nZones]);
dn = zeros([tspan/Res 16 nZones]);
Temp = zeros([(tspan/Res)+1 nZones]);

% Initialise Blast Furnace simulation:

for j = CokeOreRatio
    [RateCoef,Rate,Hr,SpeciesVolume,nTotal,CpnAv,VolUsed,VolUsedFraction,VolAvailable,VolAvailableFraction,HrNet,dTdt,dhWall,dhProducts,dhBurdenDescentNet,dhGasAscentNet,dhHotBlast,HNet,FeedRate,dnOut,dnBurdenDescent,dnGasAscent,VolUsedNew,VolAvailableBurden,n,dn,Temp,t] = CHblastfurnace(ChemData,ReactionData,tspan,j,HearthDiameter,HearthHeight,WallThickness,ThermalConductivity,HotBlastRate,HotBlastTemp,OxygenEnrichment,Temp_ext,Temp_initial,SilicaImpurity,nZones,Res,HydrogenInflow);
end

% Tata charge in terms of tonnes. Express feedrate in terms of that.

% For reduction rate, use grey colormap.

% Plot stuff

Results = zeros(length(0:Res:tspan),length(1:1:nZones));

ReductionRate = n(:,10,:)./(2*n(:,7,:)+3*n(:,8,:)+n(:,9,:)+n(:,10,:));

Results(:,:) = n(:,7,:);
%imagesc(0:Res:tspan,1:1:nZones,Results);hold on;colorbar;colormap(hot);xlabel('Time (s)'),ylabel('Zone');title('Blast Furnace zone Air reduction rate'),ylabel(colorbar, 'reduction rate');

ReductionRate(:,:) = n(:,10,:)./(2*n(:,7,:)+3*n(:,8,:)+n(:,9,:)+n(:,10,:));

%imagesc(0:Res:tspan,1:1:nZones,Temp');hold on;colorbar;colormap(hot);xlabel('Time (s)'),ylabel('Zone');title('Blast Furnace zone temperature vs. time'),ylabel(colorbar, 'Temperature (°C)');

% FINAL CALCULATIONS FOR PLOTTING

% Calculate Feed Rate in terms of mass:

FeedRateMass = (ChemData.Mmass(1:16).*FeedRate(:,:)')';

% Calculate percentage of reduction in each zone:

%PercentReduction(:,:) = (n(:,10,:)./(n(:,7,:)+n(:,8,:)+n(:,9,:)+n(:,10,:)))*100;

% Total percentage reduction is calculated in final column of PercentReduction:

%PercentReduction(:,nZones+1) = sum(PercentReduction,2)/nZones;

% Calculate dhTopGas = dhGasAscent from top zone = 1:

dhTopGas = dhGasAscentNet(:,1);

% Calculate average blast furnace temperature:

TempAvg = sum(Temp,2)/nZones;

% Cumulative Heat calculations:

HrNetCumulative = cumsum(HrNet);
dhWallCumulative = cumsum(dhWall);
dhProductsCumulative = cumsum(dhProducts(:,nZones));
dhTopGasCumulative = cumsum(dhTopGas(:,1));
dhHotBlastCumulative = cumsum(dhHotBlast(:,nZones));
HNetCumulative = cumsum(HNet);
FeoutCumulative = cumsum(dnOut(:,10));
SlagoutCumulative = cumsum(dnOut(:,15));
sprintf('Job done :-)')

plot(0:Res:tspan,FeoutCumulative);

toc;