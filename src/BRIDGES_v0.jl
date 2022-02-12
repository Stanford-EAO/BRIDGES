################################################################################
# Building Resilient Integrated, Decarbonized Gas-Electric Systems (BRIDGES)
# Version: 0.0
# Author: Gregory Von Wald
# Objective: To permit integrated design and operations optimization for 
# system planning of integrated public-interest (i.e., gas and electric) energy systems
# across time horizons with declining constraints on greenhouse gas emissions.
# Model minimizes the present value of total system costs (investment plus operations)
# across a multi-year investment time horizon, with co-optimized system operations 
# simulated for a set of weighted, representative periods (e.g., days) 
# composed of operational time steps (e.g., hours)
#
################################################################################
import Pkg
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Clustering")
# Pkg.add("JuMP")
# Pkg.add("Distances")
# Pkg.add("Gurobi")
# Pkg.add("Tables")
#Pkg.add("DelimitedFiles")
using DataFrames, CSV, Tables, Clustering, JuMP, Distances, Gurobi, Dates

# creates output folder automatically
function mk_dirs()
    timestamp = Dates.format(now(), "YYYYmmdd-HHMMSS")
    top_level_path = abspath(joinpath(@__DIR__, ".."))
    
    data_dir *= top_level_path * "data/"
    out_dir = joinpath(top_level_path, "output", "$timestamp")
    
    @assert !ispath(out_dir) "File name already taken!"
    mkpath(out_dir)
    return out_dir * "/", data_dir *
end

out_dir, data_dir *= mk_dirs()

# Toggles to turn on/off different model relaxations and functionality

# Indicate whether you want to fix gas flow in the direction specified in the GasTransmission and ElecTransmission csv files (=0) or introduce binary variables to allow bi-directional flow directions (=1)
GASFLOW_DIRECTIONS = 0

# Indiciate whether gas distribution retirement should be contemplated (with the associated savings)
gasdistretirement_allowed = 0

# Indicate whether transmission expansion/retirement should be explicitly considered
TRANSMISSION_EXPANSION = 0
# Indicate whether to include constraints that link representative time periods
# for tracking storage state of charge (if = 0, periodicity constraints are imposed for each rep. period)
LINKED_PERIODS_STORAGE = 1
# for generator operations such as min up/down times and ramp rates (if = 0, constraints only apply within each rep. period)
LINKED_PERIODS_GENOPS = 1

# Indicate whether steady-state physics should be simulated for the electric and gas systems
# If = 0, the flows of power and gas will be governed by simple transport models
STEADYSTATE_ELEC = 1
STEADYSTATE_GAS = 1

appliance_decisions = 1
hybrids_allowed = 0
liquids_allowed = 0


# Clustering parameters
T_inv = 1               # Number of investment time periods modeled
N_Periods = 1         # Number of representative operational time slices modeled for each investment period
HOURS_PER_PERIOD = 24   # Number of hourly time steps in each rep. op. time slice
# Clustering technique to use for generating representative days
# Options include:
# (a) "average",
# (b) "ward",
# (c) "kmeans"
clustering_case = "average"     

BaseYear = 2019                             # Initial year
Years = [2020,2025,2030,2035,2040]          # Modeled investment years
AllYears = [2019,2020,2025,2030,2035,2040]  # List of all calendar years

# For file import, and to permit many sensitivity scenarios,
# each configuration can be specified by a system and a region
# which corresponds to the necessary import files
system = "Network"
region = "Cold"
num = ""
bounding_steady_states = 0

# Base case is mid, low, low
biomethane = "Mid"
industrials = "Low"
buildingretrofits = "Low"

CleanElecCosts = "Mid"
CleanGasCosts = "Mid"

retirements_case = "NoGasDistRet"
gasdistretirement_forced = [0,0,0,0,0]

offsets_case = "NoOffsets"
# No = 0
# Low = 0.1
# Mid = 0.25
# High = 0.5
# Unlimited = 1.0
maxOffsets = 0.0*ones(T_inv)                  # % of gross emissions
NETSCost = "Mid"

offsets_Cost = [650, 550, 450, 350, 250]                        # $/tCO2e
if NETSCost == "High"
    global offsets_Cost = [650, 600, 550, 500, 450]                        # $/tCO2e
end

case = "WithAppDecisions"
if appliance_decisions == 0
    global case = "NoAppDecisions"
end

GasQuality = "Nodal" # "Annual", "No"

EITrajectory = "MidEI"
buildrate = "MidBuild"

loadgrowthscen = "NoGrowth"

br = 1.0                        # build rate multipliuer
transmission_multiplier = 1.0   # electric transmission rating multiplier

max_biomethane_share = 0.30     # annual system-wide limitation on biomethane production (as a share of initial core gas demands)
if biomethane == "High"
    global max_biomethane_share = 0.50
end
if biomethane == "Low"
    global max_biomethane_share = 0.1
end

# Set to 10% of baseline electricity demands
baselinegasdemand_multiplier = 1
if industrials == "No"
    global baselinegasdemand_multiplier = 0.001
end
if industrials == "Mid"
    global baselinegasdemand_multiplier = 2.5
end
if industrials =="High"
    global baselinegasdemand_multiplier = 5
end


WACC = 0.07                         # Weighted average cost of capital (WACC) applied to annualize utility-scale capital investments
WACC_APPLIANCES = 0.15              # Weighted average cost of capital (WACC) applied to annualize customer-scale appliance investments
societal_discounting = 0.01         # discount rate [%] applied to discount future costs to present value
LoadGrowthRate = 0.0                # [%/year] of baseline electricity/gas demand growth 


# Specify emissions intensity targets for the electricity sector and gas sector
# with Slow and Fast sensitivity scenarios possible.
# For reference, a natural gas-fired generator will yield ~500kg/MWh elec.
# coal-fired generators will yield ~1000kg/MWh elec.
# fossil natural gas delivered for direct-use will release ~181kg/MWh thermal
EI_ElecSector = [500,250,75,50,0]  # kg/MWh electricity generated
EI_GasSector = [200,150,50,15,0]   # kg/MWh gas delivered (to core customers)
if EITrajectory == "SlowEI"
    global EI_ElecSector = [500,500,250,250,0]  # kg/MWh electricity generated
    global EI_GasSector = [200,200,90,90,0.0]   # kg/MWh gas delivered (to core customers)
end
if EITrajectory == "FastEI"
    global EI_ElecSector = [1000,500,100,10,0]  # kg/MWh electricity generated
    global EI_GasSector = [200,100,20,2,0.0]    # kg/MWh gas delivered (to core customers)
end


# Print outs just to confirm the proper scenario is running
println("$(system) $(num)")
println(region)
println(case)
println(EITrajectory)
println(buildrate)
println("$(GasQuality) GasQuality")
println("$(biomethane) Biomethane")
println("$(industrials) Industrials")
println("$(buildingretrofits) buildingretrofits")
println(offsets_case)
println(retirements_case)
println("$(CleanElecCosts) Cost Clean Elec")
println("$(CleanGasCosts) Cost Clean Gas")


## Compute the discounting factor for each investment period's annualized costs based on the number of years represented by each period
# Eq. 2.72 in Von Wald thesis.
EndOfCostHorizon = 2050             # Specifies the horizon over which societal costs should be included in objective function
discountfactor = zeros(T_inv)
for i = 1:T_inv
    if i < T_inv
        for j = 1:Int(Years[i+1]-Years[i])
            discountfactor[i] = discountfactor[i] + 1/((1+societal_discounting)^(Years[i]-BaseYear+j-1))
        end
    end
    if i == T_inv
        for j = 1:Int(EndOfCostHorizon - Years[i] - 1)
            discountfactor[i] = discountfactor[i] + 1/((1+societal_discounting)^(Years[i]-BaseYear+j-1))
        end
    end
end

# Specify the number of residential and commercial customers on each distribution system
N_ResCust = 50000
N_CommCust = 2000
# Residential system costs are estimated at $350/customer-year
# Commercial system costs are estimated at $1200/cust.-year
Costs_GasDistSys = 350*N_ResCust + 1200*N_CommCust # Dollars per distsyst. per year

# Cost of electric distribution infrastructure is based on (Fares, et. al, )
# as a function of the peak electrical demand
Cost_DistributionInfrastructure = 73        # $/kW peak

# Cost values for transmission expansion/retirement modeling
# Not currently employed in v.0.0
ElecTransmissionCapitalCosts = 0 # $/MW
ElecTransmissionOperatingCosts = 0 # $/MW
GasTransmissionCapitalCosts = 0 # $/km
GasTransmissionOperatingCosts = 0 # $/km


T_ops = N_Periods                                           # Number of operational periods simulated for each investment year
t_ops = HOURS_PER_PERIOD                                    # Number of hours simulated for each operational period
HOURS_PER_YEAR = 8760   # hours/year
Periods_Per_Year = Int(HOURS_PER_YEAR/HOURS_PER_PERIOD)     # Number of rep. operational periods per calendar year


## Conversion constants
SEC_PER_HOUR = 3600     # sec/hour
MJ_PER_MWh = 3600       # MJ/MWh
MWh_PER_MMBTU = 0.293   # MWh/MMBtu
EF_NG = 53.1/MWh_PER_MMBTU/1000    #tCO2/MWh NG     # Per EPA emissions inventory (only CO2, no CH4 leakage)
HHV_H2 = 12.7  # MJ/standard m3
HHV_CH4 = 37.7  # MJ/standard m3
LHV_H2 = 10.24 # MJ/standard m3
LHV_CH4 = 33.9  # MJ/standard m3

## Gas/Power flow parameters
SLACK_BUS = 1
BASEMVA = 100
SLACK_NODE = 1
PRESSURE_MIN = (3447380/10^6)^2 # 500 psi squared
PRESSURE_MAX = (5515808/10^6)^2 # 800 psi squared
SLACK_NODE_Pressure = PRESSURE_MAX

################################################################################
# Baseline Energy Demands
################################################################################
# Import the baseline electrical demands [MWh/hr] across all system nodes
D_Elec2 = CSV.read(data_dir * "BaselineElectricDemands$(system)$(region).csv",DataFrame)

NODES_ELEC = length(D_Elec2[1,:])       # Number of electrical nodes is specified based on the number of columns in D_Elec2
# Set up a new array to hold electrical demand info.
D_Elec = zeros(8760,NODES_ELEC)         
for n = 1:NODES_ELEC
    D_Elec[:,n] = D_Elec2[:,n]
end

# Import the baseline gas demands [MWh/hr] across all system nodes
D_Gas2 = CSV.read(data_dir * "BaselineGasDemands$(system)$(region).csv",DataFrame)

NODES_GAS = length(D_Gas2[1,:])         # Number of gas nodes is specified based on the number of columns in D_Gas2
# Set up a new array to hold gas demand info.
D_Gas = zeros(8760,NODES_GAS)
for n = 1:NODES_GAS
    D_Gas[:,n] = baselinegasdemand_multiplier*D_Gas2[:,n]
end

# Names/Numbers of nodes to facilitate look-ups
REGIONS_ELEC = String.(names(D_Elec2))
REGIONS_GAS = String.(names(D_Gas2))

# Maximum fossil gas supply at the boundary/slack node
MAXSLACK = zeros(1,NODES_GAS)   # Set slack supply to zero everywhere except for the boundary node
MAXSLACK[SLACK_NODE] = 10000    # [MW] (Arbitrarily large, but not so large as to trigger numerical issues in optimization)


################################################################################
# End-use appliance demands
################################################################################
EndUseAppliances = CSV.read(data_dir * "EndUseAppliances$(system)$(num).csv",DataFrame)
APPLIANCES = length(EndUseAppliances[:, :1])        # Number of appliance classes modeled
ApplianceServices = EndUseAppliances[:,5]           # End-use services satisfied by each appliance class
PrimeMover_APPLIANCES = EndUseAppliances[:,6]       # Technology type for each appliance class
InitialAppliancePopulation = EndUseAppliances[:,7]  # Initial appliance population [no. units]
ApplianceLifetime = EndUseAppliances[:,8]           # Expected appliance lifetime [years]
IS_HYBRID = EndUseAppliances[:,9]                   # Indicator for whether the appliance is hybrid gas-electric
upgrade_cost = EndUseAppliances[:,10]               # Building infrastructure upgrade costs associated with transitioning to this appliance [$]
if buildingretrofits == "High"
    global upgrade_cost = EndUseAppliances[:,11]
end
CRF_APPLIANCES = (WACC_APPLIANCES.*(1+WACC_APPLIANCES).^ApplianceLifetime)./((1+WACC_APPLIANCES).^ApplianceLifetime .- 1)   # Capital recovery factor [yr^-1] for annualizing appliance investments

## Create a matrix that maps each appliance to the energy service that it satisfies
################################################################################
SERVICES = length(unique(EndUseAppliances[:,5]))
ServiceList = unique(EndUseAppliances[:,5]) # List of all energy services (i.e., residential space heating, residential water heating, commercial space heating, commercial water heating, etc.)
AppliancesToServices = zeros(APPLIANCES,SERVICES)
# For each appliance, a, put a 1 in the column corresponding to its energy service s
for a = 1:APPLIANCES
    AppliancesToServices[a,findfirst(occursin.([(ApplianceServices[a])],ServiceList))] = 1
end

# Pre-compute the cumulative failure fraction for each appliance in each investment period
# See Eq. 2.9 in Von Wald thesis
################################################################################
cumulativefailurefrac = zeros(APPLIANCES,T_inv,T_inv)
failureProb = zeros(APPLIANCES,150)
failureArchive = CSV.read(data_dir ** * "failureProb.csv",DataFrame)
# First, calculate failure probabilities for each appliance class in each year of its lifetime from 1 to 50.
for a = 1:APPLIANCES
    for i = 1:50
#       Using Poisson probability distribution to assess failure fractions
#       failureProb[a,i] = exp(-ApplianceLifetime[a])*(ApplianceLifetime[a]^(i))/factorial(i) 
#       However, the factorial function in Julia won't go over 20! which limits our ability to model long-lived equipment
#       Instead, we use an exogenous file generated using python's factorial function.
        failureProb[a,i] = failureArchive[Int(ApplianceLifetime[a]),i]
    end
    # Ensures that the sum across each row equals 1 (i.e., no appliance lasts longer than 50 years)
    failureProb[a,50] = 1 - sum(failureProb[a,1:49])
end
# Second, compute the cumulative failure fraction for each appliance type, in each investment year
# i.e., cumulativefailurefrac[a,v,t] corresponds to the cumulative failure fraction of appliances of type a
# that were installed in investment period v, that will fail by investment period t.
for a = 1:APPLIANCES
    for v = 1:T_inv
        for t = 1:T_inv
            cumulativefailurefrac[a,v,t] = round(sum(failureProb[a,1:max(Years[t]-Years[v],1)]),digits = 4)  # rounding to avoid numerical issues in the optimization program due to small coefficients
            if t == v
                cumulativefailurefrac[a,v,t] = 0.0
            end
        end
    end
end

## Appliance level energy demand profiles (hourly)
# In MWh/hr per unit, for each hour in a typical year; then must be clustered down
################################################################################
ApplianceProfilesGAS2 = CSV.read(data_dir * "ApplianceProfiles_GAS$(system)$(region).csv",DataFrame)
ApplianceProfilesELEC2 = CSV.read(data_dir * "ApplianceProfiles_ELEC$(system)$(region).csv",DataFrame)
#ApplianceProfilesLIQ2 = CSV.read("ApplianceProfiles_LPG$(system)$(region).csv",DataFrame)

#ApplianceProfilesLIQ = zeros(8760,length(ApplianceProfilesLIQ2[1,:]))
ApplianceProfilesGAS = zeros(8760,length(ApplianceProfilesGAS2[1,:]))
ApplianceProfilesELEC = zeros(8760,length(ApplianceProfilesELEC2[1,:]))
# Rounded to avoid introducing numerical issues
for i = 1:length(ApplianceProfilesGAS2[1,:])
    ApplianceProfilesGAS[:,i] = round.(ApplianceProfilesGAS2[:,i], digits = 8)
    ApplianceProfilesELEC[:,i] = round.(ApplianceProfilesELEC2[:,i], digits = 8)
#    ApplianceProfilesLIQ[:,i] = round.(ApplianceProfilesLIQ2[:,i], digits = 8)
end

## Growth rates used for forecasting and back-casting appliance sales
# Set all growth rates to zero
################################################################################
ServicesGrowthRate = zeros(SERVICES,1)      # %\year
HistoricalGrowthRate = zeros(APPLIANCES,1)  # %\year  
ForecastGrowthRate = zeros(APPLIANCES,1)    # %\year    

# Distribution systems are set up to potentially exist at the sub-transmission nodal level
# i.e., multiple distribution systems may exist and operate independently at the
# same transmission node.
################################################################################
DISTSYS_ELEC = (unique(EndUseAppliances[:,3]))
DISTSYS_GAS  = (unique(EndUseAppliances[:,4]))
DIST_ELEC = length(DISTSYS_ELEC)
DIST_GAS = length(DISTSYS_GAS)

# APP_DistSystemLoc_GAS to tie appliances to gas distribution systems
# APP_DistSystemLoc_ELEC to tie appliances to electric distribution systems
APP_DistSystemLoc_ELEC = zeros(DIST_ELEC, APPLIANCES)
APP_DistSystemLoc_GAS = zeros(DIST_GAS, APPLIANCES)

Loc_ELEC = EndUseAppliances[:,3]
Loc_GAS = EndUseAppliances[:,4]
for a = 1:APPLIANCES
    APP_DistSystemLoc_ELEC[findfirst(occursin.([string(Loc_ELEC[a])],string.(DISTSYS_ELEC))),a] = 1
    APP_DistSystemLoc_GAS[findfirst(occursin.([string(Loc_GAS[a])],string.(DISTSYS_GAS))),a] = 1
end

APPLIANCES_NodalLoc_ELEC = zeros(NODES_ELEC, APPLIANCES)
APPLIANCES_NodalLoc_GAS = zeros(NODES_GAS, APPLIANCES)

Loc_ELEC = EndUseAppliances[:,1]
Loc_GAS = EndUseAppliances[:,2]
for a = 1:APPLIANCES
 APPLIANCES_NodalLoc_ELEC[findfirst(occursin.([string(Loc_ELEC[a])],REGIONS_ELEC)),a] = 1
 APPLIANCES_NodalLoc_GAS[findfirst(occursin.([string(Loc_GAS[a])],REGIONS_GAS)),a] = 1
end



## Gas distribution utility financial assumptions
################################################################################
# For gas distribution retirement evaluation, we need to estimate the potential avoided costs
# of gas system maintenance and reinvestment. To do this, we use the estimated revenue requirement
# and how it evolves across the planning time horizon using a simplified set of assumptions.
RR_est = Costs_GasDistSys     # Each distribution system has an associated total revenue requirement [$/year]
# Here, we assess the potential annual costs of gas system maintenance, depreciation, and reinvestment
# for two cases: business as usual (BAU) and Accelerated Depreciation (AccDep).
BAUGasSyst_FixedCosts = zeros(T_inv)
AccDepGasSyst_FixedCosts = zeros(T_inv,T_inv)
# Using a simplified set of financial assumptions
equity = 0.10               # return on equity afforded to utility shareholders [%]
debt = 0.04                 # interest rate on debt associated with securitization of the gas system  [%]
shareFOM = 0.1              # share of total annual revenue requirement that is fixed operating costs (as opposed to capital investment)
ReinvestmentRate = 0.025    # % of reinvestment 
AvgDepreciation = 0.03      # average % of depreciation per year
println("Estimated Revenue Requirement = $(RR_est) per year")
RB_est = ((1-shareFOM)*RR_est/(equity + AvgDepreciation))
println("Estimated Ratebase = $(RB_est)")
for i = 2:T_inv
    depTimeHorizon = Years[i] - Years[1]
    syd = depTimeHorizon*(depTimeHorizon+1)/2
    nb = RB_est
    depTimeHorizon_remaining = Years[i] - Years[1]
    for j = 1:T_inv-1
        cost = 0
        if j < i
            for y = 1:(Years[j+1]-Years[j])
                cost = cost + depTimeHorizon_remaining/syd*RB_est + nb*debt
                nb = nb - depTimeHorizon_remaining/syd*RB_est
                depTimeHorizon_remaining = depTimeHorizon_remaining - 1
            end
            AccDepGasSyst_FixedCosts[i,j] = cost/(Years[j+1]-Years[j])
        end
        if j >= i
            AccDepGasSyst_FixedCosts[i,j] = 0
        end
    end
    BAUGasSyst_FixedCosts[i] = sum(RR_est*shareFOM + (equity+AvgDepreciation)*((1-AvgDepreciation+ReinvestmentRate)^y)*((1-shareFOM)*RR_est/(equity + AvgDepreciation)) for y = (AllYears[i]-BaseYear+1):(AllYears[i+1]-BaseYear))/(AllYears[i+1]-AllYears[i])
    println("BAU Gas Costs = $(BAUGasSyst_FixedCosts[i])")
    println("ShutDown Gas Costs = $(AccDepGasSyst_FixedCosts[i,:])")
end

BAUGasSyst_FixedCosts[1] = RR_est
# If you retire the gas system in investment period 1, then you must pay off the entire rate base
# in this year
AccDepGasSyst_FixedCosts[1,1] = RB_est



################################################################################
# Transmission interchanges
################################################################################
TransmissionLinks_ELEC = CSV.read(data_dir * "ElecTransmission$(system).csv",DataFrame)
EDGES_ELEC = length(TransmissionLinks_ELEC[:,1])
MAXFLOW_ELEC = TransmissionLinks_ELEC[:,3].*transmission_multiplier
Line_Rating = TransmissionLinks_ELEC[:,4].*transmission_multiplier
Line_Reactance = TransmissionLinks_ELEC[:,5]
ExistingUnits_ElecTrans = TransmissionLinks_ELEC[:,6]
MaxNewUnits_ElecTrans = TransmissionLinks_ELEC[:,7]
CAPEX_ELECTrans = ElecTransmissionCapitalCosts.*Line_Rating
FOM_ELECTrans = ElecTransmissionOperatingCosts.*Line_Rating

TransmissionLinks_GAS = CSV.read(data_dir * "GasTransmission$(system).csv",DataFrame)
EDGES_GAS = length(TransmissionLinks_GAS[:,1])
MAXFLOW_GAS = TransmissionLinks_GAS[:,3]./10
Diameter_Pipes = TransmissionLinks_GAS[:,4]  #[m]
Length_Pipes = TransmissionLinks_GAS[:,5]    #[m]
FrictionFactor_Pipes = TransmissionLinks_GAS[:,6]
ExistingUnits_GasTrans = TransmissionLinks_GAS[:,7]
MaxNewUnits_GasTrans = TransmissionLinks_GAS[:,8]
CompressionRatio_MAX_Branch = TransmissionLinks_GAS[:,9]
CAPEX_GASTrans = GasTransmissionCapitalCosts.*Length_Pipes./1000
FOM_GASTrans = GasTransmissionOperatingCosts.*Length_Pipes./1000

## Parameters for gas pipeline flow simulation
################################################################################
Temp_GAS = 300 # [K]
Temp_N = 298.15  # [K]
Pressure_N = 101325   # [Pa]
pi = 3.14
SpecGravity = 0.64
Compressibility = 0.96
# Pressure is in Pascals, which puts the actual pressure variables in Pa^2
# Compressibility and specific gravity will vary depending on actual injection
# of alternative fuels, but for the purposes of flow evauluation we assume them
# constant to avoid a fully nonlinear problem
K = zeros(size(Diameter_Pipes))
K1 = zeros(size(Diameter_Pipes))
K2 = zeros(size(Diameter_Pipes))
V = zeros(size(Diameter_Pipes))
C = zeros(size(Diameter_Pipes))

# Here, we present three different approaches to the gas flow equation:
# (1) General flow equation
for e = 1:EDGES_GAS
    K[e] = 1/((13.2986*Temp_N/Pressure_N)^2*Diameter_Pipes[e]^5/(Length_Pipes[e]*SpecGravity*Temp_GAS*Compressibility*FrictionFactor_Pipes[e]))
    V[e] = pi/4*Diameter_Pipes[e]^2*Length_Pipes[e]       # m3
    C[e] = V[e]*Temp_N/Pressure_N/Compressibility/Temp_GAS
end
# (2) Weymouth equation
for e = 1:EDGES_GAS
    K[e] = 1/((137.2364*Temp_N/Pressure_N)^2*Diameter_Pipes[e]^5.33/(Length_Pipes[e]*SpecGravity*Temp_GAS*Compressibility))
end
M_CH4 = 16/1000                         # kg/mol
UnivGasConstant = 8.314                 # J/mol-K
UnivDensity_GAS = 1/0.024465            # moles/m3
GasConstant = 8.314/M_CH4               # J/kg-K
Density_GAS = 1/0.024465*M_CH4          # kg/m3 (converted from moles/m3)
# (3) Per Correa-Posada, Carlos M., and Pedro Sanchez-Martin. "Integrated power and natural gas model for energy adequacy in short-term operation." IEEE Transactions on Power Systems 30.6 (2014): 3347-3355.
for e = 1:EDGES_GAS
    K1[e] = (pi/4)*Diameter_Pipes[e]^2/GasConstant/Temp_GAS/Compressibility/Density_GAS
    K2[e] = (pi/4)^2*Diameter_Pipes[e]^5/FrictionFactor_Pipes[e]/GasConstant/Temp_N/Compressibility/Density_GAS^2
end

# Correcting all constants to bring pressures up to MPa
K1 = K1.*10^6
K2 = K2.*10^12
K = K./10^12
C = C.*10^6

################################################################################
### Import set of energy supply/storage/demand units
################################################################################
Generators = CSV.read(data_dir * "Generators$(system).csv",DataFrame)
HourlyVRE2 = CSV.read(data_dir * "HourlyVRE$(system)$(region).csv",DataFrame)
HourlyVRE = zeros(8760,length(HourlyVRE2[1,:]))
for i = 1:length(HourlyVRE2[1,:])
    # Capacity factors must be greater than 0
    HourlyVRE[:,i] = max.(HourlyVRE2[:,i],0)
end

GEN = length(Generators[:, :1])
PrimeMover_GEN = Generators[:,4]
Fuel_GEN = Generators[:,5]
NumUnits_GEN = Generators[:,6]                  # [units]
UnitSize_GEN = Generators[:,7]                  # [MW]
MaxNewUnitsAnnual_GEN = Generators[:,8].*br     # [units/year]
MaxNewUnitsTotal_GEN = Generators[:,9].*br      # [units]
Pmin_GEN = Generators[:,10]                     # [p.u.]
Pmax_GEN = Generators[:,11]                     # [p.u.]
RampDownRate_GEN = Generators[:,12]             # [p.u.]
RampUpRate_GEN = Generators[:,13]               # [p.u.]
MinUpTime_GEN = Generators[:,14]                # [hours]
MinDownTime_GEN = Generators[:,15]              # [hours]
IS_RENEWABLE = Generators[:,16]                 # [bin.]
HeatRate = Generators[:,17]                     # [MMBtu fuel/MWh elec.]
NG_fueled = Generators[:,18]                    # [bin.]
emissions_factors = Generators[:,19]./1000      # [tCO2/MMBtu fuel]
StartUpCosts = Generators[:,20]                 # [$/start]
EconomicLifetime_GEN = Generators[:,21]         # [years]
Lifetime_GEN = Generators[:,22]                 # [years]
StartupFuel = Generators[:,23]                  # [MMBtu/start]
RetirementYear_GEN = min.(Generators[:,24]+Lifetime_GEN,Generators[:,25])
CRF_GEN = (WACC.*(1+WACC).^EconomicLifetime_GEN)./((1+WACC).^EconomicLifetime_GEN .- 1)

PowerToGas = CSV.read(data_dir * "PowerToGas$(system).csv",DataFrame)
P2G = length(PowerToGas[:, :1])
PrimeMover_P2G = PowerToGas[:,4]
NumUnits_P2G = PowerToGas[:,5]                  # [units]
UnitSize_P2G = PowerToGas[:,6]                  # [MW]
MaxNewUnitsAnnual_P2G = PowerToGas[:,7].*br     # [units/year]
MaxNewUnitsTotal_P2G = PowerToGas[:,8].*br      # [units]
Pmin_P2G = PowerToGas[:,9]                      # [p.u.]
Pmax_P2G = PowerToGas[:,10]                     # [p.u.]
RampDownRate_P2G = PowerToGas[:,11]             # [p.u.]
RampUpRate_P2G = PowerToGas[:,12]               # [p.u.]
MinUpTime_P2G = PowerToGas[:,13]                # [hours]
MinDownTime_P2G = PowerToGas[:,14]              # [hours]
eta_P2G = PowerToGas[:,15]                      # [MJ gas/MJ elec.]
eta_P2L = PowerToGas[:,16]                      # [MJ LPG/MJ elec.]
EconomicLifetime_P2G = PowerToGas[:,17]         # [years]
Lifetime_P2G = PowerToGas[:,18]                 # [years]
ISBIOMETHANE = PowerToGas[:,19]                 # [bin.]
ISBIOMASS = PowerToGas[:,20]                    # [bin.]
MoleFracs_P2G = Matrix(PowerToGas[:,23:24])             # [%]
CRF_P2G = (WACC.*(1+WACC).^EconomicLifetime_P2G)./((1+WACC).^EconomicLifetime_P2G .- 1)
RetirementYear_P2G = min.(PowerToGas[:,21]+Lifetime_P2G, PowerToGas[:,22])

ElectricalStorage = CSV.read(data_dir * "Storage_ELEC$(system).csv",DataFrame)
STORAGE_ELEC = length(ElectricalStorage[:, :1])
PrimeMover_STORAGE_ELEC = ElectricalStorage[:,4]
NumUnits_STORAGE_ELEC = ElectricalStorage[:,5]                  # [units]
UnitSize_STORAGE_ELEC = ElectricalStorage[:,6]                  # [MW]
MaxNewUnitsAnnual_STORAGE_ELEC = ElectricalStorage[:,7].*br     # [units/year]
MaxNewUnitsTotal_STORAGE_ELEC = ElectricalStorage[:,8].*br      # [units]
duration_ELEC = ElectricalStorage[:,9]                          # [hours]
eta_charging_ELEC = ElectricalStorage[:,10]                     # [%]
eta_discharging_ELEC = ElectricalStorage[:,11]                  # [%]
eta_loss_ELEC = ElectricalStorage[:,12]                         # [%]
EconomicLifetime_STORAGE_ELEC = ElectricalStorage[:,13]         # [years]
Lifetime_STORAGE_ELEC = ElectricalStorage[:,14]                 # [years]
CRF_STORAGE_ELEC = (WACC.*(1+WACC).^EconomicLifetime_STORAGE_ELEC)./((1+WACC).^EconomicLifetime_STORAGE_ELEC .- 1)
RetirementYear_STORAGE_ELEC = min.(ElectricalStorage[:,15]+Lifetime_STORAGE_ELEC,ElectricalStorage[:,16])

GasStorage = CSV.read(data_dir * "Storage_GAS$(system).csv",DataFrame)
STORAGE_GAS = length(GasStorage[:, :1])
PrimeMover_STORAGE_GAS = GasStorage[:,4]
NumUnits_STORAGE_GAS = GasStorage[:,5]                          # [units]
UnitSize_STORAGE_GAS = GasStorage[:,6]                          # [MW]
MaxNewUnitsAnnual_STORAGE_GAS = GasStorage[:,7]                 # [units/year]
MaxNewUnitsTotal_STORAGE_GAS = GasStorage[:,8]                  # [units]
duration_GAS = GasStorage[:,9]                                  # [hours]
eta_charging_GAS = GasStorage[:,10]                             # [%]
eta_discharging_GAS = GasStorage[:,11]                          # [%]
eta_loss_GAS = GasStorage[:,12]                                 # [%]
EconomicLifetime_STORAGE_GAS = GasStorage[:,13]                 # [years]
Lifetime_STORAGE_GAS = GasStorage[:,14]                         # [years]
CRF_STORAGE_GAS = (WACC.*(1+WACC).^EconomicLifetime_STORAGE_GAS)./((1+WACC).^EconomicLifetime_STORAGE_GAS .- 1)
RetirementYear_STORAGE_GAS = min.(GasStorage[:,15]+Lifetime_STORAGE_GAS,GasStorage[:,16])
MoleFracs_STORAGE = Matrix(GasStorage[:,17:18])

# Gas storage facilities are assumed to be maintained regardless of decisions made in optimization
CAPEX_STORAGE_GAS = 0*ones(T_inv,STORAGE_GAS)                   # [$]
FOM_STORAGE_GAS = 0*ones(T_inv,STORAGE_GAS)                     # [$]


################################################################################
### Gas quality tracking information
################################################################################
GAS_COMPONENTS = 2          # Currently set up for CH4, H2
V_m = 40.87                 # moles/standard m3
MolarMass = [16, 2]         # kg/kmol
LHV = [50, 120]             # MJ/kg
MoleFrac_MAX = [1.0, 0.2]   # kmol/kmol gas
HV_MIN = 40                 # MJ/kg
HV_MAX = 120                # MJ/kg
MoleFracs_SLACK = zeros(NODES_GAS,GAS_COMPONENTS)
MoleFracs_SLACK[:,1] .= 1.0

## For each source of gas, calculate the molar mass [kg/kmol] and LHV [MJ/kg] of gas provided
MolarMass_SLACK = sum(MoleFracs_SLACK.*transpose(MolarMass), dims = 2)         # [kg/kmol gas]
MolarMass_STORAGE = sum(MoleFracs_STORAGE.*transpose(MolarMass), dims = 2)     # [kg/kmol gas]
MolarMass_P2G = sum(MoleFracs_P2G.*transpose(MolarMass), dims = 2)             # [kg/kmol gas]

LHV_SLACK = sum(MoleFracs_SLACK.*transpose(MolarMass.*LHV), dims = 2)./MolarMass_SLACK        # [MJ/kg gas]
LHV_STORAGE = sum(MoleFracs_STORAGE.*transpose(MolarMass.*LHV), dims = 2)./MolarMass_STORAGE    # [MJ/kg gas]
LHV_P2G = sum(MoleFracs_P2G.*transpose(MolarMass.*LHV), dims = 2)./MolarMass_P2G            # [MJ/kg gas]

################################################################################
### CAPEX, FOM, VOM, and fuel costs
################################################################################
CAPEXLookup = CSV.read(data_dir * "CAPEXLookup.csv",DataFrame)
FOMLookup = CSV.read(data_dir * "FOMLookup.csv",DataFrame)
VOMLookup = CSV.read(data_dir * "VOMLookup.csv",DataFrame)
FuelCostLookup = CSV.read(data_dir * "FuelCostLookUp.csv",DataFrame)

CAPEX_GEN = zeros(T_inv,GEN)
FOM_GEN = zeros(T_inv,GEN)
VOM_GEN = zeros(T_inv,GEN)
FuelCosts = zeros(T_inv,GEN)
CAPEX_P2G = zeros(T_inv,P2G)
FOM_P2G = zeros(T_inv,P2G)
VOM_P2G = zeros(T_inv,P2G)
CAPEX_STORAGE_ELEC = zeros(T_inv,STORAGE_ELEC)
FOM_STORAGE_ELEC = zeros(T_inv,STORAGE_ELEC)
CAPEX_APPLIANCES = zeros(T_inv, APPLIANCES)
FOM_APPLIANCES = zeros(T_inv, APPLIANCES)

## Assign the appropriate cost scenario based on CleanElecCosts and CleanGasCosts
################################################################################
CostScenarios = CSV.read(data_dir * "CostScenarios.csv",DataFrame)
if CleanElecCosts =="High"
    if CleanGasCosts == "Low"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_HighElecLowGas.csv",DataFrame)
    end
    if CleanGasCosts == "High"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_HighElecHighGas.csv",DataFrame)
    end
    if CleanGasCosts == "Mid"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_HighElecMidGas.csv",DataFrame)
    end
end

if CleanElecCosts =="Low"
    if CleanGasCosts == "Low"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_LowElecLowGas.csv",DataFrame)
    end
    if CleanGasCosts == "High"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_LowElecHighGas.csv",DataFrame)
    end
    if CleanGasCosts == "Mid"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_LowElecMidGas.csv",DataFrame)
    end
end

if CleanElecCosts =="Mid"
    if CleanGasCosts == "Low"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_MidElecLowGas.csv",DataFrame)
    end
    if CleanGasCosts == "High"
        global CostScenarios = CSV.read(data_dir * "CostScenarios_MidElecHighGas.csv",DataFrame)
    end
    if CleanGasCosts == "Mid"
        global CostScenarios = CSV.read(data_dir * "CostScenarios.csv",DataFrame)
    end
end

# Look up each technology, the associated calendar year in the data tables and assign
# it a cost value
################################################################################
for i = 1:T_inv
    for g = 1:GEN
        subset = findall(in([PrimeMover_GEN[g]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_GEN[i,g] = CAPEXLookup[index, Int(Years[i]-2015)]
        subset = findall(in([PrimeMover_GEN[g]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_GEN[i,g] = FOMLookup[index, Int(Years[i]-2015)]
        subset = findall(in([PrimeMover_GEN[g]]),VOMLookup.Technology)
        scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        VOM_GEN[i,g] = VOMLookup[index, Int(Years[i]-2015)]
        subset = findall(in([Fuel_GEN[g]]),FuelCostLookup.Fuel)
        scen = findall(in([Fuel_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FuelCostLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FuelCosts[i,g] = FuelCostLookup[index, Int(Years[i]-2015)]
    end
    for d = 1:P2G
        subset = findall(in([PrimeMover_P2G[d]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_P2G[i,d] = CAPEXLookup[index, Int(Years[i]-2015)]
        subset = findall(in([PrimeMover_P2G[d]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_P2G[i,d] = FOMLookup[index, Int(Years[i]-2015)]
        subset = findall(in([PrimeMover_P2G[d]]),VOMLookup.Technology)
        scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        VOM_P2G[i,d] = VOMLookup[index, Int(Years[i]-2015)]
    end
    for s = 1:STORAGE_ELEC
        subset = findall(in([PrimeMover_STORAGE_ELEC[s]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_STORAGE_ELEC[s]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_STORAGE_ELEC[i,s] = CAPEXLookup[index, Int(Years[i]-2015)]
        subset = findall(in([PrimeMover_STORAGE_ELEC[s]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_STORAGE_ELEC[s]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_STORAGE_ELEC[i,s] = FOMLookup[index, Int(Years[i]-2015)]
    end
    for a = 1:APPLIANCES
        subset = findall(in([PrimeMover_APPLIANCES[a]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_APPLIANCES[a]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_APPLIANCES[i,a] = CAPEXLookup[index, Int(Years[i]-2015)]
        subset = findall(in([PrimeMover_APPLIANCES[a]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_APPLIANCES[a]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_APPLIANCES[i,a] = FOMLookup[index, Int(Years[i]-2015)]
    end
end

# Separately specify the cost of commodity natural gas for core demands
# using the same assumption applied to gas-fired electricity generators from FuelCosts.csv
################################################################################
CommodityCost_NG = 0.0*ones(T_inv,1)      # $/MWh NG
for i = 1:T_inv
    CommodityCost_NG[i] = sum(FuelCosts[i,:].*NG_fueled)/sum(NG_fueled)/MWh_PER_MMBTU
end


################################################################################
### Clustering time series data for operational simulations
################################################################################
# Here we hope to allow for computationally tractable simulations of grid operations
# by providing the option to cluster the full 8760 hourly profile down into
# a set of representative days. K-mediods is employed.
# Reshape your load vector for your zone of interest into lengths of HOURS_PER_PERIOD

DemandClustering = copy(D_Elec)
for n = 1:NODES_ELEC
    DemandClustering[:,n] = DemandClustering[:,n] +  sum(APPLIANCES_NodalLoc_ELEC[n,a]*ApplianceProfilesELEC[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)
end
DemandClustering = sum(DemandClustering, dims = 2)
DemandClustering = (DemandClustering.- minimum(DemandClustering))./(maximum(DemandClustering) - minimum(DemandClustering))
LOAD_ELEC = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))
DemandClustering = copy(D_Gas)
for n = 1:NODES_GAS
    DemandClustering[:,n] =DemandClustering[:,n] +  sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)
end
DemandClustering = sum(DemandClustering, dims = 2)
DemandClustering = (DemandClustering.- minimum(DemandClustering))./(maximum(DemandClustering) - minimum(DemandClustering))
LOAD_GAS = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))
ProfilesClustering = unique(HourlyVRE, dims = 2)
HourlyVREProfilesClustering = reshape(ProfilesClustering, (HOURS_PER_PERIOD, Periods_Per_Year,length(ProfilesClustering[1,:])))
global ClusteringData = vcat(LOAD_ELEC,  LOAD_GAS)
for x = 1:length(ProfilesClustering[1,:])
    global ClusteringData = vcat(ClusteringData, HourlyVREProfilesClustering[:,:,x])
end

medoids = zeros(T_inv, T_ops)
RepDays = zeros(T_inv, Periods_Per_Year)
weights = zeros(T_inv, T_ops)
m = zeros(length(ClusteringData[:,1]),1)

D = pairwise(Euclidean(), ClusteringData, dims = 2)
H = hclust(D, linkage = :average)
if clustering_case =="ward"
    global H = hclust(D, linkage = :ward)
end
a = cutree(H,k = T_ops)

if clustering_case =="kmeans"
    global H = kmeans(ClusteringData, T_ops)
    global a = assignments(H)
end

for i = 1:T_inv
    for c = 1:(T_ops)
        sub = ClusteringData[:,a .== c]
        FindMedoids = pairwise(Euclidean(),sub,dims = 2)
        Klustering = kmedoids(FindMedoids, 1)
        m[:,1] = sub[:,Klustering.medoids[1]]
        C2 = pairwise(Euclidean(),ClusteringData,m)
        medoids[i,c] = argmin(C2,dims =1)[1][1]
        weights[i,c] = length(sub[1,:])/Periods_Per_Year
    end
    RepDays[i,:] = copy(a)
end


## After identifying the representative days, the raw data are
# re-shaped into the required indexing format for optimization
#################################################################
BaselineDemand_ELEC = zeros(T_inv, T_ops, t_ops, NODES_ELEC)
BaselineDemand_GAS = zeros(T_inv, T_ops, t_ops, NODES_GAS)
HourlyVRE_full = copy(HourlyVRE)
HourlyVRE = zeros(T_ops, t_ops, GEN)
ApplianceProfiles_GAS = zeros(T_ops, t_ops, APPLIANCES)
ApplianceProfiles_ELEC = zeros(T_ops, t_ops, APPLIANCES)
#ApplianceProfiles_LPG = zeros(T_ops, t_ops, APPLIANCES)
for t = 1:T_inv
    for i = 1:T_ops
        start_hour = Int((medoids[1,i]-1)*HOURS_PER_PERIOD+1)
        end_hour = start_hour + HOURS_PER_PERIOD-1
        for n = 1:NODES_ELEC
            BaselineDemand_ELEC[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Elec[start_hour:end_hour,n]
        end
        for n = 1:NODES_GAS
            BaselineDemand_GAS[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Gas[start_hour:end_hour,n]
        end
        for g = 1:GEN
            HourlyVRE[i,:,g] = HourlyVRE_full[start_hour:end_hour,g]
        end
        for a =1:APPLIANCES
            ApplianceProfiles_ELEC[i,:,a] = round.(ApplianceProfilesELEC[start_hour:end_hour,a], digits = 8)
            ApplianceProfiles_GAS[i,:,a] = round.(ApplianceProfilesGAS[start_hour:end_hour,a], digits = 8)
#            ApplianceProfiles_LPG[i,:,a] = round.(ApplianceProfilesLIQ[start_hour:end_hour,a], digits = 8)
        end
    end
end

################################################################################
### Produce residence matrices to identify which nodes contain which resources
################################################################################
# For each set of resources (GEN, DEM, STORAGE_ELEC, STORAGE_GAS, RNG, FS, CCS)
# produce a matrix that is of dimensions [NODES_ELEC x X] and [NODES_GAS x X]
# containing a 1 where each resource resides.
GEN_NodalLoc_ELEC = zeros(NODES_ELEC, GEN)
GEN_NodalLoc_GAS = zeros(NODES_GAS, GEN)
P2G_NodalLoc_ELEC = zeros(NODES_ELEC, P2G)
P2G_NodalLoc_GAS = zeros(NODES_GAS, P2G)
STORAGE_ELEC_NodalLoc_ELEC = zeros(NODES_ELEC, STORAGE_ELEC)
STORAGE_ELEC_NodalLoc_GAS = zeros(NODES_GAS, STORAGE_ELEC)
STORAGE_GAS_NodalLoc_ELEC = zeros(NODES_ELEC, STORAGE_GAS)
STORAGE_GAS_NodalLoc_GAS = zeros(NODES_GAS, STORAGE_GAS)
APPLIANCES_NodalLoc_ELEC = zeros(NODES_ELEC, APPLIANCES)
APPLIANCES_NodalLoc_GAS = zeros(NODES_GAS, APPLIANCES)

NodalLoc_ELEC = Generators[:,1]
NodalLoc_GAS = Generators[:,2]
for g = 1:GEN
    GEN_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[g])],REGIONS_ELEC)),g] = 1
    GEN_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[g])],REGIONS_GAS)),g] = 1
end
NodalLoc_ELEC = PowerToGas[:,1]
NodalLoc_GAS = PowerToGas[:,2]
for d = 1:P2G
 P2G_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[d])],REGIONS_ELEC)),d] = 1
 P2G_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[d])],REGIONS_GAS)),d] = 1
end
NodalLoc_ELEC = ElectricalStorage[:,1]
NodalLoc_GAS = ElectricalStorage[:,2]
for s = 1:STORAGE_ELEC
 STORAGE_ELEC_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[s])],REGIONS_ELEC)),s] = 1
 STORAGE_ELEC_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[s])],REGIONS_GAS)),s] = 1
end
NodalLoc_ELEC = ElectricalStorage[:,1]
NodalLoc_GAS = ElectricalStorage[:,2]
for s = 1:STORAGE_GAS
 STORAGE_GAS_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[s])],REGIONS_ELEC)),s] = 1
 STORAGE_GAS_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[s])],REGIONS_GAS)),s] = 1
end
NodalLoc_ELEC = EndUseAppliances[:,1]
NodalLoc_GAS = EndUseAppliances[:,2]
for a = 1:APPLIANCES
 APPLIANCES_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[a])],REGIONS_ELEC)),a] = 1
 APPLIANCES_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[a])],REGIONS_GAS)),a] = 1
end

################################################################################
### Produce topology matrices to identify which nodes are connected by the edges
################################################################################
# For each Flow (EDGES_ELEC and EDGES_GAS) produce an edge-nodal-incidence matrix A
# If simulating a single-node system, set A to all zeros, there is no transfer outside of the modeled region
A_ELEC = zeros(NODES_ELEC, EDGES_ELEC)
for e = 1:EDGES_ELEC
    A_ELEC[findfirst(occursin.([string(TransmissionLinks_ELEC[e,1])],REGIONS_ELEC)),e] = 1
    A_ELEC[findfirst(occursin.([string(TransmissionLinks_ELEC[e,2])],REGIONS_ELEC)),e] = -1
end
if NODES_ELEC == 1
    A_ELEC = A_ELEC*0
end
A_GAS = zeros(NODES_GAS, EDGES_GAS)
for e = 1:EDGES_GAS
    A_GAS[findfirst(occursin.([string(TransmissionLinks_GAS[e,1])],REGIONS_GAS)),e] = 1
    A_GAS[findfirst(occursin.([string(TransmissionLinks_GAS[e,2])],REGIONS_GAS)),e] = -1
end
if NODES_GAS == 1
    A_GAS = A_GAS*0
end



################################################################################
### Biomethane limitations
################################################################################
# Compute the maximum biomethane and bio-energy use [MWh/year] as a share of initial core gas demands
maxBiomethane = max_biomethane_share*sum(sum(sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES) + D_Gas[:,n] for n = 1:NODES_GAS))*ones(T_inv)               # MWh/yr
# The more generic use of sustainable biomass allows for eventual incorporation of bio-LPG which competes with biomethane for access to limited bioenergy feedstocks
maxSustainableBiomass = max_biomethane_share*sum(sum(sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES) + D_Gas[:,n] for n = 1:NODES_GAS))*ones(T_inv)               # MWh/yr



################################################################################
################################################################################
## Optimization program
################################################################################
################################################################################
m = Model(optimizer_with_attributes(Gurobi.Optimizer,"Threads" => 128,"BarHomogeneous" => 1,"ScaleFlag"=>2, "FeasibilityTol"=> 0.005, "OptimalityTol" => 0.001, "BarConvTol"=> 0.0001, "Method"=> 2, "Crossover"=> 0))

###############################################################################
### Expansion and retirement of energy supply/demand units
# See Eq. 2.5a/2.5b/2.5c in Von Wald thesis
################################################################################
@variable(m, 0 <= unitsbuilt_GEN[I = 1:T_inv, g = 1:GEN]  <= MaxNewUnitsAnnual_GEN[g])                                      # [units]
@variable(m, 0 <= unitsretired_GEN[I = 1:T_inv, g = 1:GEN])                                                                 # [units]
@constraint(m, [g = 1:GEN], sum(unitsbuilt_GEN[i,g] for i = 1:T_inv) <= MaxNewUnitsTotal_GEN[g])                            # [units]

@variable(m, 0 <= unitsbuilt_P2G[I = 1:T_inv, d = 1:P2G]  <= MaxNewUnitsAnnual_P2G[d])                                      # [units]
@variable(m, 0 <= unitsretired_P2G[I = 1:T_inv, d = 1:P2G])                                                                 # [units]
@constraint(m, [d = 1:P2G], sum(unitsbuilt_P2G[i,d] for i = 1:T_inv) <= MaxNewUnitsTotal_P2G[d])                            # [units]

@variable(m, 0 <= unitsbuilt_STORAGE_ELEC[I = 1:T_inv, s = 1:STORAGE_ELEC]  <= MaxNewUnitsAnnual_STORAGE_ELEC[s])           # [units]
@variable(m, 0 <= unitsretired_STORAGE_ELEC[I = 1:T_inv, s = 1:STORAGE_ELEC])                                               # [units]
@constraint(m, [s = 1:STORAGE_ELEC], sum(unitsbuilt_STORAGE_ELEC[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_ELEC[s]) # [units]

@variable(m, 0 <= unitsbuilt_STORAGE_GAS[I = 1:T_inv, s = 1:STORAGE_GAS]  <= MaxNewUnitsAnnual_STORAGE_GAS[s])              # [units]
@variable(m, 0 <= unitsretired_STORAGE_GAS[I = 1:T_inv, s = 1:STORAGE_GAS])                                                 # [units]
@constraint(m, [s = 1:STORAGE_GAS], sum(unitsbuilt_STORAGE_GAS[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_GAS[s])    # [units]

@variable(m, unitsbuilt_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)       # [thousands of units]
@variable(m, unitsremaining_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)   # [thousands of units]
@variable(m, unitsretired_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)     # [thousands of units]


###############################################################################
### Retirement functions for generators, p2g, and storage units
# See Eq. 2.5d/2.5e in Von Wald thesis 
################################################################################
@constraint(m, [I = 1, g = 1:GEN], unitsretired_GEN[I,g] <= NumUnits_GEN[g] + unitsbuilt_GEN[I,g])                          
@constraint(m, [I = 1, g = 1:GEN], unitsretired_GEN[I,g] >= NumUnits_GEN[g]*max(min(Years[I] - RetirementYear_GEN[g],1),0)) 
if T_inv > 1
    @constraint(m, [I = 2:T_inv, g = 1:GEN], unitsretired_GEN[I,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i0,g] - unitsretired_GEN[i0,g]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, g = 1:GEN], unitsretired_GEN[I,g] >= NumUnits_GEN[g]*max(min(Years[I] - RetirementYear_GEN[g],1),0) + sum(unitsbuilt_GEN[i0,g]*max(min(Years[I] - (Years[i0] + Lifetime_GEN[g]),1),0) for i0 = 1:I-1) - sum(unitsretired_GEN[i0,g] for i0 = 1:I-1))
end

@constraint(m, [I = 1, d = 1:P2G], unitsretired_P2G[I,d] <= NumUnits_P2G[d] + unitsbuilt_P2G[I,d])
@constraint(m, [I = 1, d = 1:P2G], unitsretired_P2G[I,d] >= NumUnits_P2G[d]*max(min(Years[I] - RetirementYear_P2G[d],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, d = 1:P2G], unitsretired_P2G[I,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, d = 1:P2G], unitsretired_P2G[I,d] >= NumUnits_P2G[d]*max(min(Years[I] - RetirementYear_P2G[d],1),0) + sum(unitsbuilt_P2G[i0,d]*max(min(Years[I] - (Years[i0] + Lifetime_P2G[d]),1),0) for i0 = 1:I-1) - sum(unitsretired_P2G[i0,d] for i0 = 1:I-1))
end

@constraint(m, [I = 1, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] <= NumUnits_STORAGE_ELEC[s] + unitsbuilt_STORAGE_ELEC[I,s])
@constraint(m, [I = 1, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] >= NumUnits_STORAGE_ELEC[s]*max(min(Years[I] - RetirementYear_STORAGE_ELEC[s],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] <= NumUnits_STORAGE_ELEC[s] + sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] >= NumUnits_STORAGE_ELEC[s]*max(min(Years[I] - RetirementYear_STORAGE_ELEC[s],1),0) + sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min(Years[I] - (Years[i0] + Lifetime_STORAGE_ELEC[s]),1),0) for i0 = 1:I-1)  -  sum(unitsretired_STORAGE_ELEC[i0,s]  for i0 = 1:I-1))
end

@constraint(m, [I = 1, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] <= NumUnits_STORAGE_GAS[s] + unitsbuilt_STORAGE_GAS[I,s])
@constraint(m, [I = 1, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] >= NumUnits_STORAGE_GAS[s]*max(min(Years[I] - RetirementYear_STORAGE_GAS[s],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] <= NumUnits_STORAGE_GAS[s] + sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] >= NumUnits_STORAGE_GAS[s]*max(min(Years[I] - RetirementYear_STORAGE_GAS[s],1),0) + sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min(Years[I] - (Years[i0] + Lifetime_STORAGE_GAS[s]),1),0) for i0 = 1:I-1)  -  sum(unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:I-1))
end


###############################################################################
### Retirement and replacement functions for end-use appliances
# See Eq. 2.10-2.13 in Von Wald thesis
###############################################################################
BaseYear_Sales = zeros(APPLIANCES)
unitsremaining_APPS_historical = zeros(T_inv, APPLIANCES)
for a = 1:APPLIANCES
    BaseYear_Sales[a] = (InitialAppliancePopulation[a]/1000) / sum( max(1 - round(sum(failureProb[a,1:k]),digits = 4),0)/(1+HistoricalGrowthRate[a])^k   for k = 1:50)
    for I = 1:T_inv
        unitsremaining_APPS_historical[I,a] = BaseYear_Sales[a] * sum((1+HistoricalGrowthRate[a])^(-1*k) * (1 - round(sum(failureProb[a,1:k+Int(Years[I] - BaseYear)]),digits = 4))   for k = 1:50)
    end
end

# See Eq. 2.5 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, a = 1:APPLIANCES], unitsremaining_APPS[I,a] == unitsremaining_APPS_historical[I,a] + sum(unitsbuilt_APPS[i0,a] - unitsretired_APPS[i0,a] for i0 = 1:I))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, a = 1:APPLIANCES], unitsretired_APPS[I,a] >= sum(round(cumulativefailurefrac[a,v,I],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1) - sum(unitsretired_APPS[i0,a] for i0 = 1:I-1))
end

# See Eq. 2.14 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, s = 1:SERVICES, d = 1:DIST_GAS], sum(APP_DistSystemLoc_GAS[d,a]*AppliancesToServices[a,s]*(unitsremaining_APPS[I,a]) for a = 1:APPLIANCES) >= (1+ServicesGrowthRate[s])^(Years[I]-BaseYear)*sum(AppliancesToServices[a,s]*APP_DistSystemLoc_GAS[d,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)/1000)
@constraint(m, [I = 1:T_inv, s = 1:SERVICES, d = 1:DIST_ELEC], sum(APP_DistSystemLoc_ELEC[d,a]*AppliancesToServices[a,s]*(unitsremaining_APPS[I,a]) for a = 1:APPLIANCES) >= (1+ServicesGrowthRate[s])^(Years[I]-BaseYear)*sum(AppliancesToServices[a,s]*APP_DistSystemLoc_ELEC[d,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)/1000)

# See Eq. 2.15 in Von Wald thesis
###############################################################################
if appliance_decisions == 0
    @constraint(m, [I = 1:T_inv, a = 1:APPLIANCES], unitsremaining_APPS[I,a] >= ((1+ForecastGrowthRate[a])^(Years[I]-BaseYear))*InitialAppliancePopulation[a]/1000)
end

# See Eq. 2.16 in Von Wald thesis
###############################################################################
@variable(m, Demand_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC] >= 0)
@variable(m, Demand_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS] >= 0)
#@variable(m, Demand_LPG[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS] >= 0)

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], Demand_ELEC[I,T,t,n] == BaselineDemand_ELEC[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_ELEC[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_ELEC[T,t,a] for a = 1:APPLIANCES))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], Demand_GAS[I,T,t,n] == BaselineDemand_GAS[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_GAS[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES))
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], Demand_LPG[I,T,t,n] == 1000*sum(APPLIANCES_NodalLoc_GAS[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_LPG[T,t,a] for a = 1:APPLIANCES))

# To permit sensitivity testing to disallowing hybrid appliance strategies
if hybrids_allowed == 0
    @constraint(m, [a = 1:APPLIANCES, I = 1:T_inv], unitsbuilt_APPS[I,a] <= (1-IS_HYBRID[a])*10^3)
end


###############################################################################
### Dispatch operations of electric generators and gas production
# See Eq. 2.22a in Von Wald thesis
################################################################################
@variable(m, generation[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)         # [MWh]
@variable(m, commit_GEN[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)         # [no. units]
@variable(m, startup_GEN[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)        # [no. units]
@variable(m, shutdown_GEN[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)       # [no. units]

@variable(m, P2G_dispatch[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)       # [MWh]
@variable(m, commit_P2G[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)         # [no. units]
@variable(m, startup_P2G[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)        # [no. units]
@variable(m, shutdown_P2G[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)       # [no. units]


### Startup and shutdown events
# See Eq. 2.22a/2.22e in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], commit_GEN[I,T,t,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], startup_GEN[I,T,t,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], shutdown_GEN[I,T,t,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, g = 1:GEN], commit_GEN[I,T,t,g] == commit_GEN[I,T,t-1,g] + startup_GEN[I,T,t,g] - shutdown_GEN[I,T,t,g])

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], commit_P2G[I,T,t,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i,d] - unitsretired_P2G[i,d] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], startup_P2G[I,T,t,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i,d] - unitsretired_P2G[i,d] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], shutdown_P2G[I,T,t,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i,d] - unitsretired_P2G[i,d] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, d = 1:P2G], commit_P2G[I,T,t,d] == commit_P2G[I,T,t-1,d] + startup_P2G[I,T,t,d] - shutdown_P2G[I,T,t,d])

### Min and Max generation constraints
# See Eq. 2.22b/2.22c in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], generation[I,T,t,g] >= Pmin_GEN[g]*UnitSize_GEN[g]*commit_GEN[I,T,t,g])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], generation[I,T,t,g] <= Pmax_GEN[g]*UnitSize_GEN[g]*commit_GEN[I,T,t,g])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], P2G_dispatch[I,T,t,d] >= Pmin_P2G[d]*UnitSize_P2G[d]*commit_P2G[I,T,t,d])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], P2G_dispatch[I,T,t,d] <= Pmax_P2G[d]*UnitSize_P2G[d]*commit_P2G[I,T,t,d])

### Constraints on fixed profile generation resources
# See Eq. 2.22d in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], generation[I,T,t,g] <= HourlyVRE[T,t,g]*UnitSize_GEN[g]*(NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))) # allows for curtailment
# Explicit calculation of renewable energy curtailments
@variable(m, curtailmentRE[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], curtailmentRE[I,T,t,g] == IS_RENEWABLE[g]*(HourlyVRE[T,t,g]*UnitSize_GEN[g]*(NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I)) - generation[I,T,t,g])) # allows for curtailment

### Ramping constraint
# See Eq. 2.23 in Von Wald thesis
###############################################################################
if t_ops > 1
    minimax = min.(Pmax_GEN, max.(Pmin_GEN,RampDownRate_GEN))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, g = 1:GEN], generation[I,T,t-1,g]-generation[I,T,t,g] <= RampDownRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,T,t,g]-startup_GEN[I,T,t,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*startup_GEN[I,T,t,g] + minimax[g]*UnitSize_GEN[g]*shutdown_GEN[I,T,t,g])
    minimax = min.(Pmax_GEN, max.(Pmin_GEN,RampUpRate_GEN))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, g = 1:GEN], generation[I,T,t,g]-generation[I,T,t-1,g] <= RampUpRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,T,t,g]-startup_GEN[I,T,t,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*shutdown_GEN[I,T,t,g] + minimax[g]*UnitSize_GEN[g]*startup_GEN[I,T,t,g])
end

### Min up time/down time constraints
# See Eq. 2.25/2.26 in Von Wald thesis
###############################################################################
for g = 1:GEN
    if t_ops > MinUpTime_GEN[g]
        for t = Int(MinUpTime_GEN[g]+1):t_ops
            @constraint(m, [I = 1:T_inv, T = 1:T_ops], commit_GEN[I,T,t,g] >= sum(startup_GEN[I,T,t0,g] for t0 = (t-Int(MinUpTime_GEN[g])):t))
        end
    end
end
for g = 1:GEN
    if t_ops > MinDownTime_GEN[g]
        for t = Int(MinDownTime_GEN[g]+1):t_ops
            @constraint(m, [I = 1:T_inv, T = 1:T_ops], NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I) - commit_GEN[I,T,t,g] >= sum(shutdown_GEN[I,T,t0,g] for t0 = (t-MinDownTime_GEN[g]):t))
        end
    end
end

### Generator operational constraints across linked time periods
# See Eq. 2.24/2.25/2.26 in Von Wald thesis
###############################################################################
if LINKED_PERIODS_GENOPS == 1
    for g = 1:GEN
        for i = 2:Int(Periods_Per_Year)
        # Ramping constraint
            minmax = min(Pmax_GEN[g], max(Pmin_GEN[g],RampDownRate_GEN[g]))
            @constraint(m, [I = 1:T_inv], generation[I,Int(RepDays[I,i-1]),t_ops,g]-generation[I,Int(RepDays[I,i]),1,g] <= RampDownRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,Int(RepDays[I,i]),1,g]-startup_GEN[I,Int(RepDays[I,i]),1,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*startup_GEN[I,Int(RepDays[I,i]),1,g] + minmax*UnitSize_GEN[g]*shutdown_GEN[I,Int(RepDays[I,i]),1,g])
            minmax = min(Pmax_GEN[g], max(Pmin_GEN[g],RampUpRate_GEN[g]))
            @constraint(m, [I = 1:T_inv], generation[I,Int(RepDays[I,i]),1,g]-generation[I,Int(RepDays[I,i-1]),t_ops,g] <= RampUpRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,Int(RepDays[I,i]),1,g]-startup_GEN[I,Int(RepDays[I,i]),1,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*shutdown_GEN[I,Int(RepDays[I,i]),1,g] + minmax*UnitSize_GEN[g]*startup_GEN[I,Int(RepDays[I,i]),1,g])
        end
        # Min up time/down time constraints
        firstpd_up = Int(round(MinUpTime_GEN[g]/t_ops)+2)
        X = repeat(collect(range(t_ops,stop = 1,length = t_ops)), outer = [firstpd_up+3])
        Y = repeat(collect(range(0, stop = 24, length = 25)), inner = [t_ops])
        for i = firstpd_up:Int(Periods_Per_Year)
            for t = 1:t_ops
                @constraint(m, [I = 1:T_inv], commit_GEN[I,Int(RepDays[I,i]),t,g] >= sum(startup_GEN[I,Int(RepDays[I,Int(i-Y[Int(t0-t+t_ops)])]),Int(X[Int(t0+t_ops-t)]),g] for t0 = 1:MinUpTime_GEN[g]))
            end
        end
        firstpd_down = Int(round(MinDownTime_GEN[g]/t_ops)+2)
        X = repeat(collect(range(t_ops,stop = 1,length = t_ops)), outer = [firstpd_up+3])
        Y = repeat(collect(range(0, stop = 24, length = 25)), inner = [t_ops])
        for i = firstpd_down:Int(Periods_Per_Year)
            for t = 1:HOURS_PER_PERIOD
                @constraint(m, [I = 1:T_inv], NumUnits_GEN[g] + sum(unitsbuilt_GEN[i0,g] - unitsretired_GEN[i0,g] for i0 = 1:I) - commit_GEN[I,Int(RepDays[I,i]),t,g] >= sum(shutdown_GEN[I,Int(RepDays[I,Int(i-Y[Int(t0-t+t_ops)])]),Int(X[Int(t0+t_ops-t)]),g] for t0 = 1:MinDownTime_GEN[g]))
            end
        end
    end
end


###############################################################################
### Dispatch of storage
# See Eq. 2.51-2.55 in Von Wald thesis
###############################################################################
@variable(m, storedEnergy_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_ELEC] >= 0)       # [MWh]
@variable(m, charging_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC] >= 0)             # [MW]
@variable(m, discharging_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC] >= 0)          # [MW]

@variable(m, storedEnergy_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_GAS] >= 0)         # [MWh]
@variable(m, charging_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS] >= 0)               # [MW]
@variable(m, discharging_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS] >= 0)            # [MW]


@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], storedEnergy_ELEC[I,T,t+1,s] == storedEnergy_ELEC[I,T,t,s] + eta_charging_ELEC[s]*charging_ELEC[I,T,t,s] - (1/eta_discharging_ELEC[s])*discharging_ELEC[I,T,t,s]-eta_loss_ELEC[s]*storedEnergy_ELEC[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_ELEC], storedEnergy_ELEC[I,T,t,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], charging_ELEC[I,T,t,s] <= (1/eta_charging_ELEC[s])*UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], charging_ELEC[I,T,t,s] <= duration_ELEC[s]*UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)) - storedEnergy_ELEC[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], discharging_ELEC[I,T,t,s] <= eta_discharging_ELEC[s]*UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], discharging_ELEC[I,T,t,s] <= storedEnergy_ELEC[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], discharging_ELEC[I,T,t,s]/eta_discharging_ELEC[s] + eta_charging_ELEC[s]*charging_ELEC[I,T,t,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)))

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,t+1,s] == storedEnergy_GAS[I,T,t,s] + eta_charging_GAS[s]*charging_GAS[I,T,t,s] - (1/eta_discharging_GAS[s])*discharging_GAS[I,T,t,s]-eta_loss_GAS[s]*storedEnergy_GAS[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,t,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I))*duration_GAS[s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], charging_GAS[I,T,t,s] <= (1/eta_charging_GAS[s])*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], charging_GAS[I,T,t,s] <= duration_GAS[s]*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)) - storedEnergy_GAS[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s] <= eta_discharging_GAS[s]*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s] <= storedEnergy_GAS[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s]/eta_discharging_GAS[s] + eta_charging_GAS[s]*charging_GAS[I,T,t,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)))

## To limit flexible charge/discharge associated with gas storage
# See Eq. 2.52 in Von Wald thesis
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops-1, s = 1:STORAGE_GAS], charging_GAS[I,T,t,s] == charging_GAS[I,T,t+1,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops-1, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s] == discharging_GAS[I,T,t+1,s])


#@variable(m, storedEnergy_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:2, s = 1:STORAGE_GAS] >= 0)
#@variable(m, charging_GAS[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS] >= 0)
#@variable(m, discharging_GAS[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS] >= 0)
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,t+1,s] == storedEnergy_GAS[I,T,t,s] + t_ops*eta_charging_GAS[s]*charging_GAS[I,T,s] - t_ops*(1/eta_discharging_GAS[s])*discharging_GAS[I,T,s]-eta_loss_GAS[s]*storedEnergy_GAS[I,T,t,s])
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:2, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,t,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I))*duration_GAS[s])
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1, s = 1:STORAGE_GAS], charging_GAS[I,T,s] <= (1/eta_charging_GAS[s])*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I)))
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1, s = 1:STORAGE_GAS], charging_GAS[I,T,s] <= duration_GAS[s]*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I)) - storedEnergy_GAS[I,T,t,s])
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1, s = 1:STORAGE_GAS], discharging_GAS[I,T,s] <= eta_discharging_GAS[s]*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I)))
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1, s = 1:STORAGE_GAS], discharging_GAS[I,T,s] <= storedEnergy_GAS[I,T,t,s])
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1, s = 1:STORAGE_GAS], discharging_GAS[I,T,s]/eta_discharging_GAS[s] + eta_charging_GAS[s]*charging_GAS[I,T,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I)))


### Constraints for linking energy storage across representative periods
# See Eq. 2.56-2.59 in Von Wald thesis
################################################################################
if LINKED_PERIODS_STORAGE == 0
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC], storedEnergy_ELEC[I,T,1,s] == storedEnergy_ELEC[I,T,t_ops+1,s])       # [MWh]
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,1,s] == storedEnergy_GAS[I,T,t_ops+1,s])       # [MWh]
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,1,s] == storedEnergy_GAS[I,T,2,s])       # [MWh]
end
if LINKED_PERIODS_STORAGE == 1
    @variable(m, MinSOC_ELEC[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC] >= 0)
    @variable(m, MaxSOC_ELEC[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC] >= 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC], MinSOC_ELEC[I,T,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC], MaxSOC_ELEC[I,T,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
    @variable(m, SOCTracked_ELEC[I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_ELEC] >= 0)
    @constraint(m, [I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,d,s] <=  UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], MinSOC_ELEC[I,T,s] <= storedEnergy_ELEC[I,T,t,s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], MaxSOC_ELEC[I,T,s] >= storedEnergy_ELEC[I,T,t,s])
    for i = 1:Int(Periods_Per_Year)
        @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,i,s] == storedEnergy_ELEC[I,Int(RepDays[I,1]),1,s] + sum(storedEnergy_ELEC[I,Int(RepDays[I,t]),t_ops+1,s] - storedEnergy_ELEC[I,Int(RepDays[I,t]),1,s] for t = 1:i))
        if i < Periods_Per_Year
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,i,s] + (MaxSOC_ELEC[I,Int(RepDays[I,i+1]),s] - storedEnergy_ELEC[I,Int(RepDays[I,i+1]),1,s]) <=  UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:I))*duration_ELEC[s])
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,i,s] - (storedEnergy_ELEC[I,Int(RepDays[I,i+1]),1,s] - MinSOC_ELEC[I,Int(RepDays[I,i+1]),s]) >= 0)
        end
    end
    @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,Int(Periods_Per_Year),s] == storedEnergy_ELEC[I,Int(RepDays[I,1]),1,s])
    
    @variable(m, MinSOC_GAS[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS] >= 0)
    @variable(m, MaxSOC_GAS[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS] >= 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], MinSOC_GAS[I,T,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I))*duration_GAS[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], MaxSOC_GAS[I,T,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I))*duration_GAS[s])
    @variable(m, SOCTracked_GAS[I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_GAS] >= 0)
    @constraint(m, [I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_GAS], SOCTracked_GAS[I, d, s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I))*duration_GAS[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], MinSOC_GAS[I,T,s] <= storedEnergy_GAS[I,T,t,s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], MaxSOC_GAS[I,T,s] >= storedEnergy_GAS[I,T,t,s])
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:2, s = 1:STORAGE_GAS], MinSOC_GAS[I,T,s] <= storedEnergy_GAS[I,T,t,s])
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:2, s = 1:STORAGE_GAS], MaxSOC_GAS[I,T,s] >= storedEnergy_GAS[I,T,t,s])
    for i = 1:Int(Periods_Per_Year)
#        @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] == storedEnergy_GAS[I,Int(RepDays[I,1]),1,s] + sum(storedEnergy_GAS[I,Int(RepDays[I,t]),2,s] - storedEnergy_GAS[I,Int(RepDays[I,t]),1,s] for t = 1:i))
        @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] == storedEnergy_GAS[I,Int(RepDays[I,1]),1,s] + sum(storedEnergy_GAS[I,Int(RepDays[I,t]),t_ops+1,s] - storedEnergy_GAS[I,Int(RepDays[I,t]),1,s] for t = 1:i))
        if i < Periods_Per_Year
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] + (MaxSOC_GAS[I,Int(RepDays[I,i+1]),s] - storedEnergy_GAS[I,Int(RepDays[I,i+1]),1,s]) <=  UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s] for i0 = 1:I))*duration_GAS[s])
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] - (storedEnergy_GAS[I,Int(RepDays[I,i+1]),1,s] - MinSOC_GAS[I,Int(RepDays[I,i+1]),s]) >= 0)
        end
    end
    @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,Int(Periods_Per_Year),s] == storedEnergy_GAS[I,Int(RepDays[I,1]),1,s])
end




###############################################################################
### Transmission of electric power
# See Eq. 2.18 in Von Wald thesis
###############################################################################
@variable(m, Flows_Elec[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC])

# Built to permit eventual inclusion of transmission expansion
if TRANSMISSION_EXPANSION == 0
    unitsbuilt_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    unitsretired_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
end
if TRANSMISSION_EXPANSION == 1
    @variable(m, unitsbuilt_TRANS_ELEC[I = 1:T_inv, e = 1:EDGES_ELEC], Bin)
    @variable(m, unitsretired_TRANS_ELEC[I = 1:T_inv, e = 1:EDGES_ELEC], Bin)
    @constraint(m, [e = 1:EDGES_ELEC], sum(unitsbuilt_TRANS_ELEC[I,e] for I = 1:T_inv) <=  MaxNewUnits_ElecTrans[e])
    @constraint(m, [I = 1:T_inv, e = 1:EDGES_ELEC], sum(unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) <=  ExistingUnits_ElecTrans[e] + sum(unitsbuilt_TRANS_ELEC[i0,e] for i0 = 1:I))
end

## If not simulating steady-state power flows, flows are only bound by the maximum power flow specified for the line (i.e., simple transport model)
if STEADYSTATE_ELEC == 0
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] <= (sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*MAXFLOW_ELEC[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] >= -1*(sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*MAXFLOW_ELEC[e])
end

## If simulating steady-state power flows, we introduce additional voltage angle variables and constraints governing the power flow
if STEADYSTATE_ELEC == 1
    SLACK_BUS = 1
    @variable(m, -2*pi <= BusAngle[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC] <= 2*pi)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops], BusAngle[I,T,t,SLACK_BUS] == 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] <= BASEMVA*-1/Line_Reactance[e]*sum(A_ELEC[n,e]*BusAngle[I,T,t,n] for n = 1:NODES_ELEC) + MAXFLOW_ELEC[e]*(1-sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e]))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] >= BASEMVA*-1/Line_Reactance[e]*sum(A_ELEC[n,e]*BusAngle[I,T,t,n] for n = 1:NODES_ELEC) - MAXFLOW_ELEC[e]*(1-sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e]))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] <= (sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*Line_Rating[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] >= -1*(sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*Line_Rating[e])
end

###############################################################################
### Energy Balance for electricity grid
# See Eq. 2.19 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], sum(GEN_NodalLoc_ELEC[n,g]*generation[I,T,t,g] for g = 1:GEN) + sum(-1*A_ELEC[n,e]*Flows_Elec[I,T,t,e] for e = 1:EDGES_ELEC) - sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(charging_ELEC[I,T,t,s]-discharging_ELEC[I,T,t,s]) for s = 1:STORAGE_ELEC) - Demand_ELEC[I,T,t,n] - sum(P2G_NodalLoc_ELEC[n,d]*P2G_dispatch[I,T,t,d]*(1-ISBIOMETHANE[d]) for d = 1:P2G) >= 0)


###############################################################################
### Transmission of gaseous fuel
# See Eq. 2.X in Von Wald thesis
###############################################################################
# Nodal gaseous energy demand balance imposed on average across the day
@variable(m, Flows_Gas[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])                      # [standard m3/sec]

if TRANSMISSION_EXPANSION == 0
    unitsbuilt_TRANS_GAS = zeros(T_inv,EDGES_GAS)
    unitsretired_TRANS_GAS = zeros(T_inv,EDGES_GAS)
end
if TRANSMISSION_EXPANSION == 1
    @variable(m, unitsbuilt_TRANS_GAS[I = 1:T_inv, e = 1:EDGES_GAS], Bin)
    @variable(m, unitsretired_TRANS_GAS[I = 1:T_inv, e = 1:EDGES_GAS], Bin)
    @constraint(m, [e = 1:EDGES_GAS], sum(unitsbuilt_TRANS_GAS[I,e] for I = 1:T_inv) <=  MaxNewUnits_GasTrans[e])
    @constraint(m, [I = 1:T_inv, e = 1:EDGES_GAS], sum(unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) <=  ExistingUnits_GasTrans[e] + sum(unitsbuilt_TRANS_GAS[i0,e] for i0 = 1:I))
end

if GASFLOW_DIRECTIONS == 0
    GasFlowDirection = ones(T_inv,T_ops,EDGES_GAS)
end
if GASFLOW_DIRECTIONS == 1
    @variable(m, GasFlowDirection[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Bin)
end

@constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], -1*Flows_Gas[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e] <= GasFlowDirection[I,T,e]*MAXFLOW_GAS[e])
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e] >= (GasFlowDirection[I,T,e]-1)*(MAXFLOW_GAS[e]))

if (STEADYSTATE_GAS == 1) & (NODES_GAS > 1)
    @variable(m, (PRESSURE_MIN) <= NodalPressureSqrd[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] <= (PRESSURE_MAX))
    @variable(m, PRESSURE_MIN <= CompressionPressureSqrd[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd[I,T,e] >= sum(max(A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd[I,T,e] <= CompressionRatio_MAX_Branch[e]*sum(max(A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS))

    @constraint(m, [I = 1:T_inv, T = 1:T_ops], NodalPressureSqrd[I,T,SLACK_NODE] == (SLACK_NODE_Pressure))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS],  CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS) <= GasFlowDirection[I,T,e]*((PRESSURE_MAX)-(PRESSURE_MIN)))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS) >= (1-GasFlowDirection[I,T,e])*((PRESSURE_MIN)-(PRESSURE_MAX)))
    
    
    @variable(m, 0 <= lambda[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX-PRESSURE_MIN)
    a1 = -1
    a2 = PRESSURE_MIN-(PRESSURE_MAX)
    b1 = 1
    b2 = (PRESSURE_MAX)-(PRESSURE_MIN)

    # Lambda is a tight constraint such that when y = 1, lambda = P1-P2 and when y = 0, lambda = P2-P1
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] >= a2*(2*GasFlowDirection[I,T,e]-1) + a1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - a1*a2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] >= b2*(2*GasFlowDirection[I,T,e]-1) + b1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - b1*b2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] <= b2*(2*GasFlowDirection[I,T,e]-1) + a1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - a1*b2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] <= a2*(2*GasFlowDirection[I,T,e]-1) + b1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - b1*a2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e]*Flows_Gas[I,T,e] <= K2[e]/Length_Pipes[e]*(lambda[I,T,e]))
end


# ### Gaseous fuel slack supply
# # See Eq. 2.41 in Von Wald thesis
# ###############################################################################
# @variable(m, SUPPLY_GAS_slack[I = 1:T_inv, T = 1:T_ops,  n = 1:NODES_GAS] >= 0)
# @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack[I,T,n] <= MAXSLACK[n])

# maxPipelineH2ByVol = 0.2
# maxPipelineH2ByEnergy = 0.07

# ### Gas quality constraints
# # See Eq. 2.X in Von Wald thesis
# ###############################################################################
# if GasQuality == "Nodal"
#     @variable(m, Flows_H2[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])
#     @variable(m, LocalH2Consumption[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] >= 0)

#     # Makes sure that flows are in the right direction relative to binary direction variables
#     @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_H2[I,T,e] <= GasFlowDirection[I,T,e]*(MAXFLOW_GAS[e]))
#     @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_H2[I,T,e] >= (GasFlowDirection[I,T,e]-1)*(MAXFLOW_GAS[e]))

#     # Makes sure that flows are in accordance with the maximum volumetric fraction
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_H2[I,T,e] <= maxPipelineH2ByVol*Flows_Gas[I,T,e] + (1-GasFlowDirection[I,T,e])*MAXFLOW_GAS[e])
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_H2[I,T,e] >= maxPipelineH2ByVol*Flows_Gas[I,T,e] - (GasFlowDirection[I,T,e])*MAXFLOW_GAS[e])

#     # Conservation of energy with flows disaggregated between NG and H2
# #    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*(Flows_Gas[I,T,e] - Flows_H2[I,T,e]) for e = 1:EDGES_GAS)) + LHV_H2*(sum(-1*A_GAS[n,e]*(Flows_H2[I,T,e]) for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - sum(Demand_GAS[I,T,t,n] for t = 1:t_ops)/t_ops + sum(sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops - sum(sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops)/t_ops == 0)
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*(Flows_Gas[I,T,e] - Flows_H2[I,T,e]) for e = 1:EDGES_GAS)) + LHV_H2*(sum(-1*A_GAS[n,e]*(Flows_H2[I,T,e]) for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops - sum(Demand_GAS[I,T,t,n] for t = 1:t_ops)/t_ops + sum(sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops - sum(sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops)/t_ops == 0)

#     # Conservation of hydrogen at every node
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_H2*(sum(max(-1*A_GAS[n,e],0)*Flows_H2[I,T,e] for e = 1:EDGES_GAS) - sum(max(A_GAS[n,e],0)*Flows_H2[I,T,e] for e = 1:EDGES_GAS)) - LocalH2Consumption[I,T,n] + sum(sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d]*H2FracProduced[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)
#     # Makes sure that local H2 consumption is less than the max fraction by energy
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LocalH2Consumption[I,T,n] <= maxPipelineH2ByEnergy*sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops)/t_ops)
# end
# if GasQuality == "Annual"
# #    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)
#     # Makes sure that total H2 consumption is less than the max fraction by energy (annual basis)
#     @constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*sum(sum(P2G_dispatch[I,T,t,d]*eta_P2G[d]*H2FracProduced[d] for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops) <= maxPipelineH2ByEnergy*(sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) + sum(weights[I,T]*8760/t_ops*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)))
# end
# if GasQuality == "No"
# #    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)
# end



################################################################################
### Gaseous fuel slack supply
# See Eq. 2.41 in Von Wald thesis
###############################################################################
@variable(m, SUPPLY_GAS_slack[I = 1:T_inv, T = 1:T_ops,  n = 1:NODES_GAS] >= 0)
@constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack[I,T,n] <= MAXSLACK[n])

## Introduce nominal gas component tracking formulation
# See Eq. 2.42 in Von Wald thesis
###############################################################################
@variable(m, NominalGasOfftakes[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS, g = 1:GAS_COMPONENTS] >= 0)     # kmol/sec
@variable(m, NominalGasFlows[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS])                          # kmol/sec
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] <= GasFlowDirection[I,T,e]*(MAXFLOW_GAS[e]))
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] >= (GasFlowDirection[I,T,e]-1)*(MAXFLOW_GAS[e]))

## Constrain nominal values to match physical values
# See Eq. 2.43-2.44 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], sum(MolarMass[g]*NominalGasFlows[I,T,e,g] for g = 1:GAS_COMPONENTS) == M_CH4*V_m*Flows_Gas[I,T,e])
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], sum(LHV[g]*MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) == Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN))

## Molar balance
# See Eq. 2.45 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS, g = 1:GAS_COMPONENTS], MoleFracs_SLACK[n,g]/(LHV_SLACK[n]*MolarMass_SLACK[n])*SUPPLY_GAS_slack[I,T,n] + sum(-1*A_GAS[n,e]*NominalGasFlows[I,T,e,g] for e = 1:EDGES_GAS) - sum(NominalGasOfftakes[I,T,t,n,g] for t = 1:t_ops)/t_ops - sum(sum(MoleFracs_STORAGE[s,g]/(LHV_STORAGE[s]*MolarMass_STORAGE[s])*STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops + sum(sum(MoleFracs_P2G[d,g]/(LHV_P2G[d]*MolarMass_P2G[d])*P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)

## Energy balance
# See Eq. 2.47 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], sum(sum(-1*A_GAS[n,e]*LHV[g]*MolarMass[g]*NominalGasFlows[I,T,e,g] for e = 1:EDGES_GAS) for g = 1:GAS_COMPONENTS) + SUPPLY_GAS_slack[I,T,n]  - sum(sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)


################################################################################
### Gas quality constraints
# See Eq. 2.48-2.50 in Von Wald thesis
################################################################################
if GasQuality == "Annual"
    ## Mole fraction limit imposed on annual, system-wide basis
    # See Eq. 2.48 in Von Wald thesis
    @constraint(m, [I = 1:T_inv, g = 1:GAS_COMPONENTS], sum(weights[i,T]*8760/t_ops*sum(sum(NominalGasOfftakes[I,T,t,n,g] for t = 1:t_ops) for n = 1:NODES_GAS) for T = 1:T_ops) <= MoleFrac_MAX[g]*sum(sum(weights[i,T]*8760/t_ops*sum(sum(NominalGasOfftakes[I,T,t,n,h] for t = 1:t_ops) for n = 1:NODES_GAS) for T = 1:T_ops) for h = 1:GAS_COMPONENTS))
end
if GasQuality == "Nodal"
    ## Mole fraction limit imposed on daily, nodal basis
    # See Eq. 2.49 in Von Wald thesis
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] <= MoleFrac_MAX[g]*sum(NominalGasFlows[I,T,e,h] for h = 1:GAS_COMPONENTS) + (1-GasFlowDirection[I,T,e])*(MAXFLOW_GAS[e]))
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] >= MoleFrac_MAX[g]*sum(NominalGasFlows[I,T,e,h] for h = 1:GAS_COMPONENTS) - GasFlowDirection[I,T,e]*(MAXFLOW_GAS[e]))
    ## Heating value limit imposed on daily, nodal basis
    # See Eq. 2.50 in Von Wald thesis
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, n = 1:NODES_GAS], sum(sum(LHV[g]*MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops) <= HV_MAX*sum(sum(MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops))
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, n = 1:NODES_GAS], sum(sum(LHV[g]*MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops) >= HV_MIN*sum(sum(MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops))
end


###############################################################################
### Bounding steady-states
# Not included in Von Wald thesis
###############################################################################
if bounding_steady_states == 1
    @variable(m, Flows_Gas_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])                  # [m3/sec]
    @variable(m, Flows_Gas_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])                  # [m3/sec]
    @variable(m, SUPPLY_GAS_slack_Max[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] >= 0)
    @variable(m, SUPPLY_GAS_slack_Min[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] >= 0)
    @variable(m, NetGasDemands_Max[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS])
    @variable(m, NetGasDemands_Min[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS])
    
    
    if GASFLOW_DIRECTIONS == 0
        GasFlowDirection_Min = ones(T_inv,T_ops,EDGES_GAS)
        GasFlowDirection_Max = ones(T_inv,T_ops,EDGES_GAS)
    end
    if GASFLOW_DIRECTIONS == 1
        @variable(m, GasFlowDirection_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Bin)
        @variable(m, GasFlowDirection_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Bin)
    end
    
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], -1*Flows_Gas_Min[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e] <= GasFlowDirection_Min[I,T,e]*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e] >= (GasFlowDirection_Min[I,T,e]-1)*(MAXFLOW_GAS[e]))
    
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], -1*Flows_Gas_Max[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e] <= GasFlowDirection_Max[I,T,e]*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e] >= (GasFlowDirection_Max[I,T,e]-1)*(MAXFLOW_GAS[e]))
    
    if (STEADYSTATE_GAS == 1) & (NODES_GAS > 1)
        @variable(m, (PRESSURE_MIN) <= NodalPressureSqrd_Min[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] <= (PRESSURE_MAX))
        @variable(m, PRESSURE_MIN <= CompressionPressureSqrd_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX)
        @variable(m, (PRESSURE_MIN) <= NodalPressureSqrd_Max[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] <= (PRESSURE_MAX))
        @variable(m, PRESSURE_MIN <= CompressionPressureSqrd_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Min[I,T,e] >= sum(max(A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Min[I,T,e] <= CompressionRatio_MAX_Branch[e]*sum(max(A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Max[I,T,e] >= sum(max(A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Max[I,T,e] <= CompressionRatio_MAX_Branch[e]*sum(max(A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS))
    
        @constraint(m, [I = 1:T_inv, T = 1:T_ops], NodalPressureSqrd_Min[I,T,SLACK_NODE] == (SLACK_NODE_Pressure))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS],  CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS) <= GasFlowDirection_Min[I,T,e]*((PRESSURE_MAX)-(PRESSURE_MIN)))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS) >= (1-GasFlowDirection_Min[I,T,e])*((PRESSURE_MIN)-(PRESSURE_MAX)))
    
        @constraint(m, [I = 1:T_inv, T = 1:T_ops], NodalPressureSqrd_Max[I,T,SLACK_NODE] == (SLACK_NODE_Pressure))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS],  CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS) <= GasFlowDirection_Max[I,T,e]*((PRESSURE_MAX)-(PRESSURE_MIN)))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS) >= (1-GasFlowDirection_Max[I,T,e])*((PRESSURE_MIN)-(PRESSURE_MAX)))
        @variable(m, 0 <= lambda_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX-PRESSURE_MIN)
        @variable(m, 0 <= lambda_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX-PRESSURE_MIN)
        a1 = -1
        a2 = PRESSURE_MIN-(PRESSURE_MAX)
        b1 = 1
        b2 = (PRESSURE_MAX)-(PRESSURE_MIN)
    
        # Lambda is a tight constraint such that when y = 1, lambda = P1-P2 and when y = 0, lambda = P2-P1
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] >= a2*(2*GasFlowDirection_Min[I,T,e]-1) + a1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - a1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] >= b2*(2*GasFlowDirection_Min[I,T,e]-1) + b1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - b1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] <= b2*(2*GasFlowDirection_Min[I,T,e]-1) + a1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - a1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] <= a2*(2*GasFlowDirection_Min[I,T,e]-1) + b1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - b1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] >= a2*(2*GasFlowDirection_Max[I,T,e]-1) + a1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - a1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] >= b2*(2*GasFlowDirection_Max[I,T,e]-1) + b1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - b1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] <= b2*(2*GasFlowDirection_Max[I,T,e]-1) + a1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - a1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] <= a2*(2*GasFlowDirection_Max[I,T,e]-1) + b1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - b1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e]*Flows_Gas_Min[I,T,e] <= K2[e]/Length_Pipes[e]*(lambda_Min[I,T,e]))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e]*Flows_Gas_Max[I,T,e] <= K2[e]/Length_Pipes[e]*(lambda_Max[I,T,e]))
    end
    
    
    ###############################################################################
    ### Gaseous fuel Supply = Demand
    ###############################################################################
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack_Min[I,T,n] <= MAXSLACK[n])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack_Max[I,T,n] <= MAXSLACK[n])
    
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Max[I,T,n] >= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Min[I,T,n] <= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Max[I,T,n] >= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Min[I,T,n] <= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
    
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas_Min[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack_Min[I,T,n] + NetGasDemands_Min[I,T,n] == 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas_Max[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack_Max[I,T,n] + NetGasDemands_Max[I,T,n] == 0)

end


################################################################################
### Policy-driven constraints
################################################################################

### Nominal allocation of net-zero emissions gas consumption to each sector, 
# not to exceed the amount of gaseous energy consumed by that sector
# See Eq. 2.62 in Von Wald thesis
################################################################################
@variable(m, CleanGas_gassector[I = 1:T_inv] >= 0)
@variable(m, CleanGas_powersector[I = 1:T_inv] >= 0)
@constraint(m, [I = 1:T_inv], CleanGas_gassector[I] + CleanGas_powersector[I] == sum(weights[I,T]*8760/t_ops*sum(sum(P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], CleanGas_powersector[I] <= sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], CleanGas_gassector[I] <= sum(weights[I,T]*8760/t_ops*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops))

# LPG - not included in Von Wald thesis
################################################################################
# @variable(m, CleanLPG_gassector[I = 1:T_inv] >= 0)
# @constraint(m, [I = 1:T_inv], CleanLPG_gassector[I] == sum(weights[I,T]*8760/t_ops*sum(sum(P2G_dispatch[I,T,t,d]*eta_P2L[d] for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))
# @constraint(m, [I = 1:T_inv], CleanLPG_gassector[I] >= sum(weights[I,T]*8760/t_ops*sum(sum(Demand_LPG[I,T,t,n] for t = 1:t_ops) for n = 1:NODES_GAS) for T = 1:T_ops))
# @constraint(m, [I = 1:T_inv], CleanLPG_gassector[I] <= liquids_allowed*10^4)


### Slack variables for emissions constraint violations that could represent:
#   > use of negative emissions offsets at a fixed cost 
#   > excess carbon emissions evaluated in objective function at a "social cost of carbon"
# The power and gas sector emissions that exceed their respective emissions intensity constraint are constrained by the maxOffsets share of total emissions liabilities
# See Eq. 2.61 in Von Wald thesis
################################################################################
@variable(m, excess_powerEmissions[I = 1:T_inv] >= 0)
@variable(m, excess_gasEmissions[I = 1:T_inv] >= 0)
@constraint(m, [I = 1:T_inv], excess_powerEmissions[I] <= maxOffsets[I]*(sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)))
@constraint(m, [I = 1:T_inv], excess_gasEmissions[I] <= maxOffsets[I]*(sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)))

### Constraints on the emissions intensity of the electric and gas sectors.
# Here, some emissions from gaseous fuel consumption are offset by the nominal allocation of net-zero emission gas to each sector.
# See Eq. 2.60 in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*(CleanGas_powersector[I])  <= EI_ElecSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(generation[I,T,t,g0] for g0 = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) + excess_powerEmissions[I])
@constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[I] <= EI_GasSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) + excess_gasEmissions[I])


### Maximum biomethane production and use of sustainable bio-energy. 
# Total bio-energy constraint is included in the case where net-zero emissions fuel production units can be used to generate methane or LPG fuel
# but must compete for sustainable biomass feedstocks.
# See Eq. 2.63 in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv], maxBiomethane[I] >= sum(weights[I,T]*8760/t_ops*sum(sum(ISBIOMETHANE[d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], maxSustainableBiomass[I] >= sum(weights[I,T]*8760/t_ops*sum(sum(ISBIOMASS[d]*P2G_dispatch[I,T,t,d]*(eta_P2G[d]+eta_P2L[d]) for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))


###############################################################################
### Gas distribution retirement constraint set
# See Eq. 2.69 in Von Wald thesis
###############################################################################
# User can select whether to allow the model to decide when the retire the gas system
# by setting gasdistretirement_allowed = 1.
# In this case, binary variables are introduced to indicate for each distribution system
# when that system is shut down.
if gasdistretirement_allowed == 1
    @variable(m, distSysRetirement_GAS[I = 1:T_inv, d = 1:DIST_GAS], Bin)
    @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(sum(InitialAppliancePopulation[a]/1000*ApplianceProfilesGAS[t,a] for t = 1:8760) for a = 1:APPLIANCES) >= sum(sum(sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES) for t= 1:t_ops) for T = 1:T_ops))
    @constraint(m, [d = 1:DIST_GAS], sum(distSysRetirement_GAS[j,d] for j = 1:T_inv) <= 1)
end
# If no gas distribution retirement is contemplated by the model, the distSysRetirement_GAS indicators are fixed parameters
if gasdistretirement_allowed == 0
    distSysRetirement_GAS = zeros(T_inv,DIST_GAS)
end
# Use may also select whether to force the model to shut down gas distribution systems
# by using gasdistretirement_forced with a 1 in the year where the system must be shut down.
# Currently all distribution systems must be shut down in the same year.
if sum(gasdistretirement_forced) >= 1
    distSysRetirement_GAS = zeros(T_inv,DIST_GAS)
    for ii = 1:T_inv
        distSysRetirement_GAS[ii,:] .= gasdistretirement_forced[ii]
    end
    # And delivered gas volumes are constrained to be 0 during and after this designated shut-down year:
    @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(sum(InitialAppliancePopulation[a]/1000*ApplianceProfilesGAS[t,a] for t = 1:8760) for a = 1:APPLIANCES) >= sum(sum(sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES) for t= 1:t_ops) for T = 1:T_ops))
end

## The gas distribution system fixed costs are then computed based on these shut-down decisions:
# Not explicitly included in Von Wald thesis
###############################################################################
@variable(m, gasdistsyst_Cost[I = 1:T_inv] >= 0)
# Gas distribution system cost includes several terms that will be either active or zero depending on whether the retirement decision is made:
@constraint(m, [I = 1:T_inv], gasdistsyst_Cost[I] == sum(sum(AccDepGasSyst_FixedCosts[j,I]/1000*distSysRetirement_GAS[j,d] for d = 1:DIST_GAS) for j = 1:T_inv) + sum(BAUGasSyst_FixedCosts[I]/1000*(1-sum(distSysRetirement_GAS[j,d] for j = 1:T_inv)) for d = 1:DIST_GAS))


###############################################################################
### Ancillary customer electrification costs 
# See Eq. 2.66 in Von Wald thesis
###############################################################################
@variable(m, applianceInfrastructureCosts[I = 1:T_inv, a = 1:APPLIANCES] >= 0)
@constraint(m,[I = 1, a = 1:APPLIANCES], applianceInfrastructureCosts[I,a] >= unitsbuilt_APPS[I,a]*1000*upgrade_cost[a])
if T_inv > 1
    @constraint(m,[I = 2:T_inv, a = 1:APPLIANCES], applianceInfrastructureCosts[I,a] >= 1000*(unitsbuilt_APPS[I,a] - sum(round(cumulativefailurefrac[a,v,I]-cumulativefailurefrac[a,v,I-1],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1))*upgrade_cost[a])
end

###############################################################################
# Generalized distribution capital costs associated with peak electrical demand
# See Eq. 2.67/2.68 in Von Wald thesis
# Currently implement a few different approaches to compute:
# (a) the total peak power demand at each node
# (b) the peak distribution-level demand at each node
# (c) the peak incremental distribution-level demand due to appliance electrification (i.e., above baseline demand)
# Current version uses PeakDistDemandInc in objective function, but an argument could be made
# that the peak costs should be evaluated with respect to system-wide coincident peak, as opposed to
# the sum of individual nodal peaks.
###############################################################################
@variable(m, PeakDemand[I = 1:T_inv, n = 1:NODES_ELEC] >= 0)
@variable(m, PeakDistDemand[I = 1:T_inv, n = 1:NODES_ELEC] >= 0)
@variable(m, PeakDistDemandInc[I = 1:T_inv, n = 1:NODES_ELEC] >= 0)
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], PeakDemand[I,n] >= Demand_ELEC[I,T,t,n] + sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(charging_ELEC[I,T,t,s]-discharging_ELEC[I,T,t,s]) for s = 1:STORAGE_ELEC) + sum(P2G_NodalLoc_ELEC[n,d]*P2G_dispatch[I,T,t,d]*(1-ISBIOMETHANE[d]) for d = 1:P2G))
@constraint(m, [I = 1:T_inv, t = 1:8760, n = 1:NODES_ELEC], PeakDistDemand[I,n] >= D_Elec[t,n] + 1000*sum(APPLIANCES_NodalLoc_ELEC[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfilesELEC[t,a] for a = 1:APPLIANCES))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], PeakDistDemandInc[I,n] >= Demand_ELEC[I,T,t,n] - BaselineDemand_ELEC[I,T,t,n])


###############################################################################
### Objective function = total societal costs [$/yr]
# See Eq. 2.1 in Von Wald thesis
###############################################################################
@objective(m, Min, sum(discountfactor[i]*(sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*CRF_GEN[g]*CAPEX_GEN[i0,g] for i0 = 1:i) + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*FOM_GEN[i,g] for g = 1:GEN) +  sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])/1000*generation[i,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])/1000*sum(startup_GEN[i,T,t,g] for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s]-unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC) + sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[i],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + applianceInfrastructureCosts[i0,a])  for i0 =1:i)/1000 for a = 1:APPLIANCES) + sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[i],1),0)*CRF_P2G[d]*CAPEX_P2G[i0,d] for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d]-unitsretired_P2G[i0,d] for i0 = 1:i))*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]/1000*P2G_dispatch[i,T,t,d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops) + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[i],1),0)*CRF_STORAGE_GAS[s]*CAPEX_STORAGE_GAS[i0,s] for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s]-unitsretired_STORAGE_GAS[i0,s] for i0 = 1:i))*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS) +  Cost_DistributionInfrastructure*sum(PeakDistDemandInc[i,n] for n = 1:NODES_ELEC) + sum(weights[i,T]*8760/t_ops*CommodityCost_NG[i]/1000*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - CommodityCost_NG[i]*(CleanGas_gassector[i] + CleanGas_powersector[i])/1000 + gasdistsyst_Cost[i] + offsets_Cost[i]/1000*(excess_powerEmissions[i] + excess_gasEmissions[i])) for i = 1:T_inv))

status = optimize!(m)



################################################################################
################################################################################
## Export results for visualization
################################################################################
################################################################################
unitsbuilt_GEN = JuMP.value.(unitsbuilt_GEN)
unitsretired_GEN = JuMP.value.(unitsretired_GEN)
unitsbuilt_STORAGE_ELEC = JuMP.value.(unitsbuilt_STORAGE_ELEC)
unitsretired_STORAGE_ELEC = JuMP.value.(unitsretired_STORAGE_ELEC)
unitsbuilt_P2G = JuMP.value.(unitsbuilt_P2G)
unitsretired_P2G = JuMP.value.(unitsretired_P2G)
unitsbuilt_STORAGE_GAS = JuMP.value.(unitsbuilt_STORAGE_GAS)
unitsretired_STORAGE_GAS = JuMP.value.(unitsretired_STORAGE_GAS)
unitsbuilt_TRANS_GAS = JuMP.value.(unitsbuilt_TRANS_GAS)
unitsretired_TRANS_GAS = JuMP.value.(unitsretired_TRANS_GAS)
unitsbuilt_TRANS_ELEC = JuMP.value.(unitsbuilt_TRANS_ELEC)
unitsretired_TRANS_ELEC = JuMP.value.(unitsretired_TRANS_ELEC)

Demand_GAS = JuMP.value.(Demand_GAS)
Demand_ELEC = JuMP.value.(Demand_ELEC)

if gasdistretirement_allowed == 1
    distSysRetirement_GAS = JuMP.value.(distSysRetirement_GAS)
end

CleanGas_gassector = JuMP.value.(CleanGas_gassector)
CleanGas_powersector = JuMP.value.(CleanGas_powersector)
unitsbuilt_APPS= JuMP.value.(unitsbuilt_APPS)

EmissionsAndCosts = zeros(19,T_inv)
for i = 1:T_inv
    # Terms for computing average electricity and gas rates
    xA = sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*CRF_GEN[g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*1000*CAPEX_GEN[i0,g] for i0 = 1:i) + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*1000*FOM_GEN[i,g] for g = 1:GEN) +  sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*1000*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*1000*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC) + Cost_DistributionInfrastructure*1000*sum(JuMP.value.(PeakDistDemand[i,n]) for n = 1:NODES_ELEC) + offsets_Cost[i]*JuMP.value.(excess_powerEmissions[i]) - CommodityCost_NG[i]*CleanGas_powersector[i]
    xB = CleanGas_powersector[i]
    xC = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
    xD = sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*CRF_P2G[d]*1000*CAPEX_P2G[i0,d]  for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d] for i0 = 1:i))*1000*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]*JuMP.value.(P2G_dispatch[i,T,t,d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    xE = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(P2G_dispatch[i,T,t,d])*(1-ISBIOMETHANE[d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    xF = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(P2G_dispatch[i,T,t,d])*eta_P2G[d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    xG = CommodityCost_NG[i]*(sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(Demand_GAS[i,T,t,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)) + JuMP.value.(gasdistsyst_Cost[i])*1000 + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*CRF_STORAGE_GAS[s]*1000*CAPEX_STORAGE_GAS[i0,s]  for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:i))*1000*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS) + offsets_Cost[i]*JuMP.value.(excess_gasEmissions[i])  - CommodityCost_NG[i]*CleanGas_gassector[i]
    xH = CleanGas_gassector[i]
    xI = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(Demand_GAS[i,T,t,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
    xJ = sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + offsets_Cost[i]*JuMP.value.(excess_powerEmissions[i]) - CommodityCost_NG[i]*CleanGas_powersector[i]
    
    # Average electricity rate
    EmissionsAndCosts[1,i] = (xA*xF+xB*xD)/(xF*xC-xB*xE)
    # Average cost of zero-emission gas
    EmissionsAndCosts[3,i] = (xD+xE*EmissionsAndCosts[1,i])/xF
    # Average gas rate
    EmissionsAndCosts[2,i] = (xG+xH*EmissionsAndCosts[3,i])/xI

    # Average cost of zero-emission gas (exposed to marginal cost of electricity)
    EmissionsAndCosts[19,i] =  (xC*xD+xE*xJ)/(xF*xC-xB*xE)
    # Average electricity rate (assessed w.r.t remaining electricity (not used for P2G) and after the revenues provided by P2G)
    EmissionsAndCosts[17,i] = (xA + xB*EmissionsAndCosts[19,i] - xE*(xJ+xB*EmissionsAndCosts[19,i])/xC)/(xC-xE)
    # Average gas rate
    EmissionsAndCosts[18,i] = (xG+xH*EmissionsAndCosts[19,i])/xI
        
    # Gen Capex
    EmissionsAndCosts[4,i] = sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*CRF_GEN[g]*1000*CAPEX_GEN[i0,g] for i0 = 1:i) for g = 1:GEN)
    # Gen FOM
    EmissionsAndCosts[5,i] = sum(UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*1000*FOM_GEN[i,g] for g = 1:GEN)
    # Gen VOM and fuel
    EmissionsAndCosts[6,i] = sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) - (CommodityCost_NG[i])*CleanGas_powersector[i]
    # Storage ELEC
    EmissionsAndCosts[7,i] = sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*1000*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*1000*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC)
    # T&D
    EmissionsAndCosts[8,i] = Cost_DistributionInfrastructure*1000*sum(JuMP.value.(PeakDistDemand[i,n]) for n = 1:NODES_ELEC)

    # Gas sector costs
    # Commodity
    EmissionsAndCosts[9,i] = CommodityCost_NG[i]*(sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(Demand_GAS[i,T,t,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)) - CommodityCost_NG[i]*CleanGas_gassector[i]
    # Distribution systems
    EmissionsAndCosts[10,i] = JuMP.value.(gasdistsyst_Cost[i])*1000
    # Storage GAS
    EmissionsAndCosts[11,i] = sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[i],1),0)*CRF_STORAGE_GAS[s]*1000*CAPEX_STORAGE_GAS[i0,s]  for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:i))*1000*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS)

    # Appliances
    EmissionsAndCosts[12,i] = sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[i],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + JuMP.value.(applianceInfrastructureCosts[i0,a])) for i0 =1:i) for a = 1:APPLIANCES)

    # P2G costs (without electricity costs, these are included in the electricity generation sectoral costs)
    EmissionsAndCosts[13,i] = sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[i],1),0)*CRF_P2G[d]*1000*CAPEX_P2G[i0,d] for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d] for i0 = 1:i))*1000*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum((VOM_P2G[i,d])*JuMP.value.(P2G_dispatch[i,T,t,d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)

    # Negative emissions offsets
    EmissionsAndCosts[14,i] = offsets_Cost[i]*(JuMP.value.(excess_powerEmissions[i]) + JuMP.value.(excess_gasEmissions[i]))
    
    # Emissions intensity of electricity generated and gas delivered (check to make sure constraint is satisfied)
    EmissionsAndCosts[15,i] = (sum(weights[i,T]*8760/t_ops*sum(StartupFuel[g]*emissions_factors[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g])*HeatRate[g]*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_powersector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
    EmissionsAndCosts[16,i] = (sum(weights[i,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
end

CapacityBuilt = zeros(T_inv,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1)
for i = 1:T_inv
    CapacityBuilt[i,1:GEN] = unitsbuilt_GEN[i,:].*UnitSize_GEN
    CapacityBuilt[i,GEN+1:GEN+STORAGE_ELEC] = unitsbuilt_STORAGE_ELEC[i,:].*UnitSize_STORAGE_ELEC
    CapacityBuilt[i,GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = unitsbuilt_P2G[i,:].*UnitSize_P2G
    CapacityBuilt[i,GEN+STORAGE_ELEC+P2G+1:GEN+STORAGE_ELEC+P2G+STORAGE_GAS] = unitsbuilt_STORAGE_GAS[i,:].*UnitSize_STORAGE_GAS
    CapacityBuilt[i,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1] = sum(distSysRetirement_GAS[i,d] for d = 1:DIST_GAS)
end

CapacityRetired = zeros(T_inv,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1)
for i = 1:T_inv
    CapacityRetired[i,1:GEN] = unitsretired_GEN[i,:].*UnitSize_GEN
    CapacityRetired[i,GEN+1:GEN+STORAGE_ELEC] = unitsretired_STORAGE_ELEC[i,:].*UnitSize_STORAGE_ELEC
    CapacityRetired[i,GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = unitsretired_P2G[i,:].*UnitSize_P2G
    CapacityRetired[i,GEN+STORAGE_ELEC+P2G+1:GEN+STORAGE_ELEC+P2G+STORAGE_GAS] = unitsretired_STORAGE_GAS[i,:].*UnitSize_STORAGE_GAS
    CapacityRetired[i,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1] = sum(distSysRetirement_GAS[i,d] for d = 1:DIST_GAS)
end

GenerationSave = zeros(T_inv,GEN+P2G+8)
for i = 1:T_inv
     GenerationSave[i,1:GEN] = sum(weights[i,T]*8760/t_ops*sum(JuMP.value.(generation[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+1:GEN+P2G] = sum(weights[i,T]*8760/t_ops*sum(JuMP.value.(P2G_dispatch[i,T,t,:]).*eta_P2G for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+1] = sum(weights[i,T]*8760/t_ops*sum(sum(InitialAppliancePopulation[:].*ApplianceProfiles_GAS[T,t,:]) + sum(BaselineDemand_GAS[i,T,t,:])  for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+2] = sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+3] = CleanGas_powersector[i]
     GenerationSave[i,GEN+P2G+4] = CleanGas_gassector[i]
     GenerationSave[i,GEN+P2G+5] = sum(weights[i,T]*8760/t_ops*sum(sum(InitialAppliancePopulation[:].*ApplianceProfiles_ELEC[T,t,:]) + sum(BaselineDemand_ELEC[i,T,t,:])  for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+6] = sum(weights[i,T]*8760/t_ops*sum(sum(Demand_ELEC[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+7] = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(SUPPLY_GAS_slack[i,T,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+8] = sum(weights[i,T]*8760/t_ops*sum(sum(sum(GEN_NodalLoc_ELEC[n,g]*JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) + sum(-1*A_ELEC[n,e]*JuMP.value.(Flows_Elec[i,T,t,e]) for e = 1:EDGES_ELEC) - sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(JuMP.value.(charging_ELEC[i,T,t,s])-JuMP.value.(discharging_ELEC[i,T,t,s])) for s = 1:STORAGE_ELEC) - Demand_ELEC[i,T,t,n] - sum(P2G_NodalLoc_ELEC[n,d]*JuMP.value.(P2G_dispatch[i,T,t,d])*(1-ISBIOMETHANE[d]) for d = 1:P2G) for n = 1:NODES_ELEC)  for t = 1:t_ops) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(curtailmentRE[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
end

HourlyGenFullSave = zeros(T_inv*8760,GEN+STORAGE_ELEC+P2G+1)
HourlyLoadFullSave = zeros(T_inv*8760,4)
HourlySOCFullSave = zeros(T_inv*8760,STORAGE_ELEC)
HourlyTransmissionFullSave = zeros(T_inv*8760,EDGES_ELEC)
HourlyGasSOCFullSave = zeros(T_inv*8760,STORAGE_GAS)
DailyGasTransmissionFullSave = zeros(T_inv*Periods_Per_Year,EDGES_GAS)
for i = 1:T_inv
    for c = 1:Periods_Per_Year
        j = Int(RepDays[i,c])
        count = Int((i-1)*8760+(c-1)*t_ops)+1
        HourlyGenFullSave[count:count+t_ops-1,1:GEN] = JuMP.value.(generation[i,j,:,:])
        HourlyGenFullSave[count:count+t_ops-1, GEN+1:GEN+STORAGE_ELEC] = (JuMP.value.(charging_ELEC[i,j,:,:])-JuMP.value.(discharging_ELEC[i,j,:,:]))
        HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = JuMP.value.(P2G_dispatch[i,j,:,:].*transpose(ones(P2G)-ISBIOMETHANE))
        HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+P2G+1] = sum(JuMP.value.(curtailmentRE[i,j,:,:]), dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 1] = sum(Demand_ELEC[i,j,:,:], dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 2] = sum(Demand_GAS[i,j,:,:], dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 3] = sum(BaselineDemand_ELEC[i,j,:,:], dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 4] = sum(BaselineDemand_GAS[i,j,:,:], dims = 2)
        HourlySOCFullSave[count:count+t_ops-1,:] = JuMP.value.(storedEnergy_ELEC[i,j,1:t_ops,:])
        HourlyGasSOCFullSave[count:count+t_ops-1,:] = JuMP.value.(storedEnergy_GAS[i,j,1:t_ops,:])
        HourlyTransmissionFullSave[count:count+t_ops-1,:] = JuMP.value.(Flows_Elec[i,j,:,:])
        DailyGasTransmissionFullSave[Int((i-1)*Periods_Per_Year + c),:] = JuMP.value.(Flows_Gas[i,j,:])
    end
end

ApplianceDecisions = zeros(3*T_inv,APPLIANCES)
for i = 1:T_inv
    count = (3*i-2)
    ApplianceDecisions[count,:] = unitsbuilt_APPS[i,:]
    ApplianceDecisions[count+1,:] = JuMP.value.(unitsretired_APPS[i,:])
    ApplianceDecisions[count+2,:] = JuMP.value.(unitsremaining_APPS[i,:])
end


CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityApplianceResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(ApplianceDecisions'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityEmissionsandCostResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(EmissionsAndCosts'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityCapacityBuildResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(CapacityBuilt'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityCapacityRetiredResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(CapacityRetired'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityGenerationResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(GenerationSave'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityHourlyGensFullResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(HourlyGenFullSave'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityHourlyLoadFullBaseline_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(HourlyLoadFullSave'), writeheader = true)
CSV.write(out_dir * "$(system)$(region)$(num)$(biomethane)biomethane$(industrials)industrials$(buildingretrofits)buildingretrofits$(GasQuality)GasQualityHourlyTransmissionFullResults_$(case)$(offsets_case)$(retirements_case)$(CleanElecCosts)CostElec$(CleanGasCosts)CostGas$(NETSCost)NETsCost.csv",Tables.table(HourlyTransmissionFullSave'), writeheader = true)

print("Success!")

#### Outdated code for simulating ex-post system (to check for infeasible system designs based on clustered representative days)
# ################################################################################
# ## Simulate the planned system for the final investment year
# ################################################################################
# maxOffsets = 1.0*ones(T_inv)                # % of total emissions
# offsets_Cost = 5000                         # $/tCO2e
# VOLL = 40000                                # $/MWh

# # Setting the number of units in operation at the end of the investment period
# NumUnits_GEN = NumUnits_GEN + max.(transpose(sum(unitsbuilt_GEN - unitsretired_GEN,dims = 1)),-1*NumUnits_GEN)
# NumUnits_STORAGE_ELEC = NumUnits_STORAGE_ELEC + max.(transpose(sum(unitsbuilt_STORAGE_ELEC - unitsretired_STORAGE_ELEC,dims = 1)),-1*NumUnits_STORAGE_ELEC)
# NumUnits_STORAGE_GAS = NumUnits_STORAGE_GAS + max.(transpose(sum(unitsbuilt_STORAGE_GAS - unitsretired_STORAGE_GAS,dims = 1)), -1*NumUnits_STORAGE_GAS)
# NumUnits_P2G = NumUnits_P2G + max.(transpose(sum(unitsbuilt_P2G - unitsretired_P2G,dims = 1)), -1*NumUnits_P2G)

# # Pull remaining appliance population and set up demands as fixed with full clustering
# AppliancePop = 1000 .*(JuMP.value.(remainingAppliances[T_inv,:]) + (InitialApplianceCount[:,T_inv]))

# # Set emissions intensity target to the final investment period constraint value:
# EI_ElecSector = [EIs_ElecSector[T_inv]]
# EI_GasSector = [EIs_GasSector[T_inv]]

# # Simulating a single year of operations
# T_inv = 1
# T_ops = 365
# t_ops = 24
# HOURS_PER_PERIOD = 24

# InvYears = [2040]
# Years = [2040]


# DemandClustering = copy(D_Elec)
# for n = 1:NODES_ELEC
#     DemandClustering[:,n] = DemandClustering[:,n] + sum(APPLIANCES_NodalLoc_ELEC[n,a]*ApplianceProfilesELEC[:,a]*AppliancePop[a] for a = 1:APPLIANCES)
# end
# DemandClustering = sum(DemandClustering, dims = 2)
# DemandClustering = (DemandClustering.- minimum(DemandClustering))./(maximum(DemandClustering) - minimum(DemandClustering))
# LOAD_ELEC = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))

# DemandClustering = copy(D_Gas)
# for n = 1:NODES_GAS
#     DemandClustering[:,n] = DemandClustering[:,n] + sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*AppliancePop[a] for a = 1:APPLIANCES)
# end
# DemandClustering = sum(DemandClustering, dims = 2)
# DemandClustering = (DemandClustering.- minimum(DemandClustering))./(maximum(DemandClustering) - minimum(DemandClustering))
# LOAD_GAS = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))

# ProfilesClustering = unique(FixedCF, dims = 2)
# FixedCFProfilesClustering = reshape(ProfilesClustering, (HOURS_PER_PERIOD, Periods_Per_Year,length(ProfilesClustering[1,:])))

# ClusteringData = vcat(LOAD_ELEC,  LOAD_GAS)
# for x = 1:length(ProfilesClustering[1,:])
#     global ClusteringData = vcat(ClusteringData, FixedCFProfilesClustering[:,:,x])
# end

# medoids = zeros(T_inv, T_ops)
# RepDays = zeros(T_inv, Periods_Per_Year)
# weights = zeros(T_inv, T_ops)
# m = zeros(length(ClusteringData[:,1]),1)


# D = pairwise(Euclidean(), ClusteringData, dims = 2)
# H = hclust(D, linkage = :average)
# if clustering_case =="ward"
#     global H = hclust(D, linkage = :ward)
# end
# a = cutree(H,k = T_ops)
# if clustering_case =="kmeans"
#     global H = kmeans(ClusteringData, T_ops)
#     global a = assignments(H)
# end
# for i = 1:T_inv
#     for c = 1:(T_ops)
#         sub = ClusteringData[:,a .== c]
#         FindMedoids = pairwise(Euclidean(),sub,dims = 2)
#         Klustering = kmedoids(FindMedoids, 1)
#         m[:,1] = sub[:,Klustering.medoids[1]]
#         C2 = pairwise(Euclidean(),ClusteringData,m)
#         medoids[i,c] = argmin(C2,dims =1)[1][1]
#         weights[i,c] = length(sub[1,:])/Periods_Per_Year
#     end
#     RepDays[i,:] = copy(a)
# end

# Demand_ELEC = zeros(T_inv, T_ops, t_ops, NODES_ELEC)
# Demand_GAS = zeros(T_inv, T_ops, t_ops, NODES_GAS)
# Demand_LPG = zeros(T_inv, T_ops, t_ops, NODES_GAS)
# Fixed_CF = zeros(T_ops, t_ops, GEN)
# for t = 1:T_inv
#     for i = 1:T_ops
#         start_hour = Int((medoids[1,i]-1)*HOURS_PER_PERIOD+1)
#         end_hour = start_hour + HOURS_PER_PERIOD-1
#         for n = 1:NODES_ELEC
#             Demand_ELEC[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Elec[start_hour:end_hour,n] +  sum(APPLIANCES_NodalLoc_ELEC[n,a]*ApplianceProfilesELEC[start_hour:end_hour,a]*AppliancePop[a] for a = 1:APPLIANCES)
#         end
#         for n = 1:NODES_GAS
#             Demand_GAS[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Gas[start_hour:end_hour,n] +  sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[start_hour:end_hour,a]*AppliancePop[a] for a = 1:APPLIANCES)
#         end
#         for g = 1:GEN
#             Fixed_CF[i,:,g] = FixedCF[start_hour:end_hour,g]
#         end
#     end
# end

# BaselineDemand_ELEC = zeros(T_inv, T_ops, t_ops, NODES_ELEC)
# BaselineDemand_GAS = zeros(T_inv, T_ops, t_ops, NODES_GAS)
# ApplianceProfiles_GAS = zeros(T_ops, t_ops, APPLIANCES)
# ApplianceProfiles_ELEC = zeros(T_ops, t_ops, APPLIANCES)
# ApplianceProfiles_LPG = zeros(T_ops, t_ops, APPLIANCES)
# for t = 1:T_inv
#     for i = 1:T_ops
#         start_hour = Int((medoids[1,i]-1)*HOURS_PER_PERIOD+1)
#         end_hour = start_hour + HOURS_PER_PERIOD-1
#         for n = 1:NODES_ELEC
#             BaselineDemand_ELEC[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Elec[start_hour:end_hour,n]
#         end
#         for n = 1:NODES_GAS
#             BaselineDemand_GAS[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Gas[start_hour:end_hour,n]
#         end
#         for a =1:APPLIANCES
#             ApplianceProfiles_ELEC[i,:,a] = ApplianceProfilesELEC[start_hour:end_hour,a]
#             ApplianceProfiles_GAS[i,:,a] = ApplianceProfilesGAS[start_hour:end_hour,a]
#             ApplianceProfiles_LPG[i,:,a] = ApplianceProfilesLIQ[start_hour:end_hour,a]
#         end
#     end
# end

# ################################################################################
# ## Update costs specific to 2040
# ################################################################################
# VOM_GEN = zeros(T_inv,GEN)
# FuelCosts = zeros(T_inv,GEN)
# VOM_P2G = zeros(T_inv,P2G)

# CostScenarios = CSV.read("CostScenarios.csv",DataFrame)
# if costscen =="HighCost"
#     global CostScenarios = CSV.read("CostScenarios_High.csv",DataFrame)
# end
# if costscen =="LowCost"
#     global CostScenarios = CSV.read("CostScenarios_Low.csv",DataFrame)
# end

# for i = 1:T_inv
#     for g = 1:GEN
#         subset = findall(in([PrimeMover_GEN[g]]),VOMLookup.Technology)
#         scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
#         scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
#         index = findall(in(subset),scen)
#         index = scen[index[1]]
#         VOM_GEN[i,g] = VOMLookup[index, Int(Years[i]-2015)]
#         subset = findall(in([Fuel_GEN[g]]),FuelCostLookup.Fuel)
#         scen = findall(in([Fuel_GEN[g]]),CostScenarios.Technology)
#         scen = findall(in([CostScenarios.Cost[scen[1]]]),FuelCostLookup.Cost)
#         index = findall(in(subset),scen)
#         index = scen[index[1]]
#         FuelCosts[i,g] = FuelCostLookup[index, Int(Years[i]-2015)]
#     end
#     for d = 1:P2G
#         subset = findall(in([PrimeMover_P2G[d]]),VOMLookup.Technology)
#         scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
#         scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
#         index = findall(in(subset),scen)
#         index = scen[index[1]]
#         VOM_P2G[i,d] = VOMLookup[index, Int(Years[i]-2015)]
#     end
# end




# ################################################################################
# ## Coupled optimization with only operational variables (all investment decisions are fixed)
# ################################################################################


# @variable(m, LostLoad_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC] >= 0)
# @variable(m, LostLoad_GAS[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] >= 0)
# if GasQuality == 0
#     @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LostLoad_GAS[I,T,n] + LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack[I,T,n]  - sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops >= 0)
# end
# # Energy Balance for electricity grid
# @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], LostLoad_ELEC[I,T,t,n] + sum(GEN_NodalLoc_ELEC[n,g]*generation[I,T,t,g] for g = 1:GEN) + sum(-1*A_ELEC[n,e]*Flows_Elec[I,T,t,e] for e = 1:EDGES_ELEC) - sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(charging_ELEC[I,T,t,s]-discharging_ELEC[I,T,t,s]) for s = 1:STORAGE_ELEC) - Demand_ELEC[I,T,t,n] - sum(P2G_NodalLoc_ELEC[n,d]*P2G_dispatch[I,T,t,d]*(1-ISBIOMETHANE[d]) for d = 1:P2G) >= 0)
# @objective(m, Min, sum(VOLL/1000*(sum(weights[i,T]*8760/t_ops*sum(sum(LostLoad_GAS[i,T,n] for t= 1:t_ops) for n = 1:NODES_GAS) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(LostLoad_ELEC[i,T,t,n] for t= 1:t_ops) for n = 1:NODES_ELEC) for T = 1:T_ops)) +  sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*(1-NG_fueled[g])*FuelCosts[i,g])/1000*generation[i,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*(1-NG_fueled[g])*FuelCosts[i,g])/1000*sum(startup_GEN[i,T,t,g] for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]/1000*P2G_dispatch[i,T,t,d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(CommodityCost_NG[i]/1000*sum(SUPPLY_GAS_slack[i,T,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) + offsets_Cost/1000*(excess_powerEmissions[i] + excess_gasEmissions[i]) for i = 1:T_inv))


# ################################################################################
# ## Export results
# ################################################################################

# CleanGas_gassector = JuMP.value.(CleanGas_gassector)
# CleanGas_powersector = JuMP.value.(CleanGas_powersector)

# EmissionsAndCosts = zeros(16,T_inv)
# EmissionsAndCosts[3,:] = CleanGas_powersector
# EmissionsAndCosts[4,:] = CleanGas_gassector
# EmissionsAndCosts[5,:] = JuMP.value.(excess_powerEmissions)
# EmissionsAndCosts[6,:] = JuMP.value.(excess_gasEmissions)

# for i = 1:T_inv
#     EmissionsAndCosts[1,i] = (sum(weights[i,T]*8760/t_ops*sum(StartupFuel[g]*emissions_factors[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g])*HeatRate[g]*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_powersector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
#     EmissionsAndCosts[2,i] = (sum(weights[i,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
#     # Gen
#     EmissionsAndCosts[7,i] = sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*CRF_GEN[g]*1000*CAPEX_GEN[i0,g] for i0 = 1:i) + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*1000*FOM_GEN[i,g] for g = 1:GEN) +  sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*(1-NG_fueled[g])*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*(1-NG_fueled[g])*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops)
#     # Storage
#     EmissionsAndCosts[8,i] = sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*CRF_STORAGE_ELEC[s]*1000*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*1000*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC) + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*CRF_STORAGE_GAS[s]*1000*CAPEX_STORAGE_GAS[i0,s]  for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:i))*1000*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS)
#     # T&D
#     EmissionsAndCosts[9,i] = Cost_DistributionInfrastructure*1000*sum(JuMP.value.(PeakDistDemand[i,n]) for n = 1:NODES_ELEC)
#     # Gas sector costs
#     # Commodity
#     EmissionsAndCosts[10,i] = CommodityCost_NG[i]*(sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(SUPPLY_GAS_slack[i,T,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - (0.293*sum(weights[i,T]*8760/t_ops*sum(NG_fueled[g]*StartupFuel[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + 0.293*sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(NG_fueled[g]*generation[i,T,t,g])*HeatRate[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - CleanGas_powersector[i]))
#     # Distribution systems
#     EmissionsAndCosts[11,i] = sum(GasDistSysCosts[d]*(1-distSysRetirement_GAS[i,d]) for d = 1:DIST_GAS)
#     # Appliances
#     EmissionsAndCosts[12,i] = sum(sum(CRF_APPLIANCES[a]*CAPEX_APPLIANCES[i0,a]*InitialApplianceCount[a,i0]*1000 for a = 1:APPLIANCES) for i0 =1:i)
#     # P2G costs
#     EmissionsAndCosts[13,i] = sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*CRF_P2G[d]*1000*CAPEX_P2G[i0,d]  for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d] for i0 = 1:i))*1000*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]*JuMP.value.(P2G_dispatch[i,T,t,d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
#     # Negative emissions offsets
#     EmissionsAndCosts[14,i] = offsets_Cost*(JuMP.value.(excess_powerEmissions[i]) + JuMP.value.(excess_gasEmissions[i]))/1000
#     # Average electricity rate
#     EmissionsAndCosts[15,i] = (CommodityCost_NG[i]*(sum(weights[i,T]*8760/t_ops*sum(NG_fueled[g]*StartupFuel[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(NG_fueled[g]*generation[i,T,t,g])*HeatRate[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - CleanGas_powersector[i]) + EmissionsAndCosts[7,i] +  EmissionsAndCosts[8,i] + EmissionsAndCosts[9,i] + CleanGas_powersector[i]/(CleanGas_gassector[i] + CleanGas_powersector[i])*EmissionsAndCosts[13,i] + offsets_Cost*JuMP.value.(excess_powerEmissions[i])) / sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
#     # Average gas rate
#     EmissionsAndCosts[16,i] = (CleanGas_gassector[i]/(CleanGas_gassector[i] + CleanGas_powersector[i])*EmissionsAndCosts[15,i]*sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(P2G_dispatch[i,T,t,d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops) + EmissionsAndCosts[10,i] + EmissionsAndCosts[11,i] + CleanGas_gassector[i]/(CleanGas_gassector[i] + CleanGas_powersector[i])*EmissionsAndCosts[13,i] + offsets_Cost*JuMP.value.(excess_gasEmissions[i]) ) /sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
# end


# GenerationSave = zeros(T_inv,GEN+P2G+8)
# for i = 1:T_inv
#      GenerationSave[i,1:GEN] = sum(weights[i,T]*8760/t_ops*sum(JuMP.value.(generation[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+1:GEN+P2G] = sum(weights[i,T]*8760/t_ops*sum(JuMP.value.(P2G_dispatch[i,T,t,:]).*eta_P2G for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+P2G+1] = sum(weights[i,T]*8760/t_ops*sum(sum(InitialAppliancePopulation[:].*ApplianceProfiles_GAS[T,t,:]) + sum(BaselineDemand_GAS[i,T,t,:])  for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+P2G+2] = sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+P2G+3] = CleanGas_powersector[i]
#      GenerationSave[i,GEN+P2G+4] = CleanGas_gassector[i]
#      GenerationSave[i,GEN+P2G+5] = sum(weights[i,T]*8760/t_ops*sum(sum(InitialAppliancePopulation[:].*ApplianceProfiles_ELEC[T,t,:]) + sum(BaselineDemand_ELEC[i,T,t,:])  for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+P2G+6] = sum(weights[i,T]*8760/t_ops*sum(sum(Demand_ELEC[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+P2G+7] = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(SUPPLY_GAS_slack[i,T,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
#      GenerationSave[i,GEN+P2G+8] = sum(weights[i,T]*8760/t_ops*sum(sum(sum(GEN_NodalLoc_ELEC[n,g]*JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) + sum(-1*A_ELEC[n,e]*JuMP.value.(Flows_Elec[i,T,t,e]) for e = 1:EDGES_ELEC) - sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(JuMP.value.(charging_ELEC[i,T,t,s])-JuMP.value.(discharging_ELEC[i,T,t,s])) for s = 1:STORAGE_ELEC) - Demand_ELEC[i,T,t,n] - sum(P2G_NodalLoc_ELEC[n,d]*JuMP.value.(P2G_dispatch[i,T,t,d])*(1-ISBIOMETHANE[d]) for d = 1:P2G) for n = 1:NODES_ELEC)  for t = 1:t_ops) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(curtailmentRE[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
# end

# HourlyGenFullSave = zeros(T_inv*8760,GEN+STORAGE_ELEC+P2G+2)
# HourlyLoadFullSave = zeros(T_inv*8760,4)
# HourlySOCFullSave = zeros(T_inv*8760,STORAGE_ELEC)
# HourlyTransmissionFullSave = zeros(T_inv*8760,EDGES_ELEC)
# DailyGasSOCFullSave = zeros(T_inv*Periods_Per_Year,STORAGE_GAS)
# DailyGasTransmissionFullSave = zeros(T_inv*Periods_Per_Year,EDGES_GAS)
# for i = 1:T_inv
#     for c = 1:Periods_Per_Year
#         j = Int(RepDays[i,c])
#         count = Int((i-1)*8760+(c-1)*t_ops)+1
#         HourlyGenFullSave[count:count+t_ops-1,1:GEN] = JuMP.value.(generation[i,j,:,:])
#         HourlyGenFullSave[count:count+t_ops-1, GEN+1:GEN+STORAGE_ELEC] = (JuMP.value.(charging_ELEC[i,j,:,:])-JuMP.value.(discharging_ELEC[i,j,:,:]))
#         HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = JuMP.value.(P2G_dispatch[i,j,:,:].*transpose(ones(P2G)-ISBIOMETHANE))
#         HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+P2G+1] = sum(JuMP.value.(curtailmentRE[i,j,:,:]), dims = 2)
#         HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+P2G+2] = sum(JuMP.value.(LostLoad_ELEC[i,j,:,:]), dims = 2)
#         HourlyLoadFullSave[count:count+t_ops-1, 1] = sum(Demand_ELEC[i,j,:,:], dims = 2)
#         HourlyLoadFullSave[count:count+t_ops-1, 2] = sum(Demand_GAS[i,j,:,:], dims = 2)
#         HourlyLoadFullSave[count:count+t_ops-1, 3] = sum(BaselineDemand_ELEC[i,j,:,:], dims = 2)
#         HourlyLoadFullSave[count:count+t_ops-1, 4] = sum(BaselineDemand_GAS[i,j,:,:], dims = 2)
#         HourlySOCFullSave[count:count+t_ops-1,:] = JuMP.value.(storedEnergy_ELEC[i,j,:,:])
#         HourlyTransmissionFullSave[count:count+t_ops-1,:] = JuMP.value.(Flows_Elec[i,j,:,:])
#         DailyGasSOCFullSave[Int((i-1)*Periods_Per_Year*2 + c*2)-1:Int((i-1)*Periods_Per_Year*2 + c*2),:] = JuMP.value.(storedEnergy_GAS[i,j,:,:])
#         DailyGasTransmissionFullSave[Int((i-1)*Periods_Per_Year + c),:] = JuMP.value.(Flows_Gas[i,j,:])
#     end
# end

# CSV.write("FinalSim$(EITrajectory)GasQuality$(GasQuality)$(N_Periods)$(clustering_case)$(networkscen)TimeExtEmissionsandCostResults_$(case).csv",DataFrame(EmissionsAndCosts'), writeheader = true)
# CSV.write("FinalSim$(EITrajectory)GasQuality$(GasQuality)$(N_Periods)$(clustering_case)$(networkscen)TimeExtGenerationResults_$(case).csv",DataFrame(GenerationSave'), writeheader = true)
# CSV.write("FinalSim$(EITrajectory)GasQuality$(GasQuality)$(N_Periods)$(clustering_case)$(networkscen)TimeExtHourlyGensFullResults_$(case).csv",DataFrame(HourlyGenFullSave'), writeheader = true)
# CSV.write("FinalSim$(EITrajectory)GasQuality$(GasQuality)$(N_Periods)$(clustering_case)$(networkscen)TimeExtHourlyLoadFullBaseline_$(case).csv",DataFrame(HourlyLoadFullSave'), writeheader = true)

# print("Success!")


