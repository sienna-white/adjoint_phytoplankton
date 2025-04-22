#!/usr/bin/env julia

# using Plots
using Printf
using DataStructures: OrderedDict
using NCDatasets
using Arrow, DataFrames
using CSV, DataFrames
using Colors
using Plots
using Printf
using LaTeXStrings
# using Profile
using Statistics 
using ColorSchemes

include("../model_code/calculate_physical_variables.jl") 
include("../model_code/advance_variables.jl")
include("../model_code/phytoplankton.jl")
include("../model_code/forcings.jl") 
include("../model_code/output.jl")

file_out_name = @sprintf("HYDRO.nc") 

function run_my_model(file_out_name::String)


    #***********************************************************************
    # Read in the feather file 
    data  = Arrow.Table("../model_code/interpolated_temperature_profile_aug10-16.feather")
    temp_data = DataFrame(data)

    function get_temp_field(index::Int, temp_data=temp_data)
        return collect(temp_data[index,1:end-1])
    end

    function get_unstrat_temp_field(index::Int, temp_data=temp_data)
        return fill( mean(collect(temp_data[index,1:end-1])), 60) 
    end

    #***********************************************************************
    # Read in the CIMIS data
    cimis_fn = "/global/homes/s/siennaw/scratch/siennaw/turbulence-model/data/CIMIS/PAR_on_august_9-15.csv"
    df = CSV.read(cimis_fn, DataFrame)
    par = df[!,"Sol Rad (PAR)"]
    println("Read in CIMIS data ...")

    function get_light(index::Int, par=par)
        return par[index]
    end
    #***********************************************************************


    #***********************************************************************
    # Wind time series 
    wind_fn = "/global/homes/s/siennaw/scratch/siennaw/turbulence-model/data/forcing_data/wind_on_august_10-16.csv"
    df = CSV.read(wind_fn, DataFrame)
    wind = df[!,"WindSpeed"]
    real_time = df[!,"time"]
    println("Read in wind data ...")
    
    
    function get_wind_speed(index::Int, wind=wind)
        return wind[index]
    end
    #***********************************************************************


    #********************** SPATIAL DOMAIN  ***************************
    N = 60    # number of grid points
    H = 6    # depth (meters)
    dz = H/N  # grid spacing - may need to adjust to reduce oscillations
    dt = 10   # (seconds) size of time step 
    M  = 10000 #00 #000 # 50000  #500 #

    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 1 #1000
    var2save = ["U", "C", "Kz", "L", "Q2","N_BV2"]

    create_output_dict(M, isave, var2save, N)

    # Create depth vector 
    z = collect(H:-dz:dz) .- dz/2 # depth vector
    # println("Length of z is ", length(z))

    #********************** FIXED CONSTANTS  ***************************
    rhoA = 1.23                     # Density of air, kg/m^3
    rhoW = 1000                     # Density of water, kg/m^3
    specific_heat_water = 4181      # J/kg-degC
    specific_heat_air = 1007        # J/kg-degC x RH
    c_d = 0.05                      # Drag coefficient [-]
    cm2m = 0.01
    hr2s = 1/3600

    #********************** INITIAL CONDITION ***************************
    # Initialize thermocline based on tanh curve 
    base_temp = 22
    dtemp = 1.5 
    stretch = 0.25 

    #********************** DEFINE HYDRODYNAMIC FORCINGS ***************************
    # (1) PRESSURE 
    Px0 = 2e-6          # Pressure gradient forcing
    T_Px = 12           # Period [hours] on pressure gradient forcing. Set to 0 for steady

    # (2) Wind
    # Wind = 1                       # u_star =m/s >> 0.05 is  drag coefficient, 10 is my wind speed 
    # WIND = (c_d * Wind)^2 * rhoA   # this is rho * u*^2

    # (3) Temperature
    top_temp = 33
    bottom_temp = 30
    bottom_speed = 0 
    top_speed=3.5  

    # (4) Light 
    DIURNAL_LIGHT = false  
    background_turbidity =  0.6
    I_in = 350 

    #***************************************************************************
    #   Initialize variables
    #***************************************************************************

    # Create dictionary to hold important discretization parameters
    discretization = Dict("beta" => (dt/dz^2), "dz" => dz, "dt" => dt, "N" => N, "z"=> z, "H" => H)

    # Initalize velocity
    U = similar(z) .+ 1e-1
    C = get_unstrat_temp_field(1) # zeros(N) .+ LinRange(bottom_temp, top_temp, N)  
    rho = calculate_rho(C, base_temp) 
    N_BV2 = calculate_brunt_vaisala(rho, discretization)

    Q2, Q2L, Q, L, Gh, nu_t, Kq, Kz = initialize_turbulent_functions(discretization, N_BV2)

    # Initial dictionary to store variables
    variables = Dict("U" => U, "C" => C, "N_BV2" => N_BV2, 
                    "Nu" => nu_t, "Q2" => Q2, "Q2L" => Q2L, 
                    "Kq" => Kq, "Kz" => Kz, "L" => L)

                    
    Times = collect(1:dt:(M*dt))
    # println("Times = ", Times)
    real_times_saved = []

    #***************************************************************************
    

    save2output(1, 1, "U", variables["U"])
    save2output(1, 1, "Kz", variables["Kz"])
    save2output(1, 1, "C", variables["C"])
    save2output(1, 1, "L", variables["L"])
    save2output(1, 1, "Q2", variables["Q2"])
    save2output(1, 1, "N_BV2", variables["N_BV2"])
    # push!(real_times_saved, real_time[1])

    for i in 2:M

        time = Times[i];

        # [1] Advance velocity field
        pressure = get_pressure_at_timestamp(time, Px0, T_Px)
        ustar = calculate_ustar(U)

        W0 = get_wind_speed(i)
        I0 = get_light(i)

        C = get_temp_field(i) #get_temp_field(i)              # [1] Observational, sttratified temperature field 
        rho = calculate_rho(C, base_temp)  # [2] Calculate density from temperature field  
        N_BV2 = calculate_brunt_vaisala(rho, discretization) # [3] Calculate Brunt-Vaisala frequency 
        N_BV2 = clamp.(N_BV2, -1e-3, Inf)               # [4] Prevent any unstable stratification during daylight hours
               
        # Advance velocity field 
        wind_stress = wind_speed_2_wind_stress(W0, discretization) 
        U = advance_velocity(variables, pressure, discretization, wind_stress)

        #  [2] Advance TKE / Q2 
        Q2 = advance_Q2(variables, ustar, discretization) 
        Q = @. sqrt(Q2)

        #  [3] Advance Q2*L      
        Q2L = advance_Q2L(variables, ustar, discretization)
    
        #  [4] Advance temperature 
        # C = advance_scalar(variables, discretization) 

        # [7] Semi-implicit: Calculate turbulent lengthscale
        L, Q2L = calculate_lengthscale(Q2, Q2L, N_BV2, discretization)

        # Calculate stability parameter 
        gh = calculate_Gh(N_BV2, L, Q)
        nu_t, Kq, Kz = calculate_turbulent_functions(gh, Q, L, discretization) 

        # [9] Pack variables for next timestep 
        variables["U"] = U
        variables["C"] = C
        variables["N_BV2"] = N_BV2
        variables["Nu"] = nu_t
        variables["Q2"] = Q2
        variables["Q2L"] = Q2L
        variables["Kq"] = Kq
        variables["Kz"] = Kz
        variables["L"] = L

        save2output(time, i, "U", variables["U"])
        save2output(time, i, "Kz", variables["Kz"])
        save2output(time, i, "C", variables["C"])
        save2output(time, i, "L", variables["L"])
        save2output(time, i, "Q2", variables["Q2"])
        save2output(time, i, "N_BV2", variables["N_BV2"])
        # push!(real_times_saved, real_time[i])

    end
   



    # ********************** save data ****************************
    units_dict = Dict("U" => "m/s", 
        "C" => "deg C", 
        "Kz" => "m\$^2\$ s\$^{-1}\$", 
        "algae1" => L"10$^6$/cm$^3$ cells",
        "algae2" => L"10$^6$/cm$^3$ cells",
        "L" => "Turbulent length scale", 
        "Q2" => "TKE", "Q2L" => "TKE*L",
        "N_BV2" => "Brunt-Vaisala frequency", "Kq" => "Kq", "Nu" => "Nu_t")

    var2name = Dict("U" => "Velocity", 
                "C" => "Temperature", 
                "Kz" => "Turbulent diffusivity", 
                "algae1" => "Diatom concentration",
                "algae2" => "HAB concentration",
                "L" => "Turbulent length scale", 
                "Q2" => "TKE","Q2L" => "TKE*L",
                "N_BV2" => "Brunt-Vaisala frequency", "Kq" => "Kq", "Nu" => "Nu_t")

    times_unique = unique(times) 

    
    # Print shape of 
    # println("Length of times_unique = ", length(times_unique))
    println("start + end of times unique $(times_unique[1]) $(times_unique[end])")
    println("Times unique has $(length(times_unique)) elements \n")

    ds = NCDataset(file_out_name,"c")

    # model_time = collect(1:M)
    defDim(ds, "z", length(z)) 
    defDim(ds, "time", length(times_unique))

    v = defVar(ds, "z", Float32, ("z",))
    v[:] = z

    v = defVar(ds, "time", Float32, ("time",), attrib = OrderedDict("units" => "seconds"))
    v[:] = collect(1:(length(times_unique))) #model_time

    for var in var2save
        # println(var)
        v = defVar(ds, var, Float64,("z","time"), attrib = OrderedDict(
        "units" =>  units_dict[var], "long_name" => var2name[var]))
        v[:,:] = output[var];
    end

    print("Saved $file_out_name \n")
    close(ds)
    
end 

run_my_model(file_out_name)


