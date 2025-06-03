#!/usr/bin/env julia

using Plots
using Printf
using DataStructures: OrderedDict
using NCDatasets
# using Arrow, DataFrames
using CSV, DataFrames
using Colors
using ColorSchemes
using Plots
using Printf
using LaTeXStrings
using Profile
using Statistics 

include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/calculate_physical_variables.jl") 
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/advance_variables.jl")
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/phytoplankton.jl")
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/forcings.jl") 
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/output.jl")
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/define_params.jl")


function run_forward_model(file_out_name::String, adjoint_ds::String)

    #********************** SPATIAL DOMAIN  ***************************
    N = global_params["N"]   # number of grid points
    H = global_params["H"]   # depth (meters)
    dz = global_params["dz"] # grid spacing - may need to adjust to reduce oscillations
    dt = global_params["dt"] # (seconds) size of time step
    M  = global_params["M"]  # number of time steps
    time_range = global_params["time_range"] # number of time steps

    file_out_name = "$(file_out_name)_$(time_range).nc"
    println("\n\nRunning the FORWARD PHYTOPLANKTON MODEL --> we are going forward in time")
    println("Adjusting our growth guess using the gamma from: $(adjoint_ds)")
    println("Will be saving phytoplankton output to: $(file_out_name)")


    ds = NCDataset("/pscratch/sd/s/siennaw/adjoint_phytoplankton/run_hydro/HYDRO_$time_range.nc")

    if adjoint_ds == "FIRST"
        calculate_gamma = true 
        println("First run: calculating gamma")
                
        # Read in the CIMIS data
        cimis_fn = "/pscratch/sd/s/siennaw/stockton_field_data/forcing_for_model/PAR_on_$time_range.csv"
        df = CSV.read(cimis_fn, DataFrame)
        par = df[!,"Sol Rad (PAR)"]
        println("Read in CIMIS data ...")
        function get_light(index::Int, par=par)
            return par[index]
        end
    else   
        calculate_gamma = false
        println("Using growth rate from adjoint model")
        gamma_ds = NCDataset("backward_lambda/$(adjoint_ds)_$time_range.nc")  #"../backward_lambda/adjoint_2.nc")
    end 
    #***********************************************************************






    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 1 
    var2save = ["algae1", "gamma"]

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

    # (4) Light 
    DIURNAL_LIGHT = false  
    background_turbidity =  3
    I_in = 350 

    #********************** DEFINE PHYTOPLANKTON FORCINGS ***************************
    init_algae = 0.005

    algae1 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
                "pmax" => 0.04 * hr2s,           # maximum specific growth rate [1/hour]
                "ws" => 1.38e-4, #1.38e-4,           # vertical velocity [m/s]
                "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
                "Li" => 0.005 * hr2s,             # specific loss rate [1/hour]
                "name" => "HAB",           # name of the species
                "self_shading" => true)    # self-shading effect (true/false)
    # '''
    # Diatoms ws = -1.38e-5 m/s
    # Cyanobacteria = 1.38e-4 m/s
    # '''

    #***************************************************************************
    #   Initialize variables
    #***************************************************************************

    # Create dictionary to hold important discretization parameters
    discretization = Dict("beta" => (dt/dz^2), "dz" => dz, "dt" => dt, "N" => N, "z"=> z, "H" => H)

    algae1["c"] = zeros(N) .+ init_algae 
    # algae2["c"] = zeros(N) .+ 1e-3 #@init_algae 

    Times = collect(1:dt:(M*dt))

    real_times_saved = []

    #***************************************************************************
    save2output(1, 1, "algae1", algae1["c"])
    # save2output(1, 1, "algae2", algae2["c"])

    variables = Dict("Kz" => ds["Kz"][:,1])

    for i in 2:M

        time = Times[i];
        # Hydrodynamics
        variables["Kz"] = ds["Kz"][:,i]

        # Phytoplankton
        if calculate_gamma 
            I0 = get_light(i)
            # light = self_shading(algae1, algae2, I0, background_turbidity, discretization)
            light = light_decay(I0, background_turbidity, discretization)
            gamma = calculate_net_growth(algae1, light, discretization)
        else 
            gamma = gamma_ds["gamma"][:,i]
        end 

        # Algae 
        algae1["c"] = advance_algae(variables, algae1, gamma, discretization)  # zeros(N) .+ init_algae  #
        save2output(time, i, "gamma", gamma)
        save2output(time, i, "algae1", algae1["c"])

        if algae1["c"][1] > 1
            println("Time: $(time) \t algae1: $(algae1["c"][1]) \t gamma: $(gamma[1])")
        end

    end



    # ********************** save data ****************************
    units_dict = Dict("U" => "m/s", 
        "C" => "deg C", 
        "Kz" => "m\$^2\$ s\$^{-1}\$", 
        "algae1" => L"10$^6$/cm$^3$ cells",
        "algae2" => L"10$^6$/cm$^3$ cells",
        "L" => "Turbulent length scale", 
        "Q2" => "TKE", "Q2L" => "TKE*L",
        "N_BV2" => "Brunt-Vaisala frequency", 
        "Kq" => "Kq", 
        "Nu" => "Nu_t",
        "gamma" => "Net growth rate",)

    var2name = Dict("U" => "Velocity", 
                "C" => "Temperature", 
                "Kz" => "Turbulent diffusivity", 
                "algae1" => "Diatom concentration",
                "algae2" => "HAB concentration",
                "L" => "Turbulent length scale", 
                "Q2" => "TKE",
                "Q2L" => "TKE*L",
                "N_BV2" => "Brunt-Vaisala frequency", 
                "Kq" => "Kq", 
                "Nu" => "Nu_t", 
                "gamma" => "Net growth rate [1/s]")

    times_unique = unique(Times) 
    # println("start + end of times unique $(times_unique[1]) $(times_unique[end])")
    # println("Times unique has $(length(times_unique)) elements \n")


    fout = "/pscratch/sd/s/siennaw/adjoint_phytoplankton/forward_phyto/$(file_out_name)"
    # fout = file_out_name
    # fout = "forward_phyto/$(file_out_name)"
    ds = NCDataset(fout,"c")
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

    # println("Size of saved gamma: $(size(ds["gamma"][:,:]))")

    print("Saved $file_out_name \n")
    close(ds)
    
end 



# file_out_name = "phyto_fake_truth"  
# run_forward_model(file_out_name, "FIRST")

# @profilehtml run_my_model(ws1, ws2, pmax1, pmax2, file_out_name)

# StatProfilerHTML.view()
# Profile.print() 
