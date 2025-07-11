#!/usr/bin/env julia

using Printf
using DataStructures: OrderedDict
using NCDatasets
using DataFrames
using CSV, DataFrames
using Colors
using ColorSchemes
using Plots
using Printf
using LaTeXStrings
using LinearAlgebra
# using Profile
using Statistics 

include("/global/homes/s/siennaw/scratch/siennaw/two_species/adjoint_phytoplankton/model_code/calculate_physical_variables.jl") 
include("/global/homes/s/siennaw/scratch/siennaw/two_species/adjoint_phytoplankton/model_code/advance_variables.jl")
include("/global/homes/s/siennaw/scratch/siennaw/two_species/adjoint_phytoplankton/model_code/phytoplankton.jl")
include("/global/homes/s/siennaw/scratch/siennaw/two_species/adjoint_phytoplankton/model_code/forcings.jl") 
include("/global/homes/s/siennaw/scratch/siennaw/two_species/adjoint_phytoplankton/model_code/output.jl")


using Random
Random.seed!(1234);      # Seed number 1234

function run_backward_model(file_out_name::String, algae_guess_ds:: String)

    println("\n\nRunning the ADJOINT OPERATOR MODEL --> we are going backward in time")
    println("Using the algae guess dataset: $(algae_guess_ds)")
    println("Will be saving adjoint variable output to: $(file_out_name)")

    #********************** SPATIAL DOMAIN  ***************************
    N = global_params["N"]   # number of grid points
    H = global_params["H"]   # depth (meters)
    dz = global_params["dz"] # grid spacing - may need to adjust to reduce oscillations
    dt = global_params["dt"] # (seconds) size of time step
    M  = global_params["M"]  # number of time steps
    time_range = global_params["time_range"] # number of time steps

    file_out_name = "$(file_out_name)_$(time_range).nc"

    ds_hydro = NCDataset("/pscratch/sd/s/siennaw/adjoint_phytoplankton/run_hydro/HYDRO_$time_range.nc")
    ds_algae = NCDataset("/pscratch/sd/s/siennaw/two_species/adjoint_phytoplankton/forward_phyto/$(algae_guess_ds)_$(time_range).nc")  #"../forward_phyto/phyto_GUESS.nc")

    # INITIALIZE THE ADJOINT FORCING --> DIFF BETWEEN MODEL & OBS 
    # println("Initializing with ground truth + some noise")
    ds_truth = NCDataset("/pscratch/sd/s/siennaw/adjoint_phytoplankton/forward_phyto/phyto_fake_truth_june22_august_13.nc")
    # ground_truth_w_noise =  ds_truth["algae1"][:,end] + rand(N).*1e-6
    # c_diff = 2*(ds_algae["algae1"][:,end] - ground_truth_w_noise)

    # # CREATE dictionary
    adj_forcing = Dict() 
    # adj_forcing[(M-1)] = c_diff

    # # Let's say we have observations at times 
    # for i in 1:500:M
    #     # measurement is at depth N = 20
    #     OBS_DEPTH = 20 
    #     forcing = zeros(N)
    #     forcing[OBS_DEPTH] =  2*(ds_algae["algae1"][OBS_DEPTH, i] - (ds_truth["algae1"][OBS_DEPTH, i] + rand()*1e-6))
    #     adj_forcing[i] = forcing
    # end

    use_penalty = true 
    PMAX = 2e-4 
    PMIN = -1e-5

    function calculate_penalty(gamma, use_penalty, PMAX, PMIN)
        penalty = zeros(N)
        mu = 0.01 * 10
        if use_penalty
            for i in 1:N
                if gamma[i] > PMAX
                    penalty[i] = mu * (gamma[i] - PMAX)
                elseif gamma[i] < PMIN
                    penalty[i] = mu * (gamma[i] - PMIN)
                end
            end
        end 
        
        return penalty
    end

    # CREATE dictionary
    adj_forcing = Dict() 

    df = CSV.read("/global/homes/s/siennaw/scratch/siennaw/stockton_field_data/chla/chla_august_13.csv", DataFrame)
    time_steps = df[!, "model_time"]
    chla_vals = df[!, "VALUE"]

    cost = 0 
    for i in 1:length(time_steps)  
        time_step = time_steps[i]
        time_step_int = time_step + 1 #parse(Int, time_step) # Convert to integer
        chla_val = chla_vals[i]   # Get profile data at that point 
        chla_val = chla_val .* 1e-6 
        total_modeled_algae = ds_algae["algae1"][55, time_step_int] + ds_algae["algae2"][55, time_step_int]
        difference_ = (total_modeled_algae - chla_val) 
        # println("Difference at time step $(time_step_int) is $(difference_)")
        difference = zeros(N)
        difference[50:56] .= difference_
        cost += difference_
        adj_forcing[time_step_int] = difference
    end 


    # # Read csv file 
    # df = CSV.read("/pscratch/sd/s/siennaw/stockton_field_data/profiler/profiles_cells_august_13_averaged.csv", DataFrame)
    
    # # Get list of columns
    # time_steps = names(df)
    # # println("Time steps in the DataFrame: $(time_steps)")
    # for i in 1:length(time_steps)
        
    #     time_step = time_steps[i]
    #     if time_step == "z" || time_step == "Column1" 
    #         continue
    #     end
    #     # println("Processing time step: $(time_step)")
    #     time_step_int = parse(Int, time_step) # Convert to integer
    #     profile = df[!, time_step]   # Get profile data at that point 
    #     profile = profile .* 1e-6 
    #     total_modeled_algae = ds_algae["algae1"][:, time_step_int] + ds_algae["algae2"][:, time_step_int]
    #     difference = (total_modeled_algae - profile) 
    #     # difference[1:20] .= 0 
        
    #     cost += sum(difference)
    #     adj_forcing[time_step_int] = difference * 0.25
    # end 


    # # Let's say we have observations at times 
    # cost = 0 
    # for i in 1:500:M
    #     # OBS_DEPTH = 20         # measurement is at depth N = 20
    #     forcing = zeros(N)
    #     forcing[10:end] = @. 2*(ds_algae["algae1"][10:end, i] - (ds_truth["algae1"][10:end, i]))#+ rand()*1e-6))
    #     cost += sum(forcing)
    #     adj_forcing[i] = forcing
    # end

    println("Total cost of adjoint forcing is $(cost)")
    open("cost.txt","a") do io
        println(io,"$file_out_name $cost")
    end

    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 1 #1000
    var2save = ["lambda1", "lambda2"]      # Only save growth + algae

    create_output_dict(M, isave, var2save, N)

    # Create depth vector 
    z = collect(H:-dz:dz) .- dz/2

    hr2s = 1/3600 
    algae1 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
                "pmax" => 0.005 * hr2s,           # maximum specific growth rate [1/hour]
                "ws" => 1.38e-4, #1.38e-4,           # vertical velocity [m/s]
                "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
                "Li" => 0.005 * hr2s,             # specific loss rate [1/hour]
                "name" => "HAB",           # name of the species
                "self_shading" => true)    # self-shading effect (true/false)

    algae2 = Dict("k" => 0.034,              # specific light attenuation coefficient [cm^2 / 10^6 cells]
            "pmax" => 0.05 * hr2s,           # maximum specific growth rate [1/hour]
            "ws" => -1.38e-5,           # vertical velocity [m/s]
            "Hi" => 40,                # half-saturation of light-limited growth [mu mol photons * m^2/s]
            "Li" => 0.005 * hr2s,             # specific loss rate [1/hour]
            "name" => "Diatoms",           # name of the species
            "self_shading" => true)    # self-shading effect (true/false)


    discretization = Dict("beta" => (dt/dz^2), "dz" => dz, "dt" => dt, "N" => N, "z"=> z, "H" => H)

    lambda1 = zeros(N) 
    lambda2 = zeros(N)

    L1_n = zeros(N) 
    L2_n = zeros(N)

    Times = collect(1:dt:(M*dt))
    save2output(Times[end], M, "lambda1", L1_n) 
    save2output(Times[end], M, "lambda2", L2_n)

    for i in (M-1):-1:1
        # println("Time = $i") 

        time = Times[i];
        kz   = ds_hydro["Kz"][:,i]

        cost = zeros(N)
        if haskey(adj_forcing, i)
            # println("adj has forcing @ time $(i)")
            # println("Sum of penalty for alg1 = $(sum(penalty1))")
            # println("Sum of penalty for alg2 = $(sum(penalty2))")
            # println("Sum of cost = $(sum(adj_forcing[i]))")
            cost = adj_forcing[i]
        end 

        # LAMBDA 1 (LANGRANGIAN MULTIPLIER FOR HABS /ALGAE 1)
        #  ... settling speed is UPWARD (swimming!)
        
        # aL, bL, cL, dL = initialize_abcd(N)
        # ws = algae1["ws"]  # vertical velocity [m/s]

        gamma1 = ds_algae["gamma1"][:,i] .- algae1["Li"]
        penalty1 = calculate_penalty(gamma1, use_penalty, PMAX, PMIN)

        lambda1 = advance_lagrangian_multiplier(algae1, L1_n, kz, gamma1, discretization, penalty1, cost)
        save2output(time, i, "lambda1", lambda1)
        L1_n = lambda1


        gamma2 = ds_algae["gamma2"][:,i] .- algae1["Li"]
        penalty2 = calculate_penalty(gamma2, use_penalty, PMAX, PMIN)
        lambda2 = advance_lagrangian_multiplier(algae2, L2_n, kz, gamma2, discretization, penalty2, cost)
        save2output(time, i, "lambda2", lambda2)
        L2_n = lambda2


        # if ws>0
        #     for i in 2:(N-1)
        #         aL[i] =  -ws*dt/dz - (dt/dz^2)*(1/2)*(kz[i-1] + kz[i])
        #         bL[i] = 1 + ws*dt/dz - gamma1[i]*dt + (dt/dz^2)*(1/2)*(kz[i+1] + 2*kz[i] + kz[i-1])
        #         cL[i] = - (dt/dz^2)*(1/2) * (kz[i] + kz[i+1])
        #         dL[i] = L1_n[i] + penalty[i]*dt 
        #     end
        # end 

        # # Bottom-Boundary: no flux for scalars
        # bL[1] =  1 + ws*dt/dz - (gamma[1]*dt) + (dt/dz^2)*(1/2)*(kz[1] + kz[2]) 
        # cL[1] =  -ws*dt/dz - (dt/dz^2) * (1/2) * (kz[1] + kz[2])
        # dL[1] =  L_n[1] + penalty[1]*dt

        # # Top-Boundary: no flux for scalars
        # aL[end] = - (dt/dz^2)*(1/2)* (kz[end] + kz[end-1])
        # bL[end] = ws*dt/dz + 1 - gamma[end]*dt + (dt/dz^2)*(1/2)*(kz[end] + kz[end-1])  
        # dL[end] = L_n[end] + penalty[end]*dt

        # # initial condition 
        # if haskey(adj_forcing, i)
        #     # println("adj has forcing @ time $(i)")
        #     println("Sum of penalty = $(sum(penalty))")
        #     println("Sum of cost = $(sum(adj_forcing[i]))")
        #     dL = dL .+ adj_forcing[i].*dt 
        # end 
   
        # L_nminus1 = TDMA(aL, bL, cL, dL, N) 
        # save2output(time, i, "lambda", L_nminus1)
        # L_n = L_nminus1
    end 

#     # ********************** save data ****************************
    units_dict = Dict("lambda1" => "[-]", "lambda2" => "[-]")
    var2name = Dict("lambda1" => "Lagrangian multiplier", "lambda2" => "Lagrangian multiplier")


    fout = "backward_lambda/$(file_out_name)"
    ds = NCDataset(fout,"c")
    nt = div(M,isave) + 1 
    defDim(ds, "z", length(z)) 
    defDim(ds, "time", nt)

    v = defVar(ds, "z", Float32, ("z",))
    v[:] = z

    v = defVar(ds, "time", Float32, ("time",), attrib = OrderedDict("units" => "seconds"))
    v[:] = collect(1:nt)

    for var in var2save
        v = defVar(ds, var, Float64,("z","time"), attrib = OrderedDict(
        "units" =>  units_dict[var], "long_name" => var2name[var]))
        v[:,:] = output[var];
    end


    # n = norm(grad)
    # println("Norm is $(n)")
    eps = 15 #2e4 * n 
    # [1] Adjust growth rates 
    grad = ds_algae["gamma1"][:,:] .* output["lambda1"]
    new_gamma = @. ds_algae["gamma1"][:,:] - grad*eps 
    v = defVar(ds, "gamma1", Float64,("z","time"), attrib = OrderedDict(
        "units" =>  "-", "long_name" => "gradient descent parameterized growth"))
    v[:,:] = new_gamma;



    # println("Norm of $grad is $(norm(grad))")

     # [2] Adjust growth rates 
    grad = ds_algae["gamma2"][:,:] .* output["lambda2"]
    new_gamma = @. ds_algae["gamma2"][:,:] - grad*eps 
    v = defVar(ds, "gamma2", Float64,("z","time"), attrib = OrderedDict(
        "units" =>  "-", "long_name" => "gradient descent parameterized growth"))
    v[:,:] = new_gamma;


    print("Saved $file_out_name \n")
    close(ds)

    # tdiff = abs(sum(sum(grad.*eps))) 
    # # println("Increment size is $(tdiff)") 
    # if tdiff < 1e-10 #1e-6 # changed from 2 
    #     println("HITTING BELOW THE THRESHOLD!!!!!")
    #     println("STOPPING...")
    #     # Stop the julia script
    #     exit(0)

        
    # end 


end 

# run_my_model(file_out_name)


 # else
    #     for i in 2:(N-1)
    #         aA[i]  = -wsdtdz - beta/2 * (Kz_past[i-1]+ Kz_past[i]) 
    #         bA[i]  = 1 + wsdtdz - gamma[i]*dt  + beta/2*(Kz_past[i+1] + 2*Kz_past[i] + Kz_past[i-1]) 
    #         cA[i]  = -beta/2 * (Kz_past[i] + Kz_past[i+1])
    #         dA[i]  = A_past[i]
    #     end 
           
    #     # Bottom-Boundary: no flux for scalars
    #     bA[1] =  1 + wsdtdz - (gamma[1]*dt) + beta/2*(Kz_past[2] + Kz_past[1]) 
    #     cA[1] =  -wsdtdz -beta/2 * (Kz_past[2] + Kz_past[1])
    #     dA[1] =  A_past[1]

    #     # Top-Boundary: no flux for scalars
    #     aA[end] =  -beta/2 * (Kz_past[end] + Kz_past[end-1])
    #     bA[end] =  1 - gamma[end]*dt + beta/2 * (Kz_past[end] + Kz_past[end-1]) + wsdtdz # okay adding this here 
    #     dA[end] = A_past[end]   

        # if ws>0
        #     for i in 2:(N-1)
        #         aL[i] = ws*dt/dz - (dt/dz^2)*(1/2)*(kz[i-1] + kz[i])
        #         bL[i] = 1 - ws*dt/dz - gamma[i]*dt + (dt/dz^2)*(1/2)*(kz[i+1] + 2*kz[i] + kz[i-1])
        #         cL[i] = - (dt/dz^2)*(1/2) * (kz[i] + kz[i+1])
        #         dL[i] = L_n[i] 
        #     end
        # end 

        # # Bottom-Boundary: no flux for scalars
        # bL[1] =  1 - (gamma[1]*dt) + (dt/dz^2)*(1/2)*(kz[1] + kz[2]) 
        # cL[1] =  -(dt/dz^2) * (1/2) * (kz[1] + kz[2])
        # dL[1] =  L_n[1]

        # # Top-Boundary: no flux for scalars
        # aL[end] = ws*dt/dz - (dt/dz^2)*(1/2)* (kz[end] + kz[end-1])
        # bL[end] = 1 - gamma[end]*dt + (dt/dz^2)*(1/2)*(kz[end] + kz[end-1])  
        # dL[end] = L_n[end]