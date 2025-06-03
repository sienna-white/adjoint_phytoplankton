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
# using Profile
using Statistics 

include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/calculate_physical_variables.jl") 
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/advance_variables.jl")
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/phytoplankton.jl")
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/forcings.jl") 
include("/pscratch/sd/s/siennaw/adjoint_phytoplankton/model_code/output.jl")


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
    ds_algae = NCDataset("/pscratch/sd/s/siennaw/adjoint_phytoplankton/forward_phyto/$(algae_guess_ds)_$(time_range).nc")  #"../forward_phyto/phyto_GUESS.nc")

    # INITIALIZE THE ADJOINT FORCING --> DIFF BETWEEN MODEL & OBS 
    # println("Initializing with ground truth + some noise")
    ds_truth = NCDataset("/pscratch/sd/s/siennaw/adjoint_phytoplankton/forward_phyto/phyto_fake_truth_august_13.nc")
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

    # CREATE dictionary
    adj_forcing = Dict() 

    # # Read csv file 
    # df = CSV.read("/pscratch/sd/s/siennaw/stockton_field_data/profiler/profiles_cells_august_13_CLEANED.csv", DataFrame)
    
    # # Get list of columns
    # time_steps = names(df)
    # for i in 1:length(time_steps)
        
    #     time_step = time_steps[i]
    #     if time_step == "z"
    #         continue
    #     end
    #     time_step_int = parse(Int, time_step) # Convert to integer
    #     profile = df[!, time_step]   # Get profile data at that point 
    #     profile = profile .* 1e-6 
    #     difference = 2* (ds_algae["algae1"][:, time_step_int] - profile) 
    #     difference[1:20] .= 0 
        
    #     println("Found a profile at time step $(time_step)\n")
    #     # println("Profile at time step $(time_step) is $(profile)\n")
    #     println("Algae at time step $(time_step) is $(ds_algae["algae1"][1:3, time_step_int])\n")
    #     print("Difference is $(difference[end-2:end])\n")
    #     adj_forcing[time_step_int] = difference 
    # end 

    # print(adj_forcing)
    # println("Column names: ", col_names)
    # assert(false)
    # par = df[!,"Sol Rad (PAR)"]

    # # Let's say we have observations at times 
    for i in 1:500:M
        # OBS_DEPTH = 20         # measurement is at depth N = 20
        forcing = zeros(N)
        forcing[10:end] = @. 2*(ds_algae["algae1"][10:end, i] - (ds_truth["algae1"][10:end, i] ))#+ rand()*1e-6))
        adj_forcing[i] = forcing
    end

    # println("Size of ds_algae gamma: $(size(ds_algae["gamma"][:,:]))")

    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 1 #1000
    var2save = ["lambda"]

    create_output_dict(M, isave, var2save, N)

    # Create depth vector 
    z = collect(H:-dz:dz) .- dz/2

    # Swimming speed 
    ws = 1.38e-4

    lambda = zeros(N) 



    # println("C_diff = $(c_diff)")
    # c_diff = 2*(ground_truth_w_noise - ds_algae["algae1"][:,end])

    # c_diff = reverse(c_diff)
    # println(c_diff)

    # println(ds_truth['z'][:])

    # assert(false)

    L_n = zeros(N) 


    Times = collect(1:dt:(M*dt))
    save2output(Times[end], M, "lambda", L_n) 


    for i in (M-1):-1:1
        # println("Time = $i") 

        time = Times[i];

        #  If settling speed is UPWARD (swimming!)
        gamma = ds_algae["gamma"][:,i]
        kz    = ds_hydro["Kz"][:,i]
        aL, bL, cL, dL = initialize_abcd(N)
         
         if ws>0
            for i in 2:(N-1)
                aL[i] =  -ws*dt/dz - (dt/dz^2)*(1/2)*(kz[i-1] + kz[i])
                bL[i] = 1 + ws*dt/dz - gamma[i]*dt + (dt/dz^2)*(1/2)*(kz[i+1] + 2*kz[i] + kz[i-1])
                cL[i] = - (dt/dz^2)*(1/2) * (kz[i] + kz[i+1])
                dL[i] = L_n[i] 
            end
        end 

        # Bottom-Boundary: no flux for scalars
        bL[1] =  1 + ws*dt/dz - (gamma[1]*dt) + (dt/dz^2)*(1/2)*(kz[1] + kz[2]) 
        cL[1] =  -ws*dt/dz - (dt/dz^2) * (1/2) * (kz[1] + kz[2])
        dL[1] =  L_n[1]

        # Top-Boundary: no flux for scalars
        aL[end] = - (dt/dz^2)*(1/2)* (kz[end] + kz[end-1])
        bL[end] = ws*dt/dz + 1 - gamma[end]*dt + (dt/dz^2)*(1/2)*(kz[end] + kz[end-1])  
        dL[end] = L_n[end]

        # initial condition 
        if haskey(adj_forcing, i)
            # println("adj has forcing @ time $(i)")
            dL = dL .+ adj_forcing[i].*dt 
        end 

        # c_diff
        # if i == (M-1)
        #     dL = dL .+ c_diff.*dt 
        # end 
   
        L_nminus1 = TDMA(aL, bL, cL, dL, N) 

        save2output(time, i, "lambda", L_nminus1)
        L_n = L_nminus1
    end 


#     # ********************** save data ****************************
    units_dict = Dict("lambda" => "[-]")
    var2name = Dict("lambda" => "Lagrangian multiplier")

    times_unique = unique(times) 
    # println("length(times) = $(length(Times))")
    # println("start + end of times unique $(times_unique[1]) $(times_unique[end-3:end])")
    # println("Times unique has $(length(times_unique)) elements \n")


    fout = "backward_lambda/$(file_out_name)"
    ds = NCDataset(fout,"c")
    defDim(ds, "z", length(z)) 
    defDim(ds, "time", length(times_unique))

    v = defVar(ds, "z", Float32, ("z",))
    v[:] = z

    v = defVar(ds, "time", Float32, ("time",), attrib = OrderedDict("units" => "seconds"))
    v[:] = collect(1:(length(times_unique)))

    for var in var2save
        v = defVar(ds, var, Float64,("z","time"), attrib = OrderedDict(
        "units" =>  units_dict[var], "long_name" => var2name[var]))
        v[:,:] = output[var];
    end


    # println("Size of ds_algae gamma: $(size(ds_algae["gamma"][:,:]))")
    # println("Size of output lambda: $(size(output["lambda"]))")

    grad = ds_algae["gamma"][:,:] .* output["lambda"]
    println("grad: $(grad[1:5])")
    eps = 0.5

    new_gamma = @. ds_algae["gamma"][:,:] - grad*eps # plus or minus???


    v = defVar(ds, "gamma", Float64,("z","time"), attrib = OrderedDict(
        "units" =>  "-", "long_name" => "gradient descent parameterized growth"))
    v[:,:] = new_gamma;

    print("Saved $file_out_name \n")
    close(ds)

    tdiff = abs(sum(sum(grad.*eps))) 
    println("Increment size is $(tdiff)") 
    if tdiff < 1e-10 #1e-6 # changed from 2 
        println("HITTING BELOW THE THRESHOLD!!!!!")
        println("STOPPING...")
        # Stop the julia script
        exit(0)

        
    end 


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