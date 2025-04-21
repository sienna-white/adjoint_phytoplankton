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

include("../run_hydro/calculate_physical_variables.jl") 
include("../run_hydro/advance_variables.jl")
include("../run_hydro/phytoplankton.jl")
include("../run_hydro/forcings.jl") 
include("../forward_phyto/output.jl")


# file_out_name = "adjoint_2.nc" 


using Random
Random.seed!(1234);      # Seed number 1234



function run_backward_model(file_out_name::String, algae_guess_ds:: String)

    println("\n\nRunning the ADJOINT OPERATOR MODEL --> we are going backward in time")
    println("Using the algae guess dataset: $(algae_guess_ds)")
    println("Will be saving adjoint variable output to: $(file_out_name)")


    ds_hydro = NCDataset("run_hydro/HYDRO.nc")
    ds_algae = NCDataset("forward_phyto/$(algae_guess_ds)")  #"../forward_phyto/phyto_GUESS.nc")
    ds_truth = NCDataset("forward_phyto/phyto_TRUTH.nc")
    
    println("Size of ds_algae gamma: $(size(ds_algae["gamma"][:,:]))")


    # println(ds_hydro) 
    # println(ds_algae)



    #********************** SPATIAL DOMAIN  ***************************
    N = 60    # number of grid points
    H = 6    # depth (meters)
    dz = H/N  # grid spacing - may need to adjust to reduce oscillations
    dt = 10   # (seconds) size of time step 
    M  = 10000 #00 #000 # 50000  #500 #

    # Increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
    isave = 1 #1000
    var2save = ["lambda"]

    create_output_dict(M, isave, var2save, N)

    # Create depth vector 
    z = collect(H:-dz:dz) .- dz/2

    # Swimming speed 
    ws = 1.38e-4

    lambda = zeros(N) 

    println("Initializing with ground truth + some noise")
    ground_truth_w_noise =  ds_truth["algae1"][:,end] #+ rand(N).*1e-3
    c_diff = 2*(ds_algae["algae1"][:,end] - ground_truth_w_noise)
    # c_diff = reverse(c_diff)
    # println(c_diff)

    # println(ds_truth['z'][:])

    # assert(false)

    L_n = c_diff

    # print(c_diff)


    Times = collect(1:dt:(M*dt))
    # Times = Times[1:end]

    save2output(Times[end], (M), "lambda", L_n) 


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
                bL[i] = 1 + ws*dt/dz + gamma[i]*dt + (dt/dz^2)*(1/2)*(kz[i+1] + 2*kz[i] + kz[i-1])
                cL[i] = - (dt/dz^2)*(1/2) * (kz[i] + kz[i+1])
                dL[i] = L_n[i] 
            end
        end 

        # Bottom-Boundary: no flux for scalars
        bL[1] =  1 + ws*dt/dz + (gamma[1]*dt) + (dt/dz^2)*(1/2)*(kz[1] + kz[2]) 
        cL[1] =  -ws*dt/dz - (dt/dz^2) * (1/2) * (kz[1] + kz[2])
        dL[1] =  L_n[1]

        # Top-Boundary: no flux for scalars
        aL[end] = - (dt/dz^2)*(1/2)* (kz[end] + kz[end-1])
        bL[end] = ws*dt/dz + 1 + gamma[end]*dt + (dt/dz^2)*(1/2)*(kz[end] + kz[end-1])  
        dL[end] = L_n[end]


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

        L_nminus1 = TDMA(aL, bL, cL, dL, N) 

        save2output(time, i, "lambda", L_nminus1)
        L_n = L_nminus1
    end 


#     # ********************** save data ****************************
    units_dict = Dict("lambda" => "[-]")
    var2name = Dict("lambda" => "Lagrangian multiplier")

    times_unique = unique(times) 
    println("length(times) = $(length(Times))")
    println("start + end of times unique $(times_unique[1]) $(times_unique[end-3:end])")
    println("Times unique has $(length(times_unique)) elements \n")


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


    println("Size of ds_algae gamma: $(size(ds_algae["gamma"][:,:]))")
    println("Size of output lambda: $(size(output["lambda"]))")

    grad = ds_algae["gamma"][:,:] .* output["lambda"]
    eps = 1

    new_gamma = @. ds_algae["gamma"][:,:] + grad*eps 

    v = defVar(ds, "gamma", Float64,("z","time"), attrib = OrderedDict(
        "units" =>  "-", "long_name" => "gradient descent parameterized growth"))
    v[:,:] = new_gamma;

    print("Saved $file_out_name \n")
    close(ds)



end 

# run_my_model(file_out_name)
