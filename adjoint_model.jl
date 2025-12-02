#!/usr/bin/env julia

include("forward_phyto/run_forward_model.jl") 
include("backward_lambda/run_backward_model.jl")


#0.5,  5, 5, 5,  5, 5, 5,  5, 5,  5, 5,  5, 5, 5,  
step_sizes = [ 20,  10, 10, 5 , 5,5, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 ,  1, 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1 , 1,0.1,0.1, 0.1, 0.1, 0.1,  0.01,  0.01, 0.01,  0.01, 0.01,  0.01, 0.01, 0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01, 0.01, 0.01,0.001, 0.001] # Different step sizes to try
# step_sizes = [10, 5 , 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1,  0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001] # Different step sizes to try

# step_sizes = [ 20,10, 10,10, 10,10, 10,10, 10,10, 10,10, 10,10, 10,  10, 10, 5 , 5,5, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1, 1 , 1,0.1,0.1, 0.1, 0.1, 0.1,  0.01,  0.01, 0.01,  0.01, 0.01,  0.01, 0.01, 0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01,0.01, 0.01, 0.01, 0.01,0.001, 0.001] # Different step sizes to try

# adjoint_39_august_13.nc 1.6584314797922388e10
# adjoint_181_august_13.nc 1.6237770128660355e10


# First run where gamma is calculated based on physical mechanisms 
run_forward_model("forward_mc_1", "FIRST")
run_backward_model("adjoint_mc_1", "forward_mc_1", step_sizes[1])


# # Iterate! 
for i in 2:length(step_sizes)
    println("Running forward model iteration: $i")
    run_forward_model("forward_mc_$i", "adjoint_mc_$(i-1)")
    run_backward_model("adjoint_mc_$i", "forward_mc_$i", step_sizes[i])
    
end

# i = 205 
# while true
#     println("Running forward model iteration: $i")
#     run_forward_model("forward_$i", "adjoint_$(i-1)")
#     run_backward_model("adjoint_$i", "forward_$i", 3)
#     global i += 1
# end

# run_forward_model("forward_2.nc", "adjoint_1.nc")
# run_backward_model("adjoint_2.nc", "forward_2.nc")

# run_forward_model("forward_3.nc", "adjoint_2.nc")
# run_backward_model("adjoint_3.nc", "forward_3.nc")

# run_forward_model("forward_4.nc", "adjoint_3.nc")
# run_backward_model("adjoint_4.nc", "forward_4.nc")

# # Run forward model 
# forward_model_name = "forward_phyto_1.nc" 
# adjoint_model_name = "adjoint_2.nc"
# run_forward_model(forward_model_name, adjoint_model_name)

# # Run backward model
# backward_model_name = "adjoint_3.nc"
# run_backward_model(backward_model_name, forward_model_name)

# forward_model_name = "forward_phyto_2.nc"
# run_forward_model(forward_model_name, backward_model_name)

# backward_model_name = "adjoint_4.nc"
# run_backward_model(backward_model_name, forward_model_name)

# Run backward model