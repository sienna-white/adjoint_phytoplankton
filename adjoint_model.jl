#!/usr/bin/env julia



# Run forward model 

include("forward_phyto/run_forward_model.jl") 
include("backward_lambda/run_backward_model.jl")


# Run forward model 
forward_model_name = "forward_phyto_1.nc" 
adjoint_model_name = "adjoint_2.nc"
run_forward_model(forward_model_name, adjoint_model_name)

# Run backward model
backward_model_name = "adjoint_3.nc"
run_backward_model(backward_model_name, forward_model_name)

forward_model_name = "forward_phyto_2.nc"
run_forward_model(forward_model_name, backward_model_name)

backward_model_name = "adjoint_4.nc"
run_backward_model(backward_model_name, forward_model_name)

# Run backward model