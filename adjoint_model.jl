#!/usr/bin/env julia

include("forward_phyto/run_forward_model.jl") 
include("backward_lambda/run_backward_model.jl")

# First run where gamma is calculated based on physical mechanisms 
run_forward_model("forward_1", "FIRST")
run_backward_model("adjoint_1", "forward_1")


# Iterate! 
for i in 495:800
    println("Running forward model iteration: $i")
    run_forward_model("forward_$i", "adjoint_$(i-1)")
    run_backward_model("adjoint_$i", "forward_$i")
end


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