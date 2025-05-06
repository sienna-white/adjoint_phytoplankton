# adjoint_phytoplankton


Things to test:

1. What's a reasonable stopping point (how close do multiple guesses for gamma have to be?)
2. Test model with surface time series


threshold = 2e-5

--> added surface time series every 500 time steps (500*10 seconds = every 1.4 hours)

starting guess = 0.06
    *  converged in 30 steps 

starting guess = 0.05 
    * converged in 28 steps
    * 

starting guess = 0.04 
    * doesn't converge ... 

starting guess = 0.1 
    * converges in 30 

starting guess = 0.1, loss goes from "Li" => 0.004 --> "Li" => 0.003
    * diff < 3.9 e-5 after 100+ iterations 
    * 3.71912687537951e-5 after 200 
    * WOOOO WE DID IT : adjoint_697 dropped below 2e-5
    * 697 iterations !! 

starting guess = 0.11