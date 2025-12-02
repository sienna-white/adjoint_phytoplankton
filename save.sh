niter=95


cd forward_phyto/
mv forward_mc_${niter}_* ../finished_mc/
mv forward_mc_1_* ../finished_mc/

cd ../backward_lambda/
mv adjoint_mc_${niter}_* ../finished_mc/
mv adjoint_mc_1_* ../finished_mc/

# cd forward_phyto/
# mv forward_${niter}_* ../finished/
# mv forward_1_* ../finished/

# cd ../backward_lambda/
# mv adjoint_${niter}_* ../finished/
# mv adjoint_1_* ../finished/

echo "done"