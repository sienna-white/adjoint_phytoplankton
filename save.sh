niter=297


cd forward_phyto/
mv forward_${niter}_* ../finished2/
mv forward_1_* ../finished2/

cd ../backward_lambda/
mv adjoint_${niter}_* ../finished2/
mv adjoint_1_* ../finished2/

# cd forward_phyto/
# mv forward_${niter}_* ../finished/
# mv forward_1_* ../finished/

# cd ../backward_lambda/
# mv adjoint_${niter}_* ../finished/
# mv adjoint_1_* ../finished/

echo "done"