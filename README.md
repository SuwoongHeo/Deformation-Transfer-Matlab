# Deformation Transfer Matlab
The matlab implementation of Deformation transfer[1]. I had implemented this code by referring the c++ implementation code[2]. In order to run this code, please download the mesh data which can be found in project website[3]. If you have any question, please contact vereurer@gmail.com

- If you have problem with using parallel for loop(parfor), run below code in commandline
distcomp.feature( 'LocalUseMpiexec', false )

[1] : Sumner, Robert W., and Jovan Popović. "Deformation transfer for triangle meshes." ACM Transactions on Graphics (TOG) 23.3 (2004): 399-405.

[2] : https://github.com/Golevka/deformation-transfer

[3] : http://people.csail.mit.edu/sumner/research/deftransfer/
