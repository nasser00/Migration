I included 7 matlab codes in this folder

1) Gazdag_adjoint
2) Gazdag_forward
3) PSPI_adjoint
4) PSPI_adjoint_parallel
5) Split Step_adjoint
6) Split_Step_adjoint_parallel
4) test

You need to run the test.m code to check for impulse response,
linear event migration and forward application on migrated linear data
in the case of GAZDAG migration.

You also need Kevin_taper.mat file to run the programs (3-6)
For parameters you can do like this:

M=PSPI_adj_new(data_taper,dt,(0:699)*10,0,5,599*5,vel,2,3,30,1,1);

You can find some preliminary results in Results.pdf file.