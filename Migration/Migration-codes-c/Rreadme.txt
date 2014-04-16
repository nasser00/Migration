This folder contains one way wave equation shot profile migration algorithms under c programming for :

1) Phase Inversion Plus Interpolation (PSPI) 
Gazdag, J. and Sguazzero, P. (1984). �Migration of seismic data by phase shift plus interpolation.� GEOPHYSICS, 49(2), 124�131.
doi: 10.1190/1.1441643

2) Split?step Fourier migration (Split_Step)


Stoffa, P., Fokkema, J., de Luna Freire, R., and Kessinger, W. (1990). �Split?step Fourier migration.� GEOPHYSICS, 55(4), 410�421.
doi: 10.1190/1.1442850


I included post stack and pre-stack migration algorithms for Split_Step but only Post-stack migration in the case of PSPI.

The folder contains these sub_folders:

1) codes: include all codes that you need to compile (Note: you need Madagascar package ,which is a free software, to compile the two main codes that starts with big M).
  You can download madagascar from this webpage: http://www.ahay.org/wiki/Download
  
  
  
  Note: I will add part 2) and 3) later with a proper link to download the necessary files.
  2) post_stack_migration : includes SConstruct file and all the necessary input files to run the codes for post stack migration.
 After compiling the codes you can cd to this directory and just type scons view in your command prompt.
  
   3) post_stack_migration : includes SConstruct file and all the necessary input files to run the codes for prestack migration.
  after compiling the codes you can cd to this directory and just type scons view in your command prompt.
  
  If you have any questions about the package please send an email to:
  
  kazemino@ualberta.ca
  
  Nasser Kazemi
