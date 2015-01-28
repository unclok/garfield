To use this garfield++ library, you need to set an environment variable GARFIELD_HOME="path_to_trunk_in_this_directory"

All of test source using garfield++, are using both root and garfield++ library. A makefile already is in the test source directory.

In garfield/test/gemsim, there is the code gasfile.C which used for comparison with geant4. 

The file named gem in the same directory is a script for ansys batchjob to calculate an electric field of gem. You can use ansys batch in the lxplus server of cern. To use ansys batch, you need to set the environment variable PATH=$PATH:/afs/cern.ch/project/parc/bin, excute anasbatch15 gem and you can get 4 files, ELIST.lis,  MPLIST.lis,  NLIST.lis,  PRNSOL.lis. These files can be used for garfiled++ simulation.

The file named gas in the same directory is a script for ansys batchjob to calculate an electric field of vacuum. Usage of this file is same with gem.
