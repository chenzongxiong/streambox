Changes from Kaskade7.2 to Kaskade7.3 which may affect Kaskade7 user applications:
----------------------------------------------------------------------------------

Changes which have been already semi-automatically applied to every .cpp and .hh file in the
Kaskade7.3 files tree:
------------------------------------------------------------------------------------------------

1. Replaced VariableSetDescription<..>::Representation with VariableSetDescription<..>::VariableSet
2. Changed the enums PrecondType, IterateType, DirectType and MatrixProperties in the 
   utilities/emums.hh file into enum classes

Further changes:
----------------------------------------------------------------------------------------------------
1. In the tutorial examples, we use the new Kaskade7 CG solver which has been contributed by Lars Lubkoll,
   instead of the Dune::CGSolver.
2. A new tutorial example "artificial_1d_testProblem" (atp) has been added, featuring the use of the 
   Dune::OneDGrid.
3. The contents of the experiments directory and moreExamples directory has been deleted with the exception 
   of the subdirectories cmg and geomgrid, which have been moved to the tutorial directory.
   Each one, who wishes, may create a subdirectory in the experiments directory, named like her/his surname,
   and use this directory as a private directory to store own Kaskade7 related stuff.
4. The benchmarking directory has been deleted, and the contents of this directory - with the exception of
   the hlrnii subdirectory - has been moved to the tests directory. The hlrnii subdirectory has been
   deleted.
   
