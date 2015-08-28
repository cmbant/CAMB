===================
CAMB Test/Benchmark/Profiler
===================

This folder contains the CAMB Test/Benchmark/Profile (TBP) suite.

If you are running this for the first time just issue:

make create_legacy

to create the baseline results that will be used for comparison and:

make all

to run all the test, the benchmarks and the profile.

The CAMB TBP suite works with bash scripts and requires numdiff.

This code was developed by Marco Raveri (mraveri@sissa.it)

Description
=============================

Several make targets allow to:

1) make all : run the test comparison, the benchmarks and the profiling. It might take a while depending on the number of parameter files;

2) make (default) : the default make target is "all";

3) make create_spectra : creates the test spectra. Recompiles CAMB every time to ensure coherence;

4) make run_comparison : run the comparison of the test spectra with the legacy ones;

5) make run_benchmarks : run the benchmarks. Recompiles CAMB every time to ensure coherence;

6) make run_profile : run the profiler only. Recompiles CAMB every time to ensure coherence;

7) make create_legacy : run CAMB to create the baseline results against which the new ones will be compared;

8) make dir : create the directory structure;

9) make clean : clean all the results leaving legacy ones;

10) make clean_* : clean the results of a specific operation;

11) make deep_clean : remove everything and restore the TBP to the original state. Removes all results and legacy results.
