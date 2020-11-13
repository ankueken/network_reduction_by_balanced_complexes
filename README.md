# Network reduction by balanced complexes

Code for reproduction of data presented in “A structural property for reduction of biochemical networks”.

The structural property of balancing of complexes can be effectively used in model reduction under mild network structure conditions, while preserving the steady-states supported by the network.

System requirements: 
--------------------

Standard computer with any operating system, the code was tested 	on Windows 10.
To run the code R language with package igraph and R.matlab as well 	as MATLAB together with Optimization Toolbox, glpk solver and CobraToolbox 
Heirendt et al. 2019, Nature Protocols https://opencobra.github.io/cobratoolbox/stable/installation.html) for MATLAB have to be installed.

Tested for R version 4.0.2, igraph v. 1.2.5, R.matlab v. 3.6.2,
MATLAB R2019a, Optimization Toolbox Version 8.3, Cobra Toolbox 2017

Model structure requirements:
-----------------------------

Networks to analyze need to have at least the following fields (according to cobra model structure):
		.S – stoichiometric matrix (type double, size m x n),
		.mets – cell array of metabolite names (size m x 1),
		.rxns – cell array of reaction names (size n x 1),
		.c - vector of objective coefficients (size n x 1),
		.lb - lower bounds on variables (size n x 1),
		.ub - upper bounds on variables (size n x 1),
		.b – right-hand side vector (size m x 1),
		.csense – vector indicating type of constraint (size m x 1)
		(‘E’ for equality, ‘G’ for greater than or equal, ‘L’ for less than or equal)

Folders and files:
-----------------

functions/:
clean_model.m: model preprocessing 
- removal blocked reactions
- reversible reaction split into two irreversible once
- calculation of species-complex matrix Y and complex-reaction matrix A (by deficieny and get_AY_matrix.r run in R)

find_balanced_complexes.m: calculate if complex is balanced given model structure

remove_balanced_complexes_any.m: return reduced network given original network and set of balanced complexes obtained from find_balanced_complexes.m using the motif applicable for any kinetics 

remove_balanced_complexes_MA.m: return reduced network given original network and set of balanced complexes obtained from find_balanced_complexes.m using the motif applicable for mass action kinetics

balanced_complex_removal_iterative_application.m:
- find balanced complexes
- removal according to specified kinetic assumption
- additional run of balanced complex identification and removal until no further balanced complexes are identified or no balanced complex can be removed due to motif restriction under any kinetics assumption


example/: 
example network for which balanced complexes are found (functions/find_balanced_complexes.m) and removed according to any kinetics (functions/remove_balanced_complexes_any.m) or mass action kinetics (functions/remove_balanced_complexes_MA.m)

networks/: 
Ecoli_kinetic_models/: Network data in .mat format for the two E. coli kinetic models used in this study obtained from their original publications (Khodayari et al. 2014, 2016)
GEM_xml/: .xml files of twelve genome-scale metabolic networks obtained from their original publications 

Results/: 
Results for the twelve genome-scale and two kinetic models obtained running run_balanced_complex_removal_GEM.m and run_balanced_complex_removal_kinetic.m as well as .mat files including tables obtained running result_tables.m


