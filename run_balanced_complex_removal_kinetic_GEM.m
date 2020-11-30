% Script to detect and remove balanced complexes in genome-scale mass-action kinetic 
% model if E. coli
clear
addpath('functions\')
mkdir('.','Results')
mkdir('Results\','Results_irreversibility_considered')
mkdir('Results\','final_models_irreversibility_considered')

name='Ecoli_kinetic_k-ecoli457';

load('networks/Ecoli_kinetic_models/genome_scale/Data.mat')

model.S = network_data.S_f_b;
model.b = zeros(size(model.S,1),1);
model.c = zeros(size(model.S,2),1);
model.mets = network_data.metab;
model.rxns = network_data.rxn_f_b;
model.lb = zeros(size(model.S,2),1);
model.ub = ones(size(model.S,2),1)*1000;
model.rxnNames = network_data.rxn_f_b;
model.metNames = network_data.metab;
model.rev = zeros(size(model.lb));
model.csense = repmat('E',size(model.b));

% model preprocessing, calculation of A and Y matrix
[model_r,out_r]=clean_model(name,model,'irreversibility_considered','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"');

% find and remove balanced complexes according to specified assumption
% about kinetics
disp('irreversibility_considered - any kinetic')
[MODEL_r,B_r,TRIVIAL_r,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_r,'any_kinetic');
save(strcat('Results/Results_irreversibility_considered/',name,'_any_kinetic.mat'),'MODEL_r','B_r','TRIVIAL_r','B_out','TRIVIAL_out')

disp('irreversibility_considered - mass action')
[MODEL_r,B_r,TRIVIAL_r,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_r,'mass_action');
save(strcat('Results/Results_irreversibility_considered/',name,'_mass_action.mat'),'MODEL_r','B_r','TRIVIAL_r','B_out','TRIVIAL_out')
