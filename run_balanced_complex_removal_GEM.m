% Script to detect and remove balanced complexes in twelve genome-scale metabolic
% networks for organism across all kingdoms of life
%
% To run further models add the respective name of the xml file 
% located in networks/GEM_xml to species_list 

clear
changeCobraSolver('glpk')

addpath('functions\')
mkdir('.','Results')
mkdir('Results\','Results_objective')
mkdir('Results\','Results_irreversibility_considered')
mkdir('Results\','Results_all_reversible')
mkdir('Results\','final_models_objective')
mkdir('Results\','final_models_irreversibility_considered')
mkdir('Results\','final_models_all_reversible')

species_list={'A_niger_iMA871';'ArabidopsisCoreModel';'M_acetivorans_iMB745';
    'Ecoli2011-iJO1366';'M_musculus';'M_barkeri_iAF692';'T_maritima_iLJ478';'N_pharaonis';
    'C_reinhardtii_iCre1355_auto';'M_tuberculosis_iNJ661m';'P_putida_iJN746';'YeastGEM'};

for s=1:length(species_list)

    name=species_list{s}

    %% load model
    model=readCbModel(strcat('networks/GEM_xml/',name,'.xml'));

    %% create models for each scenario
        % 'all_reversible' - all reactions considered reversible
        % 'irreversibility_considered' - reaction irreversibility considered
        % 'objective' - maximize objective function 
    [model_n,BLK]=clean_model(name,model,'all_reversible','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"');
    [model_r]=clean_model(name,model,'irreversibility_considered','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"',BLK);
    [model_o]=clean_model(name,model,'objective','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"',BLK);

    %% find and remove balanced complexes for each scenario
    disp('all_reversible - any kinetic')
    [MODEL_n,B_n,TRIVIAL_n,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_n,'any_kinetic');
    save(strcat('Results/Results_all_reversible/',name,'_any_kinetic.mat'),'MODEL_n','B_n','TRIVIAL_n','B_out','TRIVIAL_out')
    
    disp('all_reversible - mass action')
    [MODEL_n,B_n,TRIVIAL_n,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_n,'mass_action');
    save(strcat('Results/Results_all_reversible/',name,'_mass_action.mat'),'MODEL_n','B_n','TRIVIAL_n','B_out','TRIVIAL_out')
         
    disp('irreversibility_considered - any kinetic')
    [MODEL_r,B_r,TRIVIAL_r,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_r,'any_kinetic');
    save(strcat('Results/Results_irreversibility_considered/',name,'_any_kinetic.mat'),'MODEL_r','B_r','TRIVIAL_r','B_out','TRIVIAL_out')
        
    disp('irreversibility_considered - mass action')
    [MODEL_r,B_r,TRIVIAL_r,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_r,'mass_action');
    save(strcat('Results/Results_irreversibility_considered/',name,'_mass_action.mat'),'MODEL_r','B_r','TRIVIAL_r','B_out','TRIVIAL_out')

    disp('objective - any kinetic')
    [MODEL_o,B_o,TRIVIAL_o,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_o,'any_kinetic');
    save(strcat('Results/Results_objective/',name,'_any_kinetic.mat'),'MODEL_o','B_o','TRIVIAL_o','B_out','TRIVIAL_out')
    
    disp('objective - mass action')
    [MODEL_o,B_o,TRIVIAL_o,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model_o,'mass_action');        
    save(strcat('Results/Results_objective/',name,'_mass_action.mat'),'MODEL_o','B_o','TRIVIAL_o','B_out','TRIVIAL_out')
  
    clearvars -except s species_list
end