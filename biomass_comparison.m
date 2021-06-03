% Script to detect and remove balanced complexes in twelve genome-scale metabolic
% networks for organism across all kingdoms of life
%
% To run further models add the respective name of the xml file
% located in networks/GEM_xml to species_list

clear
r=[];

addpath('functions\')
species_list={'A_niger_iMA871';'ArabidopsisCoreModel';'M_acetivorans_iMB745';
    'Ecoli2011-iJO1366';'M_musculus';'M_barkeri_iAF692';'T_maritima_iLJ478';'N_pharaonis';
    'C_reinhardtii_iCre1355_auto';'M_tuberculosis_iNJ661m';'P_putida_iJN746';'YeastGEM'};

pathToR = '"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"';
threshold=1e-9;

for s=1:length(species_list)
    
    name=species_list{s}
    model=readCbModel(strcat('networks/GEM_xml/',name,'.xml'));
    if s==8
        model.c(find(strcmp(model.rxns,'RM00001'))) = 1;
    end
    if all(model.c==0)
        keyboard
    end
    Results=load(strcat('Results/Results_irreversibility_considered/',name,'_any_kinetic'));
        
    %% Biomass production in original preprocesses model
    model=removeRxns(model,model.rxns(find(all(model.S==0))));
    model=removeMetabolites(model,model.mets(find(all(model.S'==0))));
    
    model.lb(model.c~=0)=0; % biomass reaction always considered irreversible
    model=convertToIrreversible(model);
    [~,Sol.f]=linprog(-model.c,[],[],model.S,model.b,model.lb,model.ub);
    
    model.lb(model.c~=0) = Sol.f*-0.1;
    [mini,maxi]=linprog_FVA(model,0.1);
    BLK=model.rxns(intersect(find(abs(mini)<threshold),find(abs(maxi)<threshold)));
    model=removeRxns(model,BLK);
    [~,Sol2.f]=linprog(-model.c,[],[],model.S,model.b,model.lb,model.ub);
    
    if abs(Sol2.f)>0.001
        
        save(strcat('Results/final_models_irreversibility_considered_biomass/',name,'_final_model.mat'))
        cd functions
        system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_irreversibility_considered_biomass/',name,'_final_model.mat')}));
        cd ..
        A=importdata(strcat('Results/final_models_irreversibility_considered_biomass/',name,'_final_model.A'));
        model.A=A.data;
        model.complexes=A.textdata(2:end);
        Y=importdata(strcat('Results/final_models_irreversibility_considered_biomass/',name,'_final_model.Y'));
        model.Y=Y.data;
        clear A Y
        save(strcat('Results/final_models_irreversibility_considered_biomass/',name,'_final_model.mat'))
        
        B=find_balanced_complexes(model,threshold);
        
        %% exclude balanced complexes including species part of biomass reaction
        species_in_biomass = find(model.S(:,model.c~=0)~=0);
        [~,complexes_with_biomass_species] = find(model.Y(species_in_biomass,:)~=0);
        B = setdiff(B,complexes_with_biomass_species);
        
        %% remove balanced complexes
        [model_any] = remove_balanced_complexes_any(model,B);
        [model_ma]  = remove_balanced_complexes_MA(model,B);
        
        size(model.S,1)-size(model_any.S,1)
        size(model.S,2)-size(model_any.S,2)
        size(model.A,1)-size(model_any.A,1)
        [~,Solany.f]=linprog(-model_any.c,[],[],model_any.S,model_any.b,model_any.lb,model_any.ub);
        [Sol2.f Solany.f]
        
        size(model.S,1)-size(model_ma.S,1)
        size(model.S,2)-size(model_ma.S,2)
        size(model.A,1)-size(model_ma.A,1)
        [~,Solma.f]=linprog(-model_ma.c,[],[],model_ma.S,model_ma.b,model_ma.lb,model_ma.ub);
        [Sol2.f Solma.f]
        
        save(strcat('Results/Results_irreversibility_considered_biomass/',name,'.mat'))
        
        clearvars -except s species_list threshold pathToR
    end 
end