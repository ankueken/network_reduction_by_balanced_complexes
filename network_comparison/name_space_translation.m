% script to translate metabolite names into MNXref name space

clear
mkdir('core_comparison\')
files=dir('../Results/Results_irreversibility_considered/*mass_action.mat');
MNXref_met=readtable('MNXref-met-translation-table.csv');

species_list = {'A.niger','A. thaliana','C.reinhardtii','E. coli','M. acetivorans','M. barkeri',...
    'M. musculus','M. tuberculosis', '','','','','N. pharaonis','P. putida','T. maritima',};

for f=[4 8 9 10 11] % bacterial models
    disp(f)
    load(strcat('../Results/Results_irreversibility_considered\',files(f).name))
    
    mets_o=MODEL_r{1}.mets;
    
    % remove compartment information from metabolite name
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'[Cytosol]','');
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'[Extra_organism]','');
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'[Mitochondria]','');
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'_mt','');
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'_biomass','');
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'[intracellular]','');
    MODEL_r{1}.mets=strrep(MODEL_r{1}.mets,'[extracellular]','');
    
    for i=1:length(MODEL_r{1}.mets)
        MODEL_r{1}.mets{i}=MODEL_r{1}.mets{i}(1:end-3);
    end
     
    MetList=cellfun(@(x) strsplit(x,','), MODEL_r{1}.mets, 'UniformOutput', false);
    MetListLong=cellfun(@(x) strsplit(x,{'-',','}), MODEL_r{1}.metNames, 'UniformOutput',false);
    
    T=readtable('translation_GEM.xlsx','Sheet',species_list{f});
    MetList=T.mets;
    
    for i=1:length(MetList)
        trList=translateIDs(MetList{i}, 'met', MNXref_met, 'BiGG', 'MNXref',false);
        trList(strcmp(trList,'nan'))=[];
        if length(trList)>1
            MetListTranslated{i,1} = strjoin(trList,'|');
        elseif length(trList)==0
            trList=translateIDs(MetListLong{i}, 'met', MNXref_met, 'NAMES', 'MNXref',false);
            trList(strcmp(trList,'nan'))=[];
            if length(trList)>1
                MetListTranslated{i,1} = strjoin(trList,'|');
            elseif length(trList)==0
                trList=translateIDs(MetList{i}, 'met', MNXref_met, 'KEGG', 'MNXref',false);
                trList(strcmp(trList,'nan'))=[];
                if length(trList)>1
                    MetListTranslated{i,1} = strjoin(trList,'|');
                elseif length(trList)==0
                    MetListTranslated{i,1} = strjoin(MetList{i},'|');
                else
                    MetListTranslated{i,1} = trList{1};
                end
            else
                MetListTranslated{i,1} = trList{1};
            end
        else
            MetListTranslated{i,1} = trList{1};
        end
    end
    
    MODEL_r{1}.translation=(length(find(contains(MetListTranslated,'MNXM')))/length(MetListTranslated))*100;
    MODEL_r{1}.met_MNXref=MetListTranslated;
    
    [~,I]=intersect(mets_o,MODEL_r{end}.mets);
    MODEL_r{end}.met_MNXref=MetListTranslated(I);
    MODEL_r{end}.translation=(length(find(contains(MetListTranslated(I),'MNXM')))/length(MetListTranslated(I)))*100;
    
    save(strcat('core_comparison/',files(f).name))
    clearvars -except f files MNXref_met species_list
end


