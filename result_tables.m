%% Results shown in Supplementary Table S1
%
% - Result_table: table including number of species, number of reactions, 
% number of complexes, number trivially balanced complexes (removed), 
% number non-trivially balanced complexes (removed),
% number trivially balanced complexes (not removed, i.e. for any kinetics
% more than one outgoing reaction), number non-trivially balanced complexes (not removed)
%
% - Removed_species: table showing the number of species removed per round
% of balanced complex removal and the number of species with stoichiometry
% larger than one
%
% - Model_reduction: table showing the total percentage of model reduction in number of species
%

filename={'Results\Results_all_reversible\*any_kinetic.mat';
    'Results\Results_all_reversible\*mass_action.mat';
    'Results\Results_irreversibility_considered\*any_kinetic.mat';
    'Results\Results_irreversibility_considered\*mass_action.mat';
    'Results\Results_objective\*any_kinetic.mat';
    'Results\Results_objective\*mass_action.mat'};

for n=1:length(filename)
    clearvars -except filename n threshold
    Result_table=table();Removed_species=table();Model_reduction=table();
    files=dir(filename{n});
    
    for f=1:length(files)
        
        clearvars -except threshold n filename f files Result_table Removed_species Model_reduction
        load(strcat(files(f).folder,'\',files(f).name))
        name=files(f).name;
        W=whos('MODEL_*');
        eval(['MODEL=', genvarname(W.name),';']);
        W=whos('B_*');
        eval(['B=', genvarname(W.name),';']);
        W=whos('TRIVIAL_*');
        eval(['TRIVIAL=', genvarname(W.name),';']);
        
        for i=1:length(MODEL)
            
            Result_table=[Result_table;table({name},size(MODEL{i}.S,1),size(MODEL{i}.S,2),size(MODEL{i}.A,1),sum(TRIVIAL{i}),sum(TRIVIAL{i}==0),sum(TRIVIAL_out{i}),sum(TRIVIAL_out{i}==0),'VariableNames',{'Organism','Spiecies','Reactions','Complexes','TriviallyBalanced_removed','NonTriviallyBalanced_removed','TriviallyBalanced_not_removed','NonTriviallyBalanced_not_removed'})];
            
            % species removal per round of complex removal
            if i>=2 % if reduction was possible
                [removed_mets,IRM]=setdiff(MODEL{i-1}.mets,MODEL{i}.mets);
                % removed species with stoichiometry > 1
                [y,~]=find(abs(MODEL{i-1}.S(IRM,:))>1);
                Removed_species=[Removed_species;table({name},length(IRM),length(unique(y)),'VariableNames',{'Organism','species_removed','species_stoichiometry_larger_one_removed'})];
            end
            
        end
        
        x=find(strcmp(Result_table.Organism,name));
        Model_reduction=[Model_reduction;table({name},(1-(Result_table.Spiecies(x(end))./Result_table.Spiecies(x(1))))*100,'VariableNames',{'Organism','Model_reduction_percent'})];
    end
    
    savename=strrep(filename{n},'\*','_');
    save(savename)
end