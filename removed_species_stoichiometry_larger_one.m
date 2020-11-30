%% Supplementary Table S2
% species removed during network reduction of stoichiometry greater than
% one for networks of organisms across kingdoms

% constraints scenario: irreversibility considered
% output in stoichiometry_g1_species
%           column 1 - species of stoichiometry > 1 removed under any
%           kinetics
%           column 2 - species of stoichiometry > 1 removed under mass
%           action kinetics
%           column 3 - species names
clear

filename={'Results\Results_irreversibility_considered\*any_kinetic.mat';
    'Results\Results_irreversibility_considered\*mass_action.mat'};

for n=1:length(filename)
    clearvars -except filename n stoichiometry_g1_species
    files=dir(filename{n});
   
    for f=1:length(files)
        
        clearvars -except n filename f files stoichiometry_g1_species
        load(strcat(files(f).folder,'\',files(f).name))
        name=files(f).name;
        W=whos('MODEL_*');
        eval(['MODEL=', genvarname(W.name)]);
        W=whos('B_*');
        eval(['B=', genvarname(W.name)]);
        W=whos('TRIVIAL_*');
        eval(['TRIVIAL=', genvarname(W.name)]);
        
        if length(MODEL)>1
            [removed_mets,IRM]=setdiff(MODEL{1}.mets,MODEL{end}.mets);
            [y,~]=find(abs(MODEL{1}.S(IRM,:))>1);
            stoichiometry_g1_species{f,n}=MODEL{1}.metNames(IRM(unique(y)));
            stoichiometry_g1_species{f,3}=files(f).name;
        end
    end
end
