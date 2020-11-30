%% Supplementary Table S3
%  Number of metabolites per model compartment in original and reduced
%  model

files=dir('Results\Results_irreversibility_considered\*mass_action.mat');

for f=1:length(files)
    
    model = load(strcat('Results\Results_irreversibility_considered\',files(f).name));
    
    if isfield(model.MODEL_r{1},'comps')
        disp(files(f).name) 
        Mets_per_compartment = nan(length(model.MODEL_r{1}.comps),2);
        for i=1:length(model.MODEL_r{1}.comps)
            
            Mets_per_compartment(i,1) = sum(contains(model.MODEL_r{1}.mets,strcat('[',model.MODEL_r{1}.comps(i),']')));
            Mets_per_compartment(i,2) = sum(contains(model.MODEL_r{end}.mets,strcat('[',model.MODEL_r{1}.comps(i),']')));
            
        end
        
        Mets_per_compartment(:,3) = (1-(Mets_per_compartment(:,2)./Mets_per_compartment(:,1)))*100;
        T{f,1}.model = files(f).name;
        T{f,1}.met_per_comp=table(model.MODEL_r{1}.compNames,Mets_per_compartment);
        
    end
end