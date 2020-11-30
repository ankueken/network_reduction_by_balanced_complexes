% script to compare the metabolites in reduced and 
% original models of Yeast and E. coli

clear
%%%%%%%%%%%%%%%%%%
% Yeast vs Ecoli %
%%%%%%%%%%%%%%%%%%

Yeast=load('Results_reversible\YeastGEM_mass_action.mat');
Ecoli=load('core_comparison\Ecoli2011-iJO1366_mass_action.mat');

length(intersect(Ecoli.MODEL_r{1}.met_MNXref,Yeast.MODEL_r{1}.metMetaNetXID))/...
  length(union(Ecoli.MODEL_r{1}.met_MNXref,Yeast.MODEL_r{1}.metMetaNetXID))

length(intersect(Ecoli.MODEL_r{end}.met_MNXref,Yeast.MODEL_r{end}.met_MNXref))/...
  length(union(Ecoli.MODEL_r{end}.met_MNXref,Yeast.MODEL_r{end}.met_MNXref))

[length(intersect(Ecoli.MODEL_r{end}.met_MNXref,Yeast.MODEL_r{end}.met_MNXref)) ...
length(intersect(Ecoli.MODEL_r{end}.met_MNXref,Yeast.MODEL_r{end}.met_MNXref))]


