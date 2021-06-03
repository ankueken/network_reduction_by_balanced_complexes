clear
load('Results\Results_irreversibility_considered\Ecoli_kinetic_k-ecoli457_mass_action.mat')

M1=find(sum(MODEL_r{1}.Y~=0,2)==1);
[~,C1]=find(MODEL_r{1}.Y(M1,:)~=0);
MODEL_r{1}.mets(M1);
MODEL_r{1}.complexes(C1);


BT = B_r{1}(TRIVIAL_r{1}==1);
BA = B_r{1};
[mets_T,comp_T]=find(MODEL_r{1}.Y(:,BT)~=0);
[mets_A,comp_A]=find(MODEL_r{1}.Y(:,BA)~=0);
[mets_T,BT(comp_T)'];
[mets_A,BA(comp_A)'];

clear
load('Results\Results_irreversibility_considered\Ecoli_kinetic_k-ecoli457_any_kinetic.mat')

MODEL_r{1}.complexes(B_r{1});

for i=1:length(B_r{1})

    rxns_in=find(MODEL_r{1}.A(B_r{1}(i),:)>0);
    rxns_out=find(MODEL_r{1}.A(B_r{1}(i),:)<0);

    [substrate_complexes,~]=find(MODEL_r{1}.A(:,rxns_in)<0);
    [product_complexes,~]=find(MODEL_r{1}.A(:,rxns_out)>0);

    if substrate_complexes~=product_complexes
        keyboard
    end
    
end