clear
files = dir('Results\Results_FVA_essential\*.mat');

for f=1:length(files)
    clearvars -except f files unaltered_rxns Equal_FVA Different_FVA Specificity ...
        Sensitivity Accuracy EN Lower_FVA Larger_FVA
    D = load(strcat(files(f).folder,'/',files(f).name));
    
    % how many reactions are unchanged?
    if isfield(D,'model_any')
        [~,idx_orig,idx_any] = intersect(D.model.rxns,model_any.rxns);
        unaltered_rxns(f,1) = length(idx_any)/length(D.model_any.rxns);
    elseif isfield(D,'R')
        [~,idx_orig,idx_any] = intersect(D.R.model.rxns,D.R.model_any.rxns);
        unaltered_rxns(f,1) = length(idx_any)/length(D.R.model_any.rxns);
    end
    
    % FVA unaltered reactions
    D.maxi = round(D.maxi,6);
    D.mini = round(D.mini,6);
    D.maxi_any = round(D.maxi_any,6);
    D.mini_any = round(D.mini_any,6);
    FVA_equal_upper_bound = D.maxi(idx_orig) == D.maxi_any(idx_any);
    FVA_equal_lower_bound = D.mini(idx_orig) == D.mini_any(idx_any);
    FVA_lower_upper_bound = D.maxi(idx_orig) > D.maxi_any(idx_any);
    FVA_lower_lower_bound = D.mini(idx_orig) > D.mini_any(idx_any);
    FVA_higher_upper_bound = D.maxi(idx_orig) < D.maxi_any(idx_any);
    FVA_higher_lower_bound = D.mini(idx_orig) < D.mini_any(idx_any);
    Equal_FVA(f,1) = sum(sum(FVA_equal_upper_bound+FVA_equal_lower_bound,2)==2)/length(FVA_equal_lower_bound);
    Larger_FVA(f,1) = length([intersect(find(FVA_lower_lower_bound),find(FVA_higher_upper_bound));...
                                intersect(find(FVA_equal_lower_bound),find(FVA_higher_upper_bound));...
                                intersect(find(FVA_lower_lower_bound),find(FVA_equal_upper_bound))])/(length(FVA_equal_lower_bound)-sum(sum(FVA_equal_upper_bound+FVA_equal_lower_bound,2)==2));
    Lower_FVA(f,1) = length([intersect(find(FVA_higher_lower_bound),find(FVA_lower_upper_bound));...
                                intersect(find(FVA_equal_lower_bound),find(FVA_lower_upper_bound));...
                                intersect(find(FVA_higher_lower_bound),find(FVA_equal_upper_bound))])/(length(FVA_equal_lower_bound)-sum(sum(FVA_equal_upper_bound+FVA_equal_lower_bound,2)==2));
    % essential reactions among unaltered reactions
    % orig_yes | orig_no
    % ------------------
    % red_yes  | red_no
    TP = length(intersect(find(D.essential(idx_orig)==1),find(D.essential_any(idx_any)==1)));
    FP = length(intersect(find(D.essential(idx_orig)==0),find(D.essential_any(idx_any)==1)));
    TN = length(intersect(find(D.essential(idx_orig)==0),find(D.essential_any(idx_any)==0)));
    FN = length(intersect(find(D.essential(idx_orig)==1),find(D.essential_any(idx_any)==0)));
    P = sum(D.essential(idx_orig)==1);
    N = sum(D.essential(idx_orig)==0);
    
    % specificity
    Specificity(f,1) = TN/N;
    % sensitivity
    Sensitivity(f,1) = TP/P;
    % accuracy
    Accuracy(f,1) = (TP+TN)/(TP+FP+TN+FN);
    
    EN(f,1)=sum(D.essential_any(setdiff(1:length(D.essential_any),idx_any))==1)/length(setdiff(1:length(D.essential_any),idx_any));
end

labels={'\itA. niger iMA871'
'\itA. thaliana AraCore'
'\itC. reinhardtii iCre1355'
'\itE. coli iJO1366'
'\itM. acetivorans iMB745'
'\itM. barkeri iAF692'
'\itM. musculus'
'\itM. tuberculosis iNJ661m'
'\itN. pharaonis'
'\itP. putida iJN746'
'\itT. maritima iLJ478'
'\itS. cerevisiae Yeast8'};

bar(Equal_FVA*100)
set(gca,'XTick',1:12,'XTickLabels',labels,'XTickLabelRotation',45)
ylabel('Percentage mapped reactions with unaltered flux range')

