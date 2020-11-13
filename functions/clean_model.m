function [model,BLK] = clean_model(name,model,constraints,pathToR,BLK,threshold)
% function to preprocess models used to identify balanced complexes
%
% [model,BLK] = clean_model(name,model,constraints,threshold,pathToR,BLK)
%
% Input
%       name: name used to save model files
%       model: struct with at least following fields
%           .c      objective vectore
%           .S      stoichiometric matrix
%           .b      right-hand side vector
%           .lb     lower bound on flux
%           .ub     upper bound on flux
%           .rxns   cell array of reaction names
%           .mets   cell array of metabolite names
%       constraints: specify scenario
%                    'all_reversible' - all reactions considered reversible
%                    'irreversibility_considered' - reaction irreversibility considered
%                    'objective' - maximize objective function 
%                                  (i.e. biomass)
%       pathToR: specify path to R in the form given below 
%                e.g. '"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"'
%       BLK: indices of reactions blocked in the original model (optional)
%       threshold: value when to consider complex balanced (e.g. 1e-9)
%                  (optional)
%
% Output
%       model: preprocessed model, extended by fields .A and .Y
%              no blocked reactions, reversible reactions split
%              constraints according to specified scenario
%              ('none','rev','obj')
%       BLK: indices of reactions blocked in the original model
%
% Step 1:
%   check if there are zero rows or columns in the stoichiometric matrix
%   and remove them
% Step 2:
%   identify and remove blocked reactions
% Step 3:
%   set case (all reversible, irreversibility considered, objective)
% Step 4:
%   convert to irreversibility
% Step 5:
%   calculate A and Y matrix
%
% -----------------------------------------------------------
if nargin<6
    threshold=1e-9;
end

% Step 1. Remove zero rows / columns
model=removeRxns(model,model.rxns(find(all(model.S==0))));
model=removeMetabolites(model,model.mets(find(all(model.S'==0))));

% Step 2. Find and remove blocked reactions
if ~exist('BLK','var')
    [mini,maxi]=linprog_FVA(model,0);
    BLK=model.rxns(intersect(find(abs(mini)<threshold),find(abs(maxi)<threshold)));
end
model=removeRxns(model,BLK);

% Step 3-4. convert to irreversibility and set constraints
if all(model.c==0) && strcmp(constraints,'objective')
    obj=find(contains(model.rxns,'bio','IgnoreCase',true));
    % N. pharaonis: model.c(132)=1; find(strcmp(model.rxns,'RM00001'))
    if isempty(obj)
        disp('Please specify flux to maximize in model.c')
        keyboard
    else
        model.c(obj(1))=1;
    end
end

lb_rev=model.lb;ub_rev=model.ub;

% Convert to irreversible
model.ub(:)=max(model.ub);
model.lb(:)=min(model.lb);
model=convertToIrreversible(model);
if min(lb_rev)<0
    model.ub=[ub_rev;abs(lb_rev)];
end

% Remove reversibility constraint
if strcmp(constraints,'all_reversible')
    model.ub(:)=max(model.ub);
    model.lb(:)=min(model.lb);
end

% Set objective to maximize
if strcmp(constraints,'objective')
    Sol=optimizeCbModel(model);
    model.lb(model.c~=0)=Sol.f*0.99;
    model.ub(model.c~=0)=Sol.f*0.99;
end

model.c(:)=0;

% Step 5. Compute A and Y matrix and save model

if strcmp(constraints,'all_reversible')
    save(strcat('Results/final_models_all_reversible/',name,'_final_model.mat'))
    cd functions
    system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_all_reversible/',name,'_final_model.mat')}));
    cd ..
    A=importdata(strcat('Results/final_models_all_reversible/',name,'_final_model.A'));
    model.A=A.data;
    model.complexes=A.textdata(2:end);
    Y=importdata(strcat('Results/final_models_all_reversible/',name,'_final_model.Y'));
    model.Y=Y.data;
    clear A Y
    save(strcat('Results/final_models_all_reversible/',name,'_final_model.mat'))
    
elseif strcmp(constraints,'irreversibility_considered')
    save(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.mat'))
    cd functions
    system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_irreversibility_considered/',name,'_final_model.mat')}));
    cd ..
    A=importdata(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.A'));
    model.A=A.data;
    model.complexes=A.textdata(2:end);
    Y=importdata(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.Y'));
    model.Y=Y.data;
    clear A Y
    save(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.mat'))
    
elseif strcmp(constraints,'objective')
    save(strcat('Results/final_models_objective/',name,'_final_model.mat'))
    cd functions
    system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_objective/',name,'_final_model.mat')}));
    cd ..
    A=importdata(strcat('Results/final_models_objective/',name,'_final_model.A'));
    model.A=A.data;
    model.complexes=A.textdata(2:end);
    Y=importdata(strcat('Results/final_models_objective/',name,'_final_model.Y'));
    model.Y=Y.data;
    clear A Y
    save(strcat('Results/final_models_objective/',name,'_final_model.mat'))
end

end