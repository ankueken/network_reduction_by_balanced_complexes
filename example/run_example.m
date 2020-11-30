%% run network reduction based on balanced complexes for example network 
%% shown in example_figure.png

addpath('..\functions\')
mkdir('..','Results')
mkdir('../Results\','Results_irreversibility_considered')
mkdir('../Results\','final_models_irreversibility_considered')

clear
% example network
%        v1 v2 v3 v4 v5 v6 v7 v8 v9 v10
model.S=[-2  2  0  0  0  0  1 -1  0  0; % A
          1 -1 -1 -1  1  2 -2  0 -2  2; % B
          0  0  1 -1  1  0  0  0  0  0; % C
          0  0  0  1 -1 -1  0  0  0  0; % D
          0  0  0  0  0  0  1 -1  0  0; % E
          0  0  0  0  0  0  0  1  1 -1];% F
model.mets={'A';'B';'C';'D';'E';'F'};
model.rxns={'2A->B';'B->2A';'B->C';'B+C->D';'D->B+C';...
    'D->2B';'2B->A+E';'A+E->F';'2B->F';'F->2B'};
model.c=zeros(size(model.rxns));
model.lb=zeros(size(model.rxns));
model.ub=ones(size(model.rxns))*1000;
model.b=zeros(size(model.mets));
model.csense=repmat('E',size(model.b));

% get A and Y matrix
cd ..
model=clean_model('example_network',model,'irreversibility_considered','"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"');

% find balanced complexes
B = find_balanced_complexes(model);
model.complexes(B) 
% .complexes includes complex names in the format 
% stoichiometry*index of metabolite
% i.e. 
% >> model.complexes
% 
% ans =
% 
%   8×1 cell array
% 
%     {'2*1'    }  --> 2A --> balanced
%     {'1*2'    }  --> B
%     {'1*3'    }  --> C
%     {'1*2+1*3'}  --> B+C
%     {'1*4'    }  --> D --> balanced
%     {'2*2'    }  --> 2B
%     {'1*1+1*5'}  --> A+E --> balanced
%     {'1*6'    }  --> F --> balanced
%      

%% any kinetics
% remove balanced complexes according to any kinetics
[model_new_any,status_any] = remove_balanced_complexes_any(model,B);
% model_new_any includes the reduced model structure 
% status_any =
% 
%      1    -1     1     1 
%  indicating that complex B(2) corresponding to complex D could not be
%  removed - a complex with two out-going reactions
B_new_any = find_balanced_complexes(model_new_any);
% no further balanced complexes then complex D which cannot be removed
% under the assumption of any kinetics

%% mass action kinetics
% remove balanced complexes according to mass action kinetics
[model_new_MA,status_MA] = remove_balanced_complexes_MA(model,B);
% model_new_MA includes the reduced model structure 
% status_MA =
% 
%      1    1     1     1 
%  indicating that all complex could be removed
B_new_MA = find_balanced_complexes(model_new_MA);
% no further balanced complexes 
