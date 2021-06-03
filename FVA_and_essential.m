files = dir('Results\Results_irreversibility_considered_biomass\*.mat');
addpath('functions\')
for f=1:length(files)
    
    R=load(strcat(files(f).folder,'/',files(f).name));
    
    [~,Sol.f] = linprog(-R.model.c,[],[],R.model.S,R.model.b,R.model.lb,R.model.ub);
    [~,Sol_ma.f] = linprog(-R.model_ma.c,[],[],R.model_ma.S,R.model_ma.b,R.model_ma.lb,R.model_ma.ub);
    [~,Sol_any.f] = linprog(-R.model_any.c,[],[],R.model_any.S,R.model_any.b,R.model_any.lb,R.model_any.ub);
    
    if abs(Sol.f)>0

        [mini,maxi]=linprog_FVA(R.model,0.99);
        
        %% essential reactions
        for r=1:length(R.model.rxns)
            
            model_temp.S = R.model.S;
            model_temp.rxns = R.model.rxns;
            model_temp.mets = R.model.mets;
            model_temp.b = R.model.b;
            model_temp.c = R.model.c;
            model_temp.lb = R.model.lb;
            model_temp.ub = R.model.ub;
            model_temp.csense = R.model.csense;

            model_temp.lb(r) = 0;
            model_temp.ub(r) = 0;

            [~,Sol.f] = linprog(-model_temp.c,[],[],model_temp.S,model_temp.b,model_temp.lb,model_temp.ub);
                
            if ~isempty(Sol.f)
                essential(r,1) = abs(Sol.f)<1e-9;
            else
                essential(r,1) = 1;
            end
        end
    end
    
    if abs(Sol_ma.f)>0

        [mini_ma,maxi_ma]=linprog_FVA(R.model_ma,0.99);    
        
        for r=1:length(R.model_ma.rxns)
            model_temp.S = R.model_ma.S;
            model_temp.rxns = R.model_ma.rxns;
            model_temp.mets = R.model_ma.mets;
            model_temp.b = R.model_ma.b;
            model_temp.c = R.model_ma.c;
            model_temp.lb = R.model_ma.lb;
            model_temp.ub = R.model_ma.ub;
            model_temp.csense = R.model_ma.csense;

            model_temp.lb(r) = 0;
            model_temp.ub(r) = 0;

            [~,Sol.f] = linprog(-model_temp.c,[],[],model_temp.S,model_temp.b,model_temp.lb,model_temp.ub);
                         
            if ~isempty(Sol.f)
                essential_ma(r,1) = abs(Sol.f)<1e-9;
            else
                essential_ma(r,1) = 1;
            end
        end
    end
     
    if abs(Sol_any.f)>0

        [mini_any,maxi_any]=linprog_FVA(R.model_any,0.99);
        
        for r=1:length(R.model_any.rxns)
            model_temp.S = R.model_any.S;
            model_temp.rxns = R.model_any.rxns;
            model_temp.mets = R.model_any.mets;
            model_temp.b = R.model_any.b;
            model_temp.c = R.model_any.c;
            model_temp.lb = R.model_any.lb;
            model_temp.ub = R.model_any.ub;
            model_temp.csense = R.model_any.csense;

            model_temp.lb(r) = 0;
            model_temp.ub(r) = 0;

            [~,Sol.f] = linprog(-model_temp.c,[],[],model_temp.S,model_temp.b,model_temp.lb,model_temp.ub);
             
            if ~isempty(Sol.f)
                essential_any(r,1) = abs(Sol.f)<1e-9;
            else
                essential_any(r,1) = 1;
            end   
        end
        
        save(strcat('Results/Results_FVA_essential/',R.name,'.mat'))
        
        clearvars -except f files
        
    end
end
