function [model_new,status] = remove_balanced_complexes_MA(model,B)
% function to remove a set of given balanced complexes according to mass
% action kinetics assumption
%
% [model_new,status] = remove_balanced_complexes_MA(model,B)
%
% Input
%   model: struct with at least following fields
%           .c      objective vectore
%           .S      stoichiometric matrix
%           .b      right-hand side vector
%           .A      complex-reaction matrix
%           .Y      species-complex matrix
%           .lb     lower bound on flux
%           .ub     upper bound on flux
%           .rxns   cell array of reaction names
%           .mets   cell array of metabolite names
%           .complexes cell array of complex names
%           .csense vectore of size .b indicating type of constraint
%                   according to cobra model definition (type 'E','G' or 'L')
%   B: indices of balanced complexes
%
% Output
%   model_new: reduced model
%   status: vector indicating if problems appeared
%            1  successful removal
%           -1  removal lead to infeasible model indicating numeric problems
%
% -------------------------------------------------

B_name=model.complexes(B);
F=fieldnames(model);
RF=setdiff(F,{'S','rxns','mets','A','Y','complexes','b','c','csense','lb','ub'});
model_new=rmfield(model,RF);
status=ones(size(B));

for b=1:length(B)
    model_new_temp=model_new; % if model is infeasible after removal we reset
    
    id=find(strcmp(model_new.complexes,B_name(b))); % find balanced complex ids in current model
    
    out_rxns = find(model_new.A(id,:)<0); % consuming reactions
    in_rxns = find(model_new.A(id,:)>0); % producing reactions
    
    % get products of out_rxns
    P.id=[]; % index
    P.stoich=[]; % stoichiometry
    P.out_rxns=[]; % track corresponding reaction
    P.lb=[];P.ub=[];
    for o=1:length(out_rxns)
        P.id=[P.id find(model_new.A(:,out_rxns(o))>0)];
        P.stoich=[P.stoich model_new.A(:,out_rxns(o))];
        P.out_rxns=[P.out_rxns out_rxns(o)];
        P.lb=[P.lb; model_new.lb(out_rxns(o))];
        P.ub=[P.ub; model_new.ub(out_rxns(o))];
    end
    
    % get substrates of in_rxns
    S.id=[]; % index
    S.stoich=[]; % stoichiometry
    S.in_rxns=[]; % track corresponding reaction
    S.lb=[];S.ub=[];
    for o=1:length(in_rxns)
        S.id=[S.id find(model_new.A(:,in_rxns(o))<0)];
        S.stoich=[S.stoich model_new.A(:,in_rxns(o))];
        S.in_rxns=[S.in_rxns in_rxns(o)];
        S.lb=[S.lb; model_new.lb(in_rxns(o))];
        S.ub=[S.ub; model_new.ub(in_rxns(o))];
    end
    
    % remove balanced complex
    model_new.Y(:,id) = [];
    model_new.complexes(id) = [];
    model_new.A(id,:) = [];
    
    % check for zero rows and columns
    m_remove = find(all(model_new.Y'==0));
    model_new.b(m_remove)=[];
    model_new.mets(m_remove)=[];
    model_new.Y(m_remove,:)=[];
    model_new.S(m_remove,:)=[];
    model_new.csense(m_remove)=[];
    
    %% model rewriting
    % internal complex 
    if ~isempty(P.id) && ~isempty(S.id)
        
        new_rxns=combvec(1:length(S.id),1:length(P.id));
        
        rxn_set=[S.id(new_rxns(1,:)); P.id(new_rxns(2,:))];
        new_rxns(:,rxn_set(1,:)==rxn_set(2,:))=[];
        rxn_set(:,rxn_set(1,:)==rxn_set(2,:))=[];
        
        if size(new_rxns,2)>1
            comb=nchoosek(1:size(new_rxns,2),2);
            remove_set=[];
            for i=1:size(comb,1)
                if all(sort(rxn_set(:,comb(i,1)))==sort(rxn_set(:,comb(i,2))))
                    remove_set= [remove_set comb(i,1)];
                end
            end
            
            new_rxns(:,remove_set)=[];
            rxn_set(:,remove_set)=[];
        end
        
        P.stoich(id,:) = [];
        S.stoich(id,:) = [];
        
        unaffected_rxns = setdiff(1:size(model_new.A,2),[S.in_rxns P.out_rxns]);
        
        A_r = model_new.A(:,unaffected_rxns);
        rxns_r = model_new.rxns(unaffected_rxns);
        lb_r = model_new.lb(unaffected_rxns);
        ub_r = model_new.ub(unaffected_rxns);
        c_r = model_new.c(unaffected_rxns);
        
        for o=1:size(new_rxns,2)
            A_r(:,end+1)=S.stoich(:,new_rxns(1,o))+P.stoich(:,new_rxns(2,o));
            rxns_r(end+1)=strcat(model_new.rxns(S.in_rxns(new_rxns(1,o))),'+',model_new.rxns(P.out_rxns(new_rxns(2,o))));
            c_r(end+1) = model_new.c(S.in_rxns(new_rxns(1,o)))+model_new.c(P.out_rxns(new_rxns(2,o)));
            lb_r(end+1) = max([S.lb(new_rxns(1,o)) P.lb(new_rxns(2,o))]);
            ub_r(end+1) = min([S.ub(new_rxns(1,o)) P.ub(new_rxns(2,o))]);
        end
        
    % complex metabolites exported to environment
    elseif isempty(P.id) && ~isempty(S.id)
        
        S.stoich(id,:) = [];
        
        unaffected_rxns = setdiff(1:size(model_new.A,2),[S.in_rxns]);
        
        A_r = model_new.A(:,unaffected_rxns);
        rxns_r = model_new.rxns(unaffected_rxns);
        lb_r = model_new.lb(unaffected_rxns);
        ub_r = model_new.ub(unaffected_rxns);
        c_r = model_new.c(unaffected_rxns);
        
        A_r(:,end+1:end+size(S.stoich,2))=S.stoich;
        rxns_r(end+1:end+size(S.stoich,2))=model_new.rxns(S.in_rxns);
        lb_r(end+1:end+size(S.stoich,2))=S.lb;
        ub_r(end+1:end+size(S.stoich,2))=S.ub;
        c_r(end+1:end+size(S.stoich,2))=model_new.c(S.in_rxns);
        
    % complex metabolites imported from environment
    elseif ~isempty(P.id) && isempty(S.id)
      
        P.stoich(id,:) = [];
        
        unaffected_rxns = setdiff(1:size(model_new.A,2),[P.out_rxns]);
        
        A_r = model_new.A(:,unaffected_rxns);
        rxns_r = model_new.rxns(unaffected_rxns);
        lb_r = model_new.lb(unaffected_rxns);
        ub_r = model_new.ub(unaffected_rxns);
        c_r = model_new.c(unaffected_rxns);
        
        A_r(:,end+1:end+size(P.stoich,2))=P.stoich;
        rxns_r(end+1:end+size(P.stoich,2))=model_new.rxns(P.out_rxns);
        lb_r(end+1:end+size(P.stoich,2))=P.lb;
        ub_r(end+1:end+size(P.stoich,2))=P.ub;
        c_r(end+1:end+size(P.stoich,2))=model_new.c(P.out_rxns);
        
    else
        A_r = model_new.A;
        rxns_r = model_new.rxns;
        lb_r = model_new.lb;
        ub_r = model_new.ub;
        c_r = model_new.c;
    end
    
    model_new.A = A_r;
    model_new.rxns = rxns_r;
    model_new.lb = lb_r;
    model_new.ub = ub_r;
    model_new.c = c_r;
    model_new.S = model_new.Y*model_new.A;
    
    %% remove zero rows and columns
    r_remove=find(all(model_new.A==0));
    c_remove=find(all(model_new.A'==0));
    m_remove=find(all(model_new.Y'==0));
    
    while ~isempty(r_remove) || ~isempty(c_remove) || ~isempty(m_remove)
        model_new.S(:,r_remove)=[];
        model_new.c(r_remove)=[];
        model_new.lb(r_remove)=[];
        model_new.ub(r_remove)=[];
        model_new.rxns(r_remove)=[];
        model_new.A(:,r_remove)=[];
        
        model_new.A(c_remove,:)=[];
        model_new.complexes(c_remove)=[];
        model_new.Y(:,c_remove)=[];
        
        m_remove=find(all(model_new.Y'==0));
   
        model_new.b(m_remove)=[];
        model_new.mets(m_remove)=[];
        model_new.Y(m_remove,:)=[];
        model_new.S(m_remove,:)=[];
        model_new.csense(m_remove)=[];
        
        r_remove=find(all(model_new.A==0));
        c_remove=find(all(model_new.A'==0));
        m_remove=find(all(model_new.Y'==0));
    end
    
end

end