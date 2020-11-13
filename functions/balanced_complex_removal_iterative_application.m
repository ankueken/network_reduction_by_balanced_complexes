function [MODEL,B,TRIVIAL,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model,assumption,threshold)
% function for iterative identification and removal of balanced complexes
%
% [MODEL,B,TRIVIAL,STATUS,B_out,TRIVIAL_out]=balanced_complex_removal_iterative_application(model,assumption,threshold)
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
%   assumption: 'any_kinetic' or 'mass_action'
%   threshold: value when to consider complex balanced (e.g. 1e-9)
%              (optional)
%
% Output
%   MODEL: cell structure of reduced models per round
%   B: cell structure of indices of balanced complexes per round
%   TRIVIAL: cell structure of length B indicating if complex is trivially or non-trivially balanced per round
%               1 - trivial, 0 - non-trivial     
%   B_out: cell structure of indices of complexes identified as balanced,
%          but not removed per round
%   TRIVIAL_out: cell structure of length B_out indicating if complex is trivially or non-trivially balanced per round
%               1 - trivial, 0 - non-trivial     
%
% ------------------------------------------------------

    if nargin<3
        threshold=1e-9;
    end
    
    status=1;count=0;
    MODEL=cell(0,0);B=cell(0,0);TRIVIAL=cell(0,0);STATUS=cell(0,0);
    B_out=cell(0,0);TRIVIAL_out=cell(0,0);
    
    while ~isempty(find(status>0))
        count=count+1;
        MODEL{count}=model;
        
        % 1. find balanced complexes
        B{count}=find_balanced_complexes(model,threshold);

        % 2. classify as trivial and non-trivial
        %    balanced complexes including species of degree 1 
        %    (appearing in a single complex) are trivially balanced
        degree_species=sum(model.Y~=0,2);
        [~,temp_var]=find(model.Y(find(degree_species==1),B{count})~=0);
        TRIVIAL{count}=zeros(size(B{count}));
        TRIVIAL{count}(temp_var)=1;

        % 3. remove balanced complexes
        if strcmp(assumption,'any_kinetic')
            % 3.1 under any kinetic
            [model,status] = remove_balanced_complexes_any(model,B{count});
            STATUS{count}=status;
        elseif strcmp(assumption,'mass_action')
            % 3.1 under the assumption of mass action
            [model,status] = remove_balanced_complexes_MA(model,B{count});
            STATUS{count}=status;
        else
            error('foo:bar','\bUnknown assumption about kinetic.\n Assumption can be: \n - any_kinetic \n - mass_action')
        end
        B_out{count}=B{count}(STATUS{count}<0);
        TRIVIAL_out{count}=TRIVIAL{count}(STATUS{count}<0);
        B{count}(STATUS{count}<0) = [];
        TRIVIAL{count}(STATUS{count}<0) = [];
    end
end