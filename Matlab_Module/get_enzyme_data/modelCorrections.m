%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modelCorrections(model)
% Corrects various issues in yeast 7.
%
% INPUT:    A yeast model as a .mat structure.
% OUTPUT:   The corrected model.
%
% Benjamín J. Sánchez. Last edited: 2016-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = modelCorrections(model)

%Delete blocked rxns (LB = UB = 0):
to_remove = boolean((model.lb == 0).*(model.ub == 0));
model     = removeRxns(model,model.rxns(to_remove));

%Correct rev vector: true if LB < 0 & UB > 0 (and it's not an exchange reaction):
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || ...
       ~isempty(strfind(model.rxnNames{i},'exchange'))
        model.rev(i) = true;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%