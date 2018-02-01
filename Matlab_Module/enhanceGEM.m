%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)
%
% Benjamin J. Sanchez. Last edited: 2017-04-12
%Ivan Domenzain.       Last edited: 2018-01-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)

%Provide your organism scientific name
org_name      = 'escherichia coli';
org_code      = 'eco';
format short e

if strcmp(toolbox,'COBRA')
   initCobraToolbox
end

%Add some RAVEN fields for easier visualization later on:
cd get_enzyme_data
model = modelCorrections(model);
model = standardizeModel(model,toolbox);
model.csense = model.csense';

%Retrieve kcats & MWs for each rxn in model:
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data, org_name);
cd ../../Models
save([org_code '_enzData.mat'],'model_data','kcats')
%Integrate enzymes in the model:
cd ../Matlab_Module/change_model
ecModel = readKcatData(model_data,kcats);

%Constrain model to batch conditions:
cd ../limit_proteins
sigma         = 0.5;      %Optimized for glucose
Ptot          = 0.55;     %Assumed constant
ecModel_batch = constrainEnzymes(ecModel,Ptot,sigma);

%Save output models:
cd ../../models
save([name '.mat'],'ecModel','model_data','kcats')
save([name '_batch.mat'],'ecModel_batch')
cd ../Matlab_Module

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%