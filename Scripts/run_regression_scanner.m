function [rsq,meanrsq] = run_regression_scanner(rmats,factors,scanner)
%UNTITLED2 Summary of this function goes here
%   rmats = individual subject corrrelation matrices (roi,roi,subject)
%   factors = vector or array with factor to include in model 
%             each row is an observation
%             each column is a variable 
%   interaction = whether or not to test interaction terms (1=yes,0=no)
%

% determine number of variables to initiate matrices
varnum = size(factors,2) + 1;

rsq = zeros(size(rmats,1),size(rmats,2));

for i = 1:(size(rmats,1)-1)
    disp(['On ROI ' num2str(i)])  
    for j = i+1:size(rmats,1)
        tbl = table(factors,scanner,squeeze(rmats(i,j,:)));
        tbl.Properties.VariableNames{3} = 'corrs';
        mdl_noint = fitlm(tbl,'corrs~factors+scanner');
        mdlint = fitlm(tbl,'corrs~factors+scanner+factors:scanner');
        rsq(i,j) = mdlint.Rsquared.Ordinary - mdl_noint.Rsquared.Ordinary;
    end
 end

% Sum down columns for composite ROI scores 

rsq = rsq + rsq';
meanrsq = mean(rsq);

end

