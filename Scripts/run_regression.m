function [betamat,tmat,pmat,rsq,summed_betas,summed_tstats,summed_pvals,corrs] = run_regression(rmats,factors,interaction)
%UNTITLED2 Summary of this function goes here
%   rmats = individual subject corrrelation matrices (roi,roi,subject)
%   factors = vector or array with factor to include in model 
%             each row is an observation
%             each column is a variable 
%   interaction = whether or not to test interaction terms (1=yes,0=no)
%

% determine number of variables to initiate matrices
varnum = size(factors,2);

betamat = zeros(size(rmats,1),size(rmats,2),varnum);
tmat = zeros(size(rmats,1),size(rmats,2),varnum);
pmat = zeros(size(rmats,1),size(rmats,2),varnum);
corrs = zeros(size(rmats,1),size(rmats,2),varnum);
rsq = zeros(size(rmats,1),size(rmats,2));

if varnum == 1 && interaction == 0
    disp('Simple Linear Regression')
    for i = 1:(size(rmats,1)-1)
        disp(['On ROI ' num2str(i)])
        for j = i+1:size(rmats,1)
            mdl = fitlm(factors,squeeze(rmats(i,j,:)));
            %rsq(i,j) = mdl.Rsquared.Ordinary;
            %betamat(i,j) = table2array(a.Coefficients(2,1)); % beta
            tmat(i,j) = table2array(mdl.Coefficients(2,3)); % t stat
            %pmat(i,j) = table2array(a.Coefficients(2,4));
            corrs(i,j) = corr(squeeze(rmats(i,j,:)),factors);
        end
    end
    
elseif varnum > 1 && interaction == 0
    disp('Multiple Linear Regression: NO interaction')
    for i = 1:(size(rmats,1)-1)
        disp(['On ROI ' num2str(i)])
        for j = i+1:size(rmats,1)
            mdl = fitlm(factors,squeeze(rmats(i,j,:)));
            rsq(i,j) = mdl.Rsquared.Ordinary;
            %[b,bint,r,rint,stats] = regress(squeeze(rmats(i,j,:)),factors);
            for v = 1:varnum
                betamat(i,j,v) = table2array(mdl.Coefficients(v+1,1)); % beta
                tmat(i,j,v) = table2array(mdl.Coefficients(v+1,3)); % t stat
                pmat(i,j,v) = table2array(mdl.Coefficients(v+1,4)); % p value 
            end

        end
    end
    
elseif varnum > 1 && interaction == 1
    disp('Multiple Linear Regression: YES interaction')
    
    for i = 1:(size(rmats,1)-1)
        disp(['On ROI ' num2str(i)])  
        for j = i+1:size(rmats,1)
            tbl = table(factors,squeeze(rmats(i,j,:)));
            mdl = fitlm(factors,squeeze(rmats(i,j,:)));
            rsq(i,j) = mdl.Rsquared.Ordinary;
            for v = 1:varnum
                betamat(i,j,v) = table2array(mdl.Coefficients(v+1,1)); % beta
                tmat(i,j,v) = table2array(mdl.Coefficients(v+1,3)); % t stat
                pmat(i,j,v) = table2array(mdl.Coefficients(v+1,4)); % p value 
            end
        end
     end
     
else
    error('Cannot run an interaction on one term!')
end

% Sum down columns for composite ROI scores 

if varnum == 1
    betamat = betamat + betamat';
    tmat = tmat + tmat';
    pmat = pmat + pmat';
    corrs = corrs + corrs';
    rsq = rsq + rsq';
    summed_betas = sum(betamat)';
    summed_tstats = sum(tmat)';
    summed_pvals = sum(pmat)';
else
    rsq = rsq + rsq';
    for v = 1:varnum
        betamat(:,:,v) = betamat(:,:,v) + betamat(:,:,v)';
        tmat(:,:,v) = tmat(:,:,v) + tmat(:,:,v)';
        pmat(:,:,v) = pmat(:,:,v) + pmat(:,:,v)';
        summed_betas(:,v) = sum(squeeze(betamat(:,:,v)))';
        summed_tstats(:,v) = sum(squeeze(tmat(:,:,v)))';
        summed_pvals(:,v) = sum(squeeze(pmat(:,:,v)))';
    end
end


end

