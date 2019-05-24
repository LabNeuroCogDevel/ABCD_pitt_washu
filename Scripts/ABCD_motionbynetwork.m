%% All

for n = 1:size(MatchingCorrMats,3)
    ThisSub = squeeze(MatchingCorrMats(:,:,n));
    Mean(n,1) = mean(mean(ThisSub(logical(triu(ThisSub)))));
end
scatter(MatchingSubsFD(FDIdx),Mean);lsline
LinearModel.fit(MatchingSubsFD(FDIdx),Mean)

%% By Network 
NetworkFCFDCorr = zeros(13,13);
for net = 1:13
    for nettwo = 1:13
        if net == nettwo
            Subs = MatchingCorrMats(NetworkIds==net,NetworkIds==nettwo,:);
            for n = 1:size(Subs,3)
                ThisSub = squeeze(Subs(:,:,n));
                Mean(n,1) = mean(mean(ThisSub(logical(triu(ThisSub)))));
            end
            %scatter(MatchingSubsFD(FDIdx),Mean);lsline
            a = LinearModel.fit(MatchingSubsFD(FDIdx),Mean);
            NetworkFCFDCorr(net,nettwo) = corr(MatchingSubsFD(FDIdx),Mean);
            NetworkFCFDpvalue(net,nettwo) = table2array(a.Coefficients(2,4));
        else
            Subs = MatchingCorrMats(NetworkIds==net,NetworkIds==nettwo,:);
            for n = 1:size(Subs,3)
                ThisSub = squeeze(Subs(:,:,n));
                Mean(n,1) = mean(mean(ThisSub));
            end
            %scatter(MatchingSubsFD(FDIdx),Mean);lsline
            a = LinearModel.fit(MatchingSubsFD(FDIdx),Mean);
            NetworkFCFDCorr(net,nettwo) = corr(MatchingSubsFD(FDIdx),Mean);
            NetworkFCFDpvalue(net,nettwo) = table2array(a.Coefficients(2,4));
        end
    end
end
figure;imagesc(NetworkFCFDCorr);colormap(jet);colorbar;caxis([-.2 .2]);
