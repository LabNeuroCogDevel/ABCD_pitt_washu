function k = make_stability_curve_allvars(NetworkCorrs,netone,nettwo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ThisNetCorrVars = squeeze(NetworkCorrs(netone,nettwo,:,:,:));
colors = zeros(size(ThisNetCorrVars,3),3);
colors(1,:) = [0 0 0];
colors(2:11,1) = [.55:.05:1];
colors(12:21,3) = [.55:.05:1];

figure;

for var = 1:size(ThisNetCorrVars,3)
    ThisNetCorr = ThisNetCorrVars(:,:,var);
    Mean_netcorr = mean(ThisNetCorr,2);
    lineProps.col = {colors(var,:)};% 170/255
    stdcorr = std(ThisNetCorr').*2;
    k = mseb([1:size(ThisNetCorrVars,1)],Mean_netcorr,stdcorr,lineProps,1);
end

set(gca,'FontWeight','bold','FontSize',14)
xlim([1 size(ThisNetCorrVars,1)]); ylim([-.5 .5])
xticks([0 12:12:size(ThisNetCorrVars,1)]) % every 300 subjects; each bin = 25 subs
xticklabels({'25'; '300';'600';'900';'1200';'1500';'1800';'2100';'2400';'2700';'3000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Brain Behavior Correlation','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
saveas(gcf,['/data/nil-bluearc/GMT/Scott/ABCD/ManhattanProject/Figures/BB_corrs_Raw/Net' num2str(netone) 'Net' num2str(nettwo) '_bbcorrs'],'epsc')

