function NetworkCorrs = abcd_corrmat_reliability_curve(rmats,binsize,iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

subnum = size(rmats,3);

NetworkCorrs = zeros(length([binsize:binsize:subnum]),iter);

Meannet = mean(rmats,3);
uidx = find(triu(Meannet,1));
vector = Meannet(uidx);

for x = [5 10 15 20 binsize:binsize:subnum]
    disp(['On sub size ' num2str(x)])
    for i = 1:iter
        % permute sib indices 
        p = randperm(size(rmats,3));
        idx = p(1:x);
        % Half Two
        
        TheseMats = mean(rmats(:,:,idx),3);
        TheseMats = TheseMats(uidx);
        NetworkCorrs(1,i) = corr(TheseMats,vector); %x/binsize
    end
    disp(['Mean is ' num2str(mean(NetworkCorrs(1,:)))])
    
end




end

