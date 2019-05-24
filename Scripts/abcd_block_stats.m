function [t,p,d] = abcd_block_stats(matrix,partitions,order)
% Calculates a mean within blocks of a matrix
%   Input: matrix = an n x n matrix; zeros or ones on the diagonal
%          partitions = vector containing start index of each block; should NOT have a 1 for first block 
%          order = resorted subjects, to match partitions 

if partitions(1) == 1
        error('Remove the 1 in your partition index!')
end

temp = ones(length(partitions)+1,1);
temp(2:end) = partitions;
partitions = temp;
partitions(end+1) = size(matrix,1);

% Reorder matrix to align with partitions
matrix = matrix(order,order);

% make diagonal nans 
matrix(matrix==0)=nan;
t = [];
p = [];
d = [];
for x = 1:length(partitions)-2
    for y = x+1:length(partitions)-1
        % on diagonal block tests
        %Net one
        Netonestartidx = partitions(x);
        Netoneendidx = partitions(x+1)-1;
        Thismat = matrix(Netonestartidx:Netoneendidx,Netonestartidx:Netoneendidx);
        uidx = find(triu(Thismat,1));
        Dataone = Thismat(uidx);
        M1 = mean(Dataone);
        SD1 = std(Dataone);

        %Net two
        Nettwostartidx = partitions(y);
        Nettwoendidx = partitions(y+1)-1;
        Thismat = matrix(Nettwostartidx:Nettwoendidx,Nettwostartidx:Nettwoendidx);
        uidx = find(triu(Thismat,1));
        Datatwo = Thismat(uidx);
        M2 = mean(Datatwo);
        SD2 = std(Datatwo);
        [h pval Ci stats] = ttest2(Dataone,Datatwo);
        p = [p;pval];
        Thist = stats.tstat;
        t = [t;Thist];

        %d 
        sdpooled = sqrt((SD1^2 + SD2^2) / 2);
        Thisd = abs((M1-M2)/sdpooled);
        d = [d;Thisd];
    
    end
end

% Back to square matrix comparing blocks 
p = squareform(p);
t = squareform(t);
d = squareform(d);

end
