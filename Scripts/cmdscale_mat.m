function [Y E h] = cmdscale_mat(mat,groups,distmet,varargin)
% This function performs multi-dimensional scaling on input matrices and
% displays a 2-dimensional plot of the relative positions of the input 
% matrices colored by group. MDS is performed using euclidean distance.
% Inputs:
% mat - node x node x subject array of matrices from all subjects
% groups - identifies which subject belongs to which group
% distmet - distance metric for comparing matrices, e.g. 'euclidean','hamming'
% varargin - group names can be specified by a third input, e.g. 
% {'group1';'group2'} and will be displayed in a legend.
% Outputs:
% Y - configuration matrix
% E - eigenvalues of Y*Y' 
% 
% TOL, 09/14
numnodes = size(mat,1);
numgroups = length(unique(groups));
colors = {'k';'y';'g';'r';'b';'c';'m';[.5 .5 .5];[.5 .9 .5];[1 .5 0]};
%colors = distinguishable_colors(numgroups);
if nargin > 3
    groupnames = varargin{1};
end

mask = ones(numnodes);
mask = triu(mask,1);
for s = 1:size(mat,3)
        temp = mat(:,:,s);
        mat_col(:,s) = temp(logical(mask));
end

% Multi-dimensional scaling
D = pdist(double(mat_col'),distmet);
[Y E] = cmdscale(D);

% Display result
h = figure('Color','white')    
hold
for c = 1:numgroups
    inds = groups==c;
    %plot(Y(inds,1),Y(inds,2),'.','Color',colors(c,:),'MarkerSize',20)
    plot3(Y(inds,1),Y(inds,2),Y(inds,3),'.','Color',colors{c},'MarkerSize',20)
end
set(gca,'FontWeight','bold','FontSize',14)
if nargin > 3
    legend(groupnames,'FontWeight','bold','FontSize',14)
end
xlabel('MDS 1','FontWeight','bold','FontSize',14)
ylabel('MDS 2','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
