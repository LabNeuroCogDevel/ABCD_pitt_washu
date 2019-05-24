function abcd_consensus_maker_knowncolors_alt_V3_marekmethod(template,ciftifile,outputfile,dconnfile,minsize)

NetworkMapping = Infomap_Assigning(template,ciftifile,outputfile);
NetworkMapping = NetworkMapping(:,[1:7 10:11]);
for v = 36000:36110
    for t = 1:size(NetworkMapping,2)
        if NetworkMapping(v,t) == 6
            NetworkMapping(v,t) = 10;
        end
    end
end
cifti_data = ft_read_cifti_mod(ciftifile);
Datalength = size(cifti_data.data,1);
% Color in all Row if any assignments

for i = 1:size(NetworkMapping,1)
    ThisRow = NetworkMapping(i,:);
    if max(ThisRow) == 4 && min(ThisRow) == 2 
        NetworkMapping(i,:) = 4;  
    end
    if min(ThisRow) == 0 && max(ThisRow) > 0
        idx = find(ThisRow~=0, 1, 'first');
        NetworkMapping(i,1:idx-1) = ThisRow(idx); 
    else
        NetworkMapping(i,:) = NetworkMapping(i,:);
    end

end

dotsloc = strfind(ciftifile,'.');
basename = ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_existingcolumnsfilledin'];
cifti_data.data = NetworkMapping;
%neighbors = cifti_neighbors(ciftifile);
if size(cifti_data.data,1) == 59412
    cifti_data.data(59413:Datalength,:) = 0;
    ft_write_cifti_mod(outname,cifti_data);
else
    ft_write_cifti_mod(outname,cifti_data);
end
set_cifti_powercolors([outname '.dtseries.nii']);


% Do consensus on filled in vertices
out_map_colored_single = zeros(size(NetworkMapping,1),1);
EmptyIdx = [];
for i = 1:size(NetworkMapping,1)
    ThisRow = NetworkMapping(i,:);
    if max(ThisRow) == 0
        EmptyIdx = [EmptyIdx; i];
    else
        if numel(unique(ThisRow))>1
            UniqueSet = unique(ThisRow)';
            %if any(UniqueSet==8)
                %out_map_colored_single(i,1) = 8;
            %else
            Output = [UniqueSet,histc(ThisRow,UniqueSet)'];
            [c d] = max(Output(:,2));
            idx = find(Output(:,2) == c);
            if length(idx) > 1
                out_map_colored_single(i,1) = max(Output(idx,1));
            else
                out_map_colored_single(i,1) = Output(d,1);
            end
            %end
        else
            out_map_colored_single(i,1) = ThisRow(1);
        end
    end

end

%OldInfoMap = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V2/infomap/' subject '/' subject '_rawassn_minsize400_regularized_recolored_cleaned.dscalar.nii']);
%ut_map_colored_single(OldInfoMap.data(1:59412) == 13) = 13;
cifti_data.data = out_map_colored_single;
dotsloc = strfind(ciftifile,'.');
basename = ciftifile(1:(dotsloc(end-1)-1));
outname = [basename 'semiclean_recolored'];
%neighbors = cifti_neighbors(ciftifile);
if size(cifti_data.data,1) == 59412
    cifti_data.data(59413:Datalength,:) = 0;
    ft_write_cifti_mod(outname,cifti_data);
else
    ft_write_cifti_mod(outname,cifti_data);
end
set_cifti_powercolors([outname '.dtseries.nii']);

% Clean up empty vertices
%Load dconn and distances
prevstring = [];
string = ['Loading correlations and distances'];
fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
prevstring = string;

dconn = ft_read_cifti_mod(dconnfile);
dconn=dconn.data(1:59412,1:59412);
dconn(dconn==0)=nan;
groupdata = ft_read_cifti_mod(template);
groupdata = groupdata.data(1:59412,1);
% Newinfomap = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/Files_forFigures/' subject '_rawassn_minsize400_regularized_recolored_cleaned.dscalar.nii']);
% Newinfomap = Newinfomap.data(1:59412);
EmptyIdx = find(out_map_colored_single == 0);
for i = 1:length(EmptyIdx)
%     out_map_colored_single(EmptyIdx(i),1) =Newinfomap(EmptyIdx(i),1);
    vert = EmptyIdx(i);
    Corrs = zeros(17,1);
    for j=1:17
        verts = j;%identify cortex
        Idx = groupdata == verts;
        Corrs(j,1) = nanmean(dconn(Idx,vert));
    end
    [a b] = nanmax(Corrs);
    out_map_colored_single(vert,1) = b;
end

%out_map_colored_single(Newinfomap == 13) = 13;
cifti_data.data = out_map_colored_single;
dotsloc = strfind(ciftifile,'.');
%basename = ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_filled_recolored_single'];
cifti_data.time = 1;
if size(cifti_data.data,1) == 59412
    cifti_data.data(59413:Datalength,1) = 0;
    ft_write_cifti_mod(outname,cifti_data);
else 
    ft_write_cifti_mod(outname,cifti_data);
end
set_cifti_powercolors([outname '.dtseries.nii']);


%% Clean up tine pieces

% Clean up tiny pieces
if logical(minsize)
    
    string = 'Cleaning maps';
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    neighbors = cifti_neighbors(ciftifile);
    
    for col = 1:size(out_map_colored_single,2)
        colmap = out_map_colored_single;
        
    
        allcolors = unique(colmap); 
        allcolors(allcolors<=0) = [];
    
    
    for color = allcolors(:)'
        clusteredmetric = zeros(size(colmap));
        thiscolorverts = find(colmap==color);
        for vertex = thiscolorverts'
            
            %find the neighbors of this vertex
            vertexneighbors = neighbors(vertex,:);
            vertexneighbors(isnan(vertexneighbors)) = [];
            
            %find which of those neighbors also pass the thresholds
            vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
            
            %find if those neighbors have already been assigned different cluster values
            uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
            uniqueneighborvals(uniqueneighborvals==0) = [];
            
            %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
            if isempty(uniqueneighborvals)
                clusteredmetric(vertexneighbors_thiscolor) = vertex;
                %if there is only one previous cluster identifier present, make all the neighbors that value
            elseif length(uniqueneighborvals)==1
                clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
                %if there are multiple cluster identifier values in the neighborhood, merge them into one
            else
                for valuenum = 2:length(uniqueneighborvals)
                    clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                end
            end
        end
        uniqueclustervals = unique(clusteredmetric);
        uniqueclustervals(uniqueclustervals==0) = [];
        
        for clusternum = uniqueclustervals'
            if nnz(clusteredmetric==clusternum) < minsize
                neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
                neighborverts(isnan(neighborverts)) = [];
                borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                borderverts(colmap(borderverts)<1) = [];
                mode_neighborval = mode(colmap(borderverts));
                colmap(clusteredmetric==clusternum) = mode_neighborval;
            end
        end
    end
    
    out_map_colored_single_cleaned = colmap;
    
    end
    
end

%out_map_colored_single_cleaned(Newinfomap == 13) = 13;
cifti_data.data = out_map_colored_single_cleaned;
dotsloc = strfind(ciftifile,'.');
%basename = ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_filled_recolored_single_cleaned'];
cifti_data.time = 1;
if size(cifti_data.data,1) == 59412
    cifti_data.data(59413:Datalength,1) = 0;
    ft_write_cifti_mod(outname,cifti_data);
else 
    ft_write_cifti_mod(outname,cifti_data);
end
set_cifti_powercolors([outname '.dtseries.nii']);


