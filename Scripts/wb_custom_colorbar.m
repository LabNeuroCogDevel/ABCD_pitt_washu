function wb_custom_colorbar(cifti,scale,varargin)

% set up colorbar & bins
numbins=100;
blues=[(1:-1/(numbins-1):0)' (1:-1/(numbins-1):0)' ones(numbins,1)];%.^.5;
reds=[ones(numbins,1) (1:-1/(numbins-1):0)' (1:-1/(numbins-1):0)'];%.^.5;
colors=[reds; blues];
needed_colors=zeros(size(colors,1),1);

% read cifti
orig=ft_read_cifti_mod(cifti);

if (nargin>2) % percentile flag
    if scale(2)>0
        positives=orig.data(orig.data>0);
        positives=sort(positives);
        scale(2)=positives(round(length(positives)*scale(2)));
    end
    if scale(1)<0
        negatives=orig.data(orig.data<0);
        negatives=sort(negatives,'descend');
        scale(1)=negatives(round(length(negatives)*-scale(1)));
    end
end

% create label
label=orig;
label.data=zeros(size(label.data));
if scale(2)>0
    poslims=0:scale(2)/numbins:scale(2);
    for i=1:numbins
        label.data(orig.data>poslims(i) & orig.data<=poslims(i+1))=i;
        needed_colors(i)=sum(orig.data>poslims(i) & orig.data<=poslims(i+1));
    end
    label.data(orig.data>poslims(end))=numbins;
    needed_colors(numbins)=needed_colors(numbins)+sum(orig.data>poslims(end));
end
if scale(1)<0
    neglims=0:scale(1)/numbins:scale(1);
    for i=1:numbins
        label.data(orig.data<neglims(i) & orig.data>neglims(i+1))=i+numbins;
        needed_colors(numbins+i)=sum(orig.data<neglims(i) & orig.data>neglims(i+1));
    end
    label.data(orig.data<neglims(end))=numbins*2;
    needed_colors(numbins*2)=needed_colors(numbins*2)+sum(orig.data<neglims(end));
end
needed_colors=logical(needed_colors);
colors=colors(needed_colors,:);

[~,name]=fileparts(cifti);
[~,name]=fileparts(name);
ft_write_cifti_mod([name '_custom_colorbar.dtseries.nii'],label);
%make_cifti_label([name '_custom_colorbar.dtseries.nii'],colors);