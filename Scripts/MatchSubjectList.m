function SubjectMatchIdx = MatchSubjectList(big_subject_list,small_subject_list)
%	Will match subject lists of ABCD 
%   Input:  big_subject_list = string containing path of a .txt file with a column
%           of subjects (whichever list is biggest)
%           small_subject_list = string containing path of a .txt file with
%           a column of subjects (should be equal size or smaller)
%   Output: SubjectMatchIdx = Binary vector of length big_subject_list
%            1: Matching between big and small list
%            0: No matching subject 


% Check out the initial ABCD batch 
ABCD_Subjects_Big = textread(big_subject_list,'%s');
ABCD_Subjects_Small = textread(small_subject_list,'%s');
%ABCD_Subjects_All_Pfactor = textread('/data/nil-bluearc/GMT/Scott/ABCD/SubList_pfactor_AllSubjects_4k.txt','%s');
%ABCD_Subjects_2kList = textread('/data/nil-bluearc/GMT/Scott/ABCD/ABCD_SubjectList_2k.txt','%s');

% Reformat lists - only want to include NDAR #,

SubListBig = [];
if any(ABCD_Subjects_Big{1}=='_') % Check input structure
    for s = 1:size(ABCD_Subjects_Big,1)
        ThisSub = ABCD_Subjects_Big{s};
        SubListBig(s,:) = ThisSub(1,[1:4 6:end]);
    end
elseif any(ABCD_Subjects_Big{1}=='-') % Check input structure
    for s = 1:size(ABCD_Subjects_Big,1)
        subject = ABCD_Subjects_Big{s};
        SubListBig(s,:) = subject(5:end);
    end
else
    error('Unknown subject format! ... aborting')
end

% Get index of each subject in big list
SubjectMatchIdx = zeros(size(ABCD_Subjects_Big,1),1);

if any(ABCD_Subjects_Small{1}=='_') % Check input structure
    for s = 1:size(ABCD_Subjects_Small,1)
        ThisSub = ABCD_Subjects_Small{s};
        subject = ThisSub(1,[1:4 6:end]);
        for t = 1:size(ABCD_Subjects_Big,1)
            if all(SubListBig(t,:) == subject)
                SubjectMatchIdx(t,1) = 1;
            end
        end
    end
elseif any(ABCD_Subjects_Small{1}=='-') % Check input structure
    for s = 1:size(ABCD_Subjects_Small,1)
        subject = ABCD_Subjects_Small{s};
        subject = subject(5:end);
        for t = 1:size(ABCD_Subjects_Big,1)
            if all(SubListBig(t,:) == subject)
                SubjectMatchIdx(t,1) = 1;
            end
        end
    end
else
    error('Unknown subject format! ... aborting')
end

disp(['You have  ' num2str(sum(SubjectMatchIdx)) ' matching subjects'])
end

