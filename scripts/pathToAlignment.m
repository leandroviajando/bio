function alignment= pathToAlignment(seq1,seq2,path,ScoringMatrix)

% Take the path through the alignment matrix and output the actual
% alignment to be used with the seqalignviewer
% inputs: seq1,seq2 are the two sequences that need to be aligned. 
%         seq1 corresponds to the first column of path.
%         ScoringMatrix is the substitution matrix to be used
%         (default=blosum50)

    if nargin<3
        error(message('IncorrectNumberOfArguments', mfilename));
    elseif nargin < 4
        ScoringMatrix=blosum50;
    end

    path=path-ones(size(path));
    % path=path(any(path>0,2),:);
    path=path(sum(path,2)>0,:);
    
    cur1=path(1,1);
    cur2=path(1,2);
    if size(path,1)>1
        for i=2:size(path,1)
            if path(i,1)==cur1
                path(i-1,1)=0;
                cur2=path(i,2);
            elseif path(i,2)==cur2
                path(i-1,2)=0;
                cur1=path(i,1);
            else
                cur1=path(i,1);
                cur2=path(i,2);
            end
        end
    end
    path = path(sum(path,2)>0,:);
    path = flipud(path);
    
    % Setting alignment size
    alignment = repmat(('- -')',1,size(path,1));

    % Adding sequence to alignment
    % alignment(1,path(:,1)>0) = seq1;
    % alignment(3,path(:,2)>0) = seq2;
    alignment(1, path(:,1)>0) = seq1(path(path(:,1)>0, 1));
    alignment(3, path(:,2)>0) = seq2(path(path(:,2)>0, 2));

    % Find locations where there are no gaps
    h=find(all(path>0,2));
    
    noGaps1=aa2int(alignment(1,h));
    noGaps2=aa2int(alignment(3,h));
    
    scoringMatrixSize=size(ScoringMatrix,1);
    
    % Erasing symbols that cannot be scored
    htodel=max([noGaps1;noGaps2])>scoringMatrixSize;
    h(htodel)=[];
    noGaps1(htodel)=[];
    noGaps2(htodel)=[];

    % Score pairs with no gap
    value = ScoringMatrix(sub2ind(size(ScoringMatrix),double(noGaps1),double(noGaps2)));

    % Insert symbols of the match string into the alignment
    alignment(2,h(value>=0)) = ':';
    alignment(2,h(noGaps1==noGaps2)) = '|';
end
