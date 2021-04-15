function [isSingular VectorsF]  = isMatrixSingular2F(Vectors)

%Vectors = EbasisMin(:,1:17)';



A = Vectors';

Nrow=size(A,1);
Ncol=size(A,2);

if any(sum(A,2) == 0)
    isZeroRow = 1;
else
    isZeroRow = 0;
end

%check if there are same vectors with single edges (?) -> not useful for cycles ? -> skip?
AsingleEdges = A(:,sum(full(A),1)==1);
AsingleEdges = unique(AsingleEdges','rows')';
NsingleEdges = size(AsingleEdges,2);
if NsingleEdges<Nrow
    isNotSingularForSure=0;
else if NsingleEdges==Nrow
        isNotSingularForSure=1;
    else
        disp('Impossible!');
        isNotSingularForSure=2;
    end
end

%sort to upper triangular like matrix
jrow=1;
while jrow <= Nrow && isZeroRow==0 && isNotSingularForSure==0
    
    [dum i_max]  = max(abs(A(jrow:Nrow, jrow)));
    i_max = i_max+jrow-1;
    A([i_max jrow], :) = A([jrow i_max ], :);
    
    %make 'upper triangular'
    for i = jrow + 1 : Nrow
        
        if A(i, jrow)~=0
            %B  = xor(A(i, :) , A(jrow, :));
            A(i, :) = xor(A(i, :) , A(jrow, :));
        else
            %B  =A(i, :);
        end
        
        %A(i, :) = B;
        
    end
    
    if any(sum(A,2) == 0)
        isZeroRow = 1;
    end
    
    jrow = jrow+1;
end

% check if matrix is diagonal
AsingleEdges = A(:,sum(full(A),1)==1);
AsingleEdges = unique(AsingleEdges','rows')';
NsingleEdges = size(AsingleEdges,2);
if NsingleEdges<Nrow
    isNotSingularForSure=0;
else if NsingleEdges==Nrow
        isNotSingularForSure=1;
    else
        disp('Impossible!');
        isNotSingularForSure=2;
    end
end




% back substitution
jcount = 1;
while isZeroRow==0 && jcount<500000 && isNotSingularForSure==0

    % find cols that have non single entry
    KcolSingleOnes = find(sum(full(A),1)==1);
    %[dum KcolFirstOnes] = max(full(A'));
    [dum KcolFirstOnes] = max(A,[],2);
    KcolFirstOnes = full(KcolFirstOnes);
    KcolFirstOnesNotSingle = my_setdiff(KcolFirstOnes, KcolSingleOnes);
    
    if ~isempty(KcolFirstOnesNotSingle)
        columnLastOne =max(KcolFirstOnesNotSingle);
        
        rowLastOne = find(KcolFirstOnes==columnLastOne);
        
        rowsToConsid = find(A(:,columnLastOne));
        
        rowsForSubtraction = rowsToConsid(rowsToConsid~=rowLastOne(1));
        
        for js=1:length(rowsForSubtraction)
                A(rowsForSubtraction(js),:) = xor(A(rowsForSubtraction(js),:),A(rowLastOne(1),:));
        end
 
    else
        isNotSingularForSure=1;
    end
    
    isZeroRow = any(sum(full(A),2)==0);
    jcount = jcount+1;
end

if jcount>400000
    disp(['Problem: too many total iter: ' int2str(jcount)])
end


if isZeroRow==1 && isNotSingularForSure==0
    isSingular = 1;
    %disp('I think Matrix is singular')
else if (isZeroRow==1 && isNotSingularForSure==1) || (isZeroRow==0 && isNotSingularForSure==0)
        isSingular = 2;
        disp('Problem')
    else if isZeroRow==0 && isNotSingularForSure==1
            isSingular = 0;
            %disp('Not singular')
        end
    end
end

VectorsF = A';
%
% A(3,:) = xor(A(3,:),A(4,:)); full(A)
% A(5,:) = xor(A(5,:),A(6,:)); full(A)
% A(6,:) = xor(A(6,:),A(7,:)); full(A)
% A(5,:) = xor(A(5,:),A(7,:)); full(A)
% A(7,:) = xor(A(7,:),A(10,:)); full(A)
% A(6,:) = xor(A(6,:),A(10,:)); full(A)
% A(5,:) = xor(A(5,:),A(10,:)); full(A)
% A(7,:) = xor(A(7,:),A(11,:)); full(A)
% A(6,:) = xor(A(6,:),A(11,:)); full(A)
% A(5,:) = xor(A(5,:),A(11,:)); full(A)
% A(4,:) = xor(A(4,:),A(13,:)); full(A)
% A(3,:) = xor(A(3,:),A(13,:)); full(A)
% A(15,:) = xor(A(15,:),A(16,:)); full(A)
% A(11,:) = xor(A(11,:),A(16,:)); full(A)
% A(7,:) = xor(A(7,:),A(16,:)); full(A)
% A(6,:) = xor(A(6,:),A(16,:)); full(A)
% A(5,:) = xor(A(5,:),A(16,:)); full(A)
%%% BACkup
% function isSingular = isMatrixSingularF(Vectors)
%
%  isSingular = 0;
% A=Vectors;
%
% Nrow=size(A,1);
% Ncol=size(A,2);
% for jrow = 1 : Nrow
%     %Find pivot for column k:
%     [dum i_max]  = max(abs(A(jrow:Nrow, jrow)));
%     i_max = i_max+jrow-1;
%     if A(i_max, jrow) == 0
%         disp(['Matrix is singular!'])
%         isSingular = 1;
%     end
%     %swap rows(jrow, i_max)
%     A([i_max jrow], :) = A([jrow i_max ], :);
%     %Do for all rows below pivot:
%     for i = jrow + 1 : Nrow
%
%         %Do for all remaining elements in current row:
%         B = zeros(1,Ncol);
%
%
%
%         B(jrow : Ncol)  = A(i, jrow : Ncol) - A(jrow, jrow : Ncol) * (A(i, jrow) / A(jrow, jrow));
%
%         if A(i, jrow)~=0
%             B(jrow : Ncol)  = xor(~~A(i, jrow : Ncol) , ~~A(jrow, jrow : Ncol));
%         else
%             B(jrow : Ncol)  =A(i, jrow : Ncol);
%         end
%
%
%         A(i, :) = B;
%         %A(i, jrow)  = 0;  %<--------is this needed
%
%     end
%
% end
