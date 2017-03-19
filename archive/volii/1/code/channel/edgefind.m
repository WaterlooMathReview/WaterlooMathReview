%% Computes the neighbour triangles
%
% Finds the neighbouring triangles by brute force searching for the three
% (or two in the case of boundary) triangles and returns the entry for that
% cell. This could probably be done quick but it runs quick enough. 
%
%%
function N = edgefind(NumTri,t)
%sort the t values since order doesn't matter as we are only looking for
%triangles with two of the same labels
T=sort(t);
%set the matrix to zero
N=zeros(3,NumTri);


for j=1:NumTri
    for i=1:NumTri
        %check if the first entry is in the neighbouring cell
        if((T(2,j)==T(2,i) || T(2,j)==T(3,i) || T(2,j)==T(4,i)) && j~=i)
            %if the cell is in the neighbouring cell check whether the
            %other two points are the neighbour
            if((T(3,j)==T(2,i) || T(3,j)==T(3,i) || T(3,j)==T(4,i)) && j~=i)
                N(1,j)=i;
            end
            if((T(4,j)==T(2,i) || T(4,j)==T(3,i) || T(4,j)==T(4,i)) && j~=i)
                N(2,j)=i;
            end
        end
        %the first case does not consider if the second and third entry is
        %in the neighbouring cell
        if((T(3,j)==T(2,i) || T(3,j)==T(3,i) || T(3,j)==T(4,i)) && j~=i)
            if((T(4,j)==T(2,i) || T(4,j)==T(3,i) || T(4,j)==T(4,i)) && j~=i)
                N(3,j)=i;
            end
        end        
    end
end