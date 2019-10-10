function boolArray = index2boolean(indexArray, N)
% INDEX2BOOLEAN     Converts list of desired indices to boolean list where
% 1 gives indicated indices
% boolArray=INDEX2BOOLEAN(indexArray) takes as input an array of integers
% and (function will round if necessary) and desired maximum absolute index
% (since output array is assumed to be indexed in the form -N, -(N-1),...,
% 0, 1,... N) and output a boolean array indicating presence of index in
% input

    % first verify positive integers
    M = 2*N+1; 
    assert(max(abs(indexArray(:))) <= N);
    
    id = indexArray+N+1;
    
    boolArray = zeros([M,1]);
    boolArray(id)=1;
end

