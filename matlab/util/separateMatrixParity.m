function [he, ho, heo, hoe] = separateMatrixParity(M)

    isSymmetric = isequal(M, M');
    if ~isSymmetric
        fprintf('Warning: input is not symmetric.')
    end
    
    [nrows, ncols] = size(M); 
    evenRowIndex = [2:2:nrows];
    oddRowIndex = [1:2:nrows];
    evenColIndex = [2:2:ncols];
    oddColIndex = [1:2:ncols];
    
    he = M(evenRowIndex, evenColIndex);
    ho = M(oddRowIndex, oddColIndex);
    heo = M(evenRowIndex, oddColIndex);
    hoe = M(oddRowIndex, evenColIndex);
end

