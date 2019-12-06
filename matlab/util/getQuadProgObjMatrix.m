function H = getQuadProgObjMatrix(k)
% assumes constant term is NOT DESIRED. k is max degree

    [rowNum, colNum] = meshgrid(1:k);
    H = 1 ./ (rowNum+colNum+1);
end

