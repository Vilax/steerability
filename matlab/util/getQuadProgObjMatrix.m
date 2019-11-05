function H = getQuadProgObjMatrix(k)
% assumes constant term is NOT DESIRED. k is max degree

    [rowNum, colNum] = meshgrid(k:-1:1);
    H = 1 ./ (rowNum+colNum+1);
end

