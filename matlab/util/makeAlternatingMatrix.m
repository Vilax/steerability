function M = makeAlternatingMatrix(N)

    [rowNum, colNum] = meshgrid(1:N, 1:N);
    sign = (-1) .^ (rowNum + colNum);
    
    M = abs(randn(N)) .* sign;
end

