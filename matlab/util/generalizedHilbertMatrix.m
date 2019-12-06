function H = generalizedHilbertMatrix(N, p)
% entry in row i column j is 1/(i+j-1+p). For the optimization of monomials
% problem, we are interested in p=2

    if nargin < 2
        p = 2;
    end
    [rowNum, colNum] = meshgrid(1:N);
    H = 1 ./ (rowNum+colNum-1+p);
end

