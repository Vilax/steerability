% Finds even polynomials to minimize int_0^1 (1-x^2)(p1(x)^2+p2(x)^2) dx
% where p1 and p2 are both even polynomials

% we will represent the monomial integral matrix of (1-x^2) as a block
% identity matrix, and the coefficients of p1, p2 (even only) as a 2*N
% vector a 

k = 3;
N = 2*k; 

[rowNum, colNum] = meshgrid(k:-1:1, k:-1:1);
H1 = -1 ./ (2*rowNum + 2*colNum+3);
H2 = 1 ./ (2*rowNum+2*colNum+1);

subH = H2 - H1;

subZero = zeros(size(subH));

H = [subH, subZero; subZero, subH]; 

Aeq = ones(2, 2*k);
Aeq(1, k+1:end) = 0;
Aeq(2, 1:k) = 0;
beq = [1,1];

f = zeros(2*k,1);
A = [];
b = [];

alpha = quadprog(H, f, A, b, Aeq, beq) 

concentration =  alpha' * H * alpha
% denominator = getIntegralMonomialMatrix(k);
% denomConcentration = alpha' * denominator * alpha;

p1 = makeEvenPolyOrigin(alpha(1:k))
p2 = makeEvenPolyOrigin(alpha(k+1:end))

x = [0:0.02:1];
p1vals = polyval(p1, x);
p2vals = polyval(p2, x);

pvals = p1vals.^2 + p2vals.^2;
figure; plot(x, pvals);


function H = getIntegralMonomialMatrix(k)
    [rowNum, colNum] = meshgrid(k:-1:1);
    H = 1./(2*rowNum+2*colNum+1);
end


