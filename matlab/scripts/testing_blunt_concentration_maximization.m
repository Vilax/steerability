clear all; close all;
kArr = [3:6];
nRuns = numel(kArr);
concArr = zeros([nRuns,1]);
for ik = 1:nRuns
    
k = kArr(ik);    
N = 2*k;
H = getQuadProgObjMatrix(k)
f = zeros(k,1);
A = []
b = []
Aeq = ones(1,k);
beq = 1;
alpha = quadprog(H,f,A,b,Aeq,beq)
concentration =  alpha' * H * alpha
denominator = getIntegralMonomialMatrix(k);
denomConcentration = alpha' * denominator * alpha;
concArr(ik) = concentration / denomConcentration


% now plot it
joined = vertcat(alpha', zeros(size(alpha')));
joined = reshape(joined, [numel(joined), 1]);
joined(numel(joined)+1) = 0;
p = joined;
x = [0:0.02:1];
pvals = polyval(p, x);
figure; plot(x, pvals)

end

function H = getQuadProgObjMatrix(k)
    [rowNum, colNum] = meshgrid(k:-1:1);
    H = 1./(2*rowNum+2*colNum+2);
end

function H = getIntegralMonomialMatrix(k)
    [rowNum, colNum] = meshgrid(k:-1:1);
    H = 1./(2*rowNum+2*colNum+1);
end



