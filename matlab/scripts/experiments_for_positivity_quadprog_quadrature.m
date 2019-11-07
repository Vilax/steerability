
for k = 1:10
N = 5;
[H, f, A, b, Aeq, beq] = getQuadProgParams(N);
range1 = generateRandIncrSeq(N, 'ascend');
range2 = generateRandIncrSeq(N, 'descend');
range1 = randn([1,N]);
range2 = randn([1,N]);
[x,y] = meshgrid(range1, range2);

x
y
H = x+y
% H = randn(N);
alpha = quadprog(H,f,A,b,Aeq,beq)
end
%theory: H only needs to be increasing in both row and col

function [H, f, A, b, Aeq, beq] = getQuadProgParams(N)
    H = getQuadProgObjMatrix(N);
    f = zeros(N,1);
    A = [];
    b = [];
    Aeq=zeros(2,N);
    for i=1:2:N
        Aeq(1,i)=1;
    end
    for i=2:2:N
        Aeq(2,i)=1;
    end
    beq = [1,-1];
    if mod(N,2) ~= 0
        beq = -1 * beq;
    end
end

function vals = generateRandIncrSeq(N,type)
    rangevals = randn([30,1]);
    b = max(rangevals);
    a = min(rangevals);
    if a == b
        a = a -1;
    end
    vals = sort((b-a)*rand(N,1)+a,type);
end