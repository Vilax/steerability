% test if 'corresponding blocks of the cross-multiply' are always positive
% (and small) 
k = 2;
N = 2*k;
p = 2;
H = generalizedHilbertMatrix(N,p);
h = inv(H);

numPos = 0;
numNeg = 0;
valArray = zeros([lim^4,1]);
id = 0;

lim = N/2;
for i1 = 1:lim
    for j1 = 1:lim
        for i2 = 1:lim
            for j2 = 1:lim
                val = multiplyBlock(i1,j1,i2,j2,h);
                
                id = id+1;
                valArray(id) = val;
            end
        end
    end
end

function val = multiplyBlock(i1,j1,i2,j2,h)
    
    N = size(h,1);
    a = max([i1, j1, i2, j2]); 
    assert(a*2 <= N);
    
    oddi = 2*i1 - 1;
    oddj = 2*j1 - 1;
    eveni = 2*i2;
    evenj = 2*j2;
    
    val = h(oddi, oddj)* h(eveni,evenj) - h(oddi,evenj)*h(eveni,oddj);
end