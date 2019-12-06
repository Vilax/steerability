% empirically testing facts about hilbert matrices

for k=2:10
    N = 2*k
    a = invhilb(N);
    [he, ho, heo, hoe] = separateMatrixParity(a);
    ge=sum(he(:));
    go=sum(ho(:));
    geo=sum(heo(:));
    
    ge/geo
    geo/go
    ge/go
    assert(ge > abs(geo));
    assert(abs(geo) > go);
    for a = 1:k
        fprintf('row %d', a)
        sherow = sum(he(a,:));
        shorow = sum(ho(a,:));
        sheorow  = sum(heo(a,:));
        shoerow = sum(hoe(a,:));
        
        fprintf('odd-even entries over odd-odd entries')
        shoerow/shorow
        fprintf('even entries over even-odd entries')
        sherow/sheorow

    end
end

% tests formula for a^ against computed values of optimal alpha

k = 5;
N = 2*k;
p=2
[H, f, A, b, Aeq, beq] = getQuadProgParams(N);
alpha = quadprog(H, f, A, b, Aeq, beq) 
lambda = computeLagrangeCoefficients(H);

h = inv(H)
H = generalizedHilbertMatrix(N,p);
h = inv(H);
[he, ho, heo, hoe] = separateMatrixParity(h);
ge=sum(he(:));
go=sum(ho(:));
geo=sum(heo(:));
goe = sum(hoe(:));

oddrows = h([1:2:N],:);
evenrows = h([2:2:N],:);
lambdaScale = 1/(ge*go - geo^2)

B = zeros(2,N);
oddId = [1:2:N];
evenId = [2:2:N];
B(1,oddId) = 1;
B(2,evenId) = 1;

adjhe = he * lambda(2);
adjho = ho * lambda(1);
adjheo = heo * lambda(1);
adjhoe = hoe * lambda(2);

rowsums = sum(h,2)
sum(heo,2) ./ sum(he,2)
sum(ho,2) ./ sum(hoe,2)

rowRatio = sum(oddrows(:))/ sum(evenrows(:))
row = 2;
thisrow = h(row,:);
oddthis = sum(thisrow(1:2:N));
eventhis = sum(thisrow(2:2:N));
thisratio = oddthis/eventhis
lambda(2)/lambda(1)

% now testing 