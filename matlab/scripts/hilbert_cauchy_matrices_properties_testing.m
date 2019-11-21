% empirically testing facts about hilbert matrices

for k=2:10
    N = 2*k
    a = invhilb(N);
    [he, ho, heo, hoe] = separateMatrixParity(a);
    she=sum(he(:));
    sho=sum(ho(:));
    sheo=sum(heo(:));
    
    she/sheo
    sheo/sho
    she/sho
    assert(she > abs(sheo));
    assert(abs(sheo) > sho);
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

k = 3;
N = 2*k;
[H, f, A, b, Aeq, beq] = getQuadProgParams(N);
alpha = quadprog(H, f, A, b, Aeq, beq) 
lambda = computeLagrangeCoefficients(H);