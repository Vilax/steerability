
close all;
x = [0:0.01:1];
for N =4:10
    monomial = zeros([N+1,1]);
    monomial(1) = 1; 
    monovals = polyval(monomial, x);
    
    [~, ~, u, ~, ~, ~] = makeQuadPoly(N, 500);
    u = u/sum(u);
    vals = polyval(u, x);
    
    figure(2); 
    plot(x, monovals);
    hold on; 
    plot(x, vals);
    legend;
end