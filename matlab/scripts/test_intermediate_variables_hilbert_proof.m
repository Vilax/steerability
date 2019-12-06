p = 2; 
N = 8;

LHS_partial = 0;
for j =1:N
    summand1 = nchoosek(N+p+j-1, N);
    summand2 = nchoosek(N, j);
    summand = j*summand1 * summand2;
    LHS_partial = LHS_partial +  summand;
end
LHS = abs(LHS_partial) / (N*(p+N));

RHS = 0;

k=1
for j=1:N
    summand1 = nchoosek(N+p+j-1, N);
    summand2 = nchoosek(N, j);
    summand3 = j/(p+j+k-1);

    summand = summand1 * summand2 * summand3;
%     num1 = factorial(N+p+j-1);
%     num2 = factorial(p+k+j-2);
%     num3 = j * factorial(N);
%     den1 = factorial(N-j);
%     den2 = factorial(p+j+k-1);
%     den3 = factorial(p+j-1);
%     den4 = factorial(j-1);
%     den5 = j * factorial(N);
%     num = num1 * num2 * num3;
%     den = den1 * den2 * den3 * den4 * den5;
%     summand = num/den;
    RHS = RHS + summand;
end
RHS = abs(RHS)
rhs2 = 0;
for j=1:N
    scale = p+k+j-1;
    binomial1 = nchoosek(N+p+k-1,N-j);
    binomial2 = nchoosek(N+p+j-1,N-k);
    binomial3 = nchoosek(p+k+j-2,k-1);
    binomial4 = nchoosek(p+k+j-2,j-1);
    summand = scale * binomial1 * binomial2 * binomial3*binomial4;
    rhs2 = rhs2 + summand;
end
divisor1=k;
divisor2 = nchoosek(N+p+k-1,N);
divisor3 = nchoosek(N,k);
divisor = divisor1 * divisor2 * divisor3;
rhs2 = rhs2 / divisor
quotient = nchoosek(p+j+k-2, p+j-1)/nchoosek(p+j+k-1,j)
LHS
RHS
LHS/RHS
LHS < RHS
LHS < rhs2