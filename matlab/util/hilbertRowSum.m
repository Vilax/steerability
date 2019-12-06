function rowSum = hilbertRowSum(N,p,i)
    sign = (-1)^(N+i);
    binomial1 = nchoosek(N,i);
    binomial2 = nchoosek(N-1+p+i,N);
    rowSum = sign * i * binomial1 * binomial2;
%     div1 = factorial(i-1);
%     div2 = factorial(N-i);
%     prod = 1;
%     for k = 0:N-1
%         term = p+i+k;
%         prod = term * prod;
%     end
%     val = prod / (div1*div2);
%     rowSum = sign*val;
end

