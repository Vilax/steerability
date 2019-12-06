function entry = hilbertEntryVal(N,p,i,j)

    sign = (-1)^(i+j);
    scale = p+i+j-1;
    binomial1 = nchoosek(N+p+i-1,N-j);
    binomial2 = nchoosek(N+p+j-1,N-i);
    binomial3 = nchoosek(p+i+j-2,i-1);
    binomial4 = nchoosek(p+i+j-2,j-1);
    entry = sign * scale * binomial1 * binomial2 * binomial3*binomial4;
    
%     prod = 1;
%     for k = 0:N-1
%         val = (p+i+k)*(p+j+k);
%         prod = val * prod;
%     end
%     quo1 = factorial(i-1);
%     quo2 = factorial(N-i);
%     quo3 = factorial(j-1);
%     quo4 = factorial(N-j);
%     quo = quo1 * quo2 * quo3 * quo4;
%     entry = sign * prod /(scale * quo);
end

