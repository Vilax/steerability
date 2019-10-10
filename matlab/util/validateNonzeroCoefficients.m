function numNonzero = validateNonzeroCoefficients(nonzeroCoefficients, N)

    assert(numel(nonzeroCoefficients) == 2*N+1);
    firstHalf = nonzeroCoefficients(1:N);
    secondHalf = nonzeroCoefficients(N+2:end);
    assert(isequal(firstHalf(:), flipud(secondHalf(:))));
    
    numNonzero = sum(nonzeroCoefficients(:));
end

