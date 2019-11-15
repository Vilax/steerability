function p = makeEvenPolyOrigin(pEven)
% Given even coefficients of a polynomial that passes through the origin
% (omitting constant term), returns full polynomial. pEven should refer to
% coefficients of descending order

    nCoeff = numel(pEven);
    pEven = reshape(pEven, [1, nCoeff]);
    joined = vertcat(pEven, zeros(size(pEven)));
    joined = reshape(joined, [numel(joined(:)), 1]);
    joined(numel(joined)+1) = 0;
    p = joined;
end

