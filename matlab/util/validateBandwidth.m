function L = validateBandwidth(L)
% perform validation for bandwidth/order L. First ensures it is even, and
% returns error if modified version is less than 0.
    if mod(L,2)~=0
        L=L-1;
    end
    assert(L>=0, 'L invalid');
end