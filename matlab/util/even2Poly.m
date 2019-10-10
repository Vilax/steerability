function p = even2Poly(a)
% EVEN2POLY convert even polynomial to standard matlab polynomial format

    p = reshape([a'; zeros(size(a'))],[2*numel(a), 1]);
    p = p(1:end-1);
end

