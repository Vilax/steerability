function powers = dirCosPowers3D(N)
% DIRCOSPOWERS3D    All combinations of putting N balls in 3 bins

%     combinations = nchoosek(0:N, 3);
%     rowsum = sum(combinations,2);
%     powersOneIter = combinations(rowsum==N,:);
%     nrows = size(powersOneIter,1);
%     powers = [];
%     for irow = 1:nrows
%         powers = [powers; perms(powersOneIter(irow,:))];
%     end

    M=(N+2)*(N+1)/2;
    powers = zeros(M, 3);
    irow = 1;
    for a = N:-1:0
        resA = N - a;
        for b = resA:-1:0
            c = resA - b;
            powers(irow,:) = [a, b, c];
            irow = irow+1;
        end
    end
end

