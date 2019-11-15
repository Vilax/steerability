function B = makeBMatrix(N)
    B = zeros([N, 2]);
    for id = 1:2:N
        B(id, 1) = 1;
        B(id+1, 2) = 1;
    end
end

