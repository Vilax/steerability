function v = computeHilbertPair(N, Theta)
    nonzeroBool = zeros([1, N+1]);
    id = N+1-[0:2:N;];
    nonzeroBool(id) = 1;
    nonzeroBool(1) = 0;
    [f,u,~,theta] = steer2dGeneral(Theta, N, nonzeroBool);
    
    u = u/sum(u);
    u = flipud(u);
    u = zeros(size(u));
    u(1) = 1;
    m = numel(u);
     % for visualizing
    pad = zeros([m,1]);
    joined = vertcat(u',pad');
    polycoeff = reshape(joined, [2*m,1]);
    a = cos(theta); 
    if mod(N,2) == 0
        polycoeff(2*m+1) = 0;
    end
    ideq = zeros(m, N);
    for id = 1:m
        coordinate = 2*id-1;
        ideq(id, coordinate) = 1;
    end
    vals = polyval(polycoeff, a);
    figure; plot(a, vals);
    % end visualizations
    
    sumeqMat=zeros(2,N);
    for i=1:2:N
        sumeqMat(1,i)=1;
    end
    for i=2:2:N
        sumeqMat(2,i)=1;
    end
    sumeqVec = [1,-1]';
    if mod(N,2) ~= 0
        sumeqVec = -1 * sumeqVec;
    end    
    H = getQuadProgObjMatrix(N);
    f = zeros([N,1]);
    Aeq = vertcat(sumeqMat, ideq);
    beq = vertcat(sumeqVec, u);
%     A = -eye(N);
%     for id = [1:2:N]
%         A(id,id) = 1;
%     end
%     b = zeros([1,N]);
    A= []; b= [];
    
    ub  = []; lb = [];
    small_constant = 1e-6;
    SMALL_CONSTANT = 1e-6;
    ub = ones([N,1]) * SMALL_CONSTANT;
    for id=[1:2:N];
        ub(id) = Inf;
    end
    v = quadprog(H,f,A,b,Aeq,beq,lb,ub)
end

