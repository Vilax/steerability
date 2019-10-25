function [f,u,bCos,phi]=  steer3dGeneral(Theta, N, nonzeroBool)
    
    % L is maximum bandwidth allowed, nonzeroBool indicates nonzero
    % coefficient to 1, z, z^2, z^3,...,z^L in that order
    if nargin < 3
        nonzeroBool = zeros([1, N+1]);
        id = N+1-[0:2:N;];
        nonzeroBool(id) = 1;
    end
    assert(numel(nonzeroBool) == N+1);
       
    
   %Create normalized cosine basis
   [bCos,sqrtSin,phi,dtheta]=getBasis(N, 400, nonzeroBool);
   %Get the Gram matrices
   G1=makeGramMatrix(bCos,sqrtSin,phi,dtheta,Theta)
   G2=makeGramMatrix(bCos,sqrtSin,phi,dtheta,pi/2)
   %Generalized eigenvector problem
    [v,d]=eigs(G1,G2);
    %Print gen eig values and 1st eigvect
    d
    u=v(:,1);
    u = u / norm(u);
    if u(1) < 0
        u = -1 * u;
    end
    u
    %Print constraints and value of obj func
   constraint=u'*G2*u
   obj=(u'*G1*u)/(u'*G2*u)
    %Make the steerable function and plot it
    f=makeSteerableFunction(v(:,1),bCos);
    fmin=min(f/f(1))
    figure(1);
    hold on;
    plot(phi,f/f(1));
    f = f/f(1);
end

function [bCos,sqrtSin,theta,dtheta]=getBasis(N, nSamples, nonzeroBool)
    dtheta=pi/(nSamples-1);
    theta=[0:dtheta:pi]';
    bCos=zeros(numel(theta),N);
    cosTheta=cos(theta);
    sqrtSin=sqrt(sin(theta));
    for j=1:N+1
        bCos(:,j)=(cosTheta.^(j-1));
    end
    bCosSin=bCos.*repmat(sqrtSin,[1 size(bCos,2)]);
    bCosNorm=(sqrt(diag(bCosSin'*bCosSin*dtheta)))';
    bCosNorm=repmat(bCosNorm,[size(bCos,1),1]);
    bCos=bCos./bCosNorm;
    bCos = bCos(:,nonzeroBool == 1);
end


function [f]=makeSteerableFunction(v,bCos)
    f=zeros(size(bCos,1),1);
    for k=1:numel(v),
        f=f+v(k)*bCos(:,k);
    end
end

function G=makeGramMatrix(bCos,sqrtSin,theta,dtheta,thetaE)
    index=find(theta<=thetaE);
    bCos=bCos(index,:);
    sqrtSin=sqrtSin(index);
    bCos=bCos.*repmat(sqrtSin,[1 size(bCos,2)]);
    G=bCos'*bCos*dtheta;
end