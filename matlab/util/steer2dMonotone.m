function  [f,u,bCos,theta] = steer2dMonotone(Theta,N, nonzeroBool)
    % N must be maximum bandwidth allowed, nonzeroBool should have length
    % N+1. If nonzeroBool not provided, will take every other power to be
    % included, starting with the MAXIMUM POWER (fills in every other entry
    % with 1, starting with the rightmost entry). This function is modified
    % from steerd2dGeneral to sacrifice concentration for  montone
    % decreasing
    if nargin < 3
        nonzeroBool = zeros([1, N+1]);
        id = N+1-[0:2:N;];
        nonzeroBool(id) = 1;
    end
    assert(numel(nonzeroBool) == N+1);
    %close all

   %Create normalized cosine basis
   [bCos,sqrtSin,theta,dtheta]=getBasis(N,400, nonzeroBool);
   %Get the Gram matrices
   G1=makeGramMatrix(bCos,theta,dtheta,Theta)
   G2=makeGramMatrix(bCos,theta,dtheta,pi/2)
   %Generalized eigenvector problem
    [v,d]=eigs(G1,G2);
    %Print gen eig values and 1st eigvect
    d
    u=v(:,1)
    u = makeMonotone(u)
    %Print constraints and value of obj func
    constraint=u'*G2*u
    obj=(u'*G1*u)/(u'*G2*u)
    %Make the steerable function and plot it
    f=makeSteerableFunction(u,bCos);
    f=f/f(1);
    fmin=min(f)
    figure(1);
    hold on;
    plot(theta,f);
    
    %Save basis and theta
%     nTerms=k;
%     nPoints=100;
%     sCoeff=u;
% %   save('steer2dBasis','f','theta','nTerms','nPoints','sCoeff');
end

function v = makeMonotone(v)
    nterms = numel(v);
    for id = 1:nterms
        if v(id) > 0
            v(id) = -1 * v(id)
        end
    end
end

function [bCos,sqrtSin,theta,dtheta]=getBasis(N,nSamples, nonzeroBool)
    dtheta=pi/(nSamples);
    theta=[0:dtheta:pi]';
    bCos=zeros(numel(theta),N+1);
    cosTheta=cos(theta);
    sqrtSin=sqrt(sin(theta));
    for j=1:N+1,
        bCos(:,j)=(cosTheta.^(j-1));
    end
    bCos=bCos;
    bCosNorm=(sqrt(diag(bCos'*bCos*dtheta)))';
    bCosNorm=repmat(bCosNorm,[size(bCos,1),1]);
    bCos=bCos./bCosNorm;
    bCos = bCos(:,nonzeroBool==1);
end


function [f]=makeSteerableFunction(v,bCos)
    f=zeros(size(bCos,1),1);
    for k=1:numel(v),
        f=f+v(k)*bCos(:,k);
    end
end

function G=makeGramMatrix(bCos,theta,dtheta,thetaE)
    G=0;
    index=find(theta<=thetaE);
    bCos=bCos(index,:);
    G=bCos'*bCos*dtheta;
end