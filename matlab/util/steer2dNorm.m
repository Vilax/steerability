function  [f,u,bCos,theta] = steer2dNorm(Theta,N)
    %close all
   k = N/2+1;

   %Create normalized cosine basis
   [bCos,sqrtSin,theta,dtheta]=getBasis(k,100);
   %Get the Gram matrices
   G1=makeGramMatrix(bCos,theta,dtheta,Theta)
   G2=makeGramMatrix(bCos,theta,dtheta,pi/2)
   %Generalized eigenvector problem
    [v,d]=eigs(G1,G2);
    %Print gen eig values and 1st eigvect
    d
    u=v(:,1)
    %Print constraints and value of obj func
    constraint=u'*G2*u
    obj=(u'*G1*u)/(u'*G2*u)
    %Make the steerable function and plot it
    f=makeSteerableFunction(v(:,1),bCos);
    f=f/f(1);
    fmin=min(f)
    figure(1);
    hold on;
    plot(theta,f);
    
    %Save basis and theta
    nTerms=k;
    nPoints=100;
    sCoeff=u;
%    save('steer2dBasis','f','theta','nTerms','nPoints','sCoeff');
end

function [bCos,sqrtSin,theta,dtheta]=getBasis(n,nSamples)
    dtheta=pi/(2*nSamples);
    theta=[0:dtheta:pi/2]';
    bCos=zeros(numel(theta),n);
    cosTheta=cos(theta);
    sqrtSin=sqrt(sin(theta));
    for k=1:n,
        bCos(:,k)=(cosTheta.^(2*(k-1)));
    end
    bCos=bCos;
    bCosNorm=(sqrt(diag(bCos'*bCos*dtheta)))';
    bCosNorm=repmat(bCosNorm,[size(bCos,1),1]);
    bCos=bCos./bCosNorm;
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