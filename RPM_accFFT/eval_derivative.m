function [ FuZ ] = eval_derivative( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,Z,dlt,eps0,INV_C )
FuZ=zeros(4*N(1)*N(2),size(Z,2));
for i=1:size(Z,2)
    FuZ(:,i)=( EM( ks1,ks2,lamda0,mu0,epsv+Z(:,i)*dlt,N,lamdad,mud,eps0,INV_C ) ...
        - EM( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,eps0,INV_C ) )/dlt;
end
end

