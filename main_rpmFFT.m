clear
clc
%INPUT discete coords
T=[1,1];              %period of the unit cell [T1,T2] 
N=[2^7,2^7];          %grids [N1,N2]
P2=10;                %parameter n_max
Nall=N(1)*N(2);
h=T./N;
fh=1./T;
xd1= 0:h(1):(N(1)-1)*h(1);
xd2= 0:h(2):(N(2)-1)*h(2);
[x1,x2]=meshgrid(xd1,xd2);
% frequency
ks1= 2*pi*fh(1) * [0: N(1)/2-1 , 0,   -(N(1)/2-1) :-1];  
ks2= 2*pi*fh(2) * [0: N(2)/2-1 , 0,   -(N(2)/2-1) :-1];
%%  dicrete Y & nu----lame's
% expression for Young's and Possion ratio Or input microstructure image
Yd=zeros(N);
nud=zeros(N);
    for i=1:N(1) 
        for j=1:N(2)
            if   (xd1(i)-0.5)^2  +  (xd2(j)-0.5)^2  >(1/32)^2
                Yd(i,j)=1;
                nud(i,j)=0.25;
            else
                Yd(i,j)=20 ;
                nud(i,j)=0.25;
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamdad=zeros(N);
mud=zeros(N);
for i=1:N(1)
    for j=1:N(2)
        lamdad(i,j)=Yd(i,j)*nud(i,j)/(1+nud(i,j))/(1-2*nud(i,j));
        mud(i,j)=0.5*Yd(i,j)/(1+nud(i,j));
    end
end
% reference lame's
lamda0=0.5*( min(min(lamdad))+max(max(lamdad)) );
mu0=0.5*( min(min(mud))+max(max(mud)) );
%% INPUT average strain
%prescribed strain 
eps00=[0;  0; 0.005] ;  % Voigt notation eps00(1)=e_11     eps00(2)=e_22      eps00(3)=e_12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps0=cell(3,1);
for i=1: 3
    eps0{i}=eps00(i)*ones(N);  % Voigt notation eps{1}=e_11     eps{2}=e_22      eps{3}=e_12
end
eps=eps0;
%% iteration
% initialization
epsv=voitrans( eps, Nall ); % convert to a vector, preparing for RPM
dlt=10^(-8);                % delta -- finite difference 
num=0;                      % # of iteration during the current unstable eigenspace
n_max=P2;
Z=zeros(4*Nall,0);
H=Z'*eval_derivative( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,Z,dlt ); % Jacobian
I=eye(size(H,2));
inv_I_H=(I-H)\I;
nepsv=MS( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud ); % FFT based fixed-point iteration
%%%%%%%%%%%%%%%%%%%%
ite_cont=0;     %iteration count
conv_test=1;    %convergence test 
tic             %time
while conv_test>10^(-4)
    z=Z'*epsv;
    zeta=Z'*nepsv;
    
    z=z+inv_I_H*(zeta-z);
    q=nepsv-Z*zeta;
    
    if num==n_max
        q3=q;
    elseif num==n_max-1
        q2=q;
    elseif num == n_max-2
        q1=q;
    end
    
    epsv=Z*z+q;
    nepsv=MS(ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud );
    ite_cont=ite_cont+1;
    num=num+1;

    if num>n_max
        D=[q3-q2,q2-q1];
        Z= incr_basis_size( Z,D,N );
        H=Z'*eval_derivative( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,Z,dlt );
        I=eye(size(H,2));
        inv_I_H=(I-H)\I;
        num=0;
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %convergence parameter   
    eps= invoitrans( nepsv, N );
    sig=epstosig( eps,lamdad,mud );
    fsig=cell(3,1);      % FFT of sig
    for i=1: 3
        fsig{i}=fft2(sig{i});
    end
    msig=[mean(mean(sig{1})), mean(mean(sig{3}));mean(mean(sig{3})), mean(mean(sig{2}))];
    con=zeros([N,2]);
    sum_div=0;
    for i=1:N(1)
        for j=1:N(2)
            con(i,j,1)=fsig{1}(i,j)*ks1(i) + fsig{3}(i,j)*ks2(j);
            con(i,j,2)=fsig{3}(i,j)*ks1(i) + fsig{2}(i,j)*ks2(j);
            sum_div=sum_div+norm(  [con(i,j,1);  con(i,j,2)]  )^2;
        end
    end
    conv_test=sqrt( sum_div / Nall) /norm(msig);
end
it_time=toc;
it_num=ite_cont;