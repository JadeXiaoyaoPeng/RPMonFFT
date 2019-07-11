clear
clc
%INPUT discete coords
T=[1,1];              %period of the unit cell [T1,T2] 
N=[2^8,2^8];              %grids [N1,N2]
Nall=N(1)*N(2);
h=T./N;
fh=1./T;
xd1= 0:h(1):(N(1)-1)*h(1);
xd2= 0:h(2):(N(2)-1)*h(2);
[x1,x2]=meshgrid(xd1,xd2);
% frequency
ks1= 2*pi*fh(1) * [0: N(1)/2-1 , 0,   -(N(1)/2-1) :-1];  
ks2= 2*pi*fh(2) * [0: N(2)/2-1 , 0,   -(N(2)/2-1) :-1];
%% change parameters
p1=[100]; %contrast K 5*10^2 10^3 5*10^3 10^4 5*10^4 10^5 5*10^5 ,10^(-5),10^(-6),10^(-1),10^(-2),10^(-3),10^(-4),10^(-5),10^(-6),10^(-7),10^(-8)
p2=[10];  %n_max
it_time=zeros(size(p1,2), size(p2,2));
it_num=zeros(size(p1,2), size(p2,2));
for b=1:size(p1,2)
    P1=p1(b);
    for a=1 : size(p2,2)
        P2=p2(a);
        %%  dicrete Y & nu----lame's
        % INPUT expression for Young's and Possion or specify Yd and nud
        Yd=zeros(N);
        nud=zeros(N);
        for i=1:N(1)
            for j=1:N(2)
                if   (xd1(i)-0.5)^2  +  (xd2(j)-0.5)^2  >(1/32)^2
                    Yd(i,j)=1;
                    nud(i,j)=0.25;
                else
                    Yd(i,j)=P1 ;
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
%         lamda0=0.5*( min(min(lamdad))+max(max(lamdad)) );
%         mu0=0.5*( min(min(mud))+max(max(mud)) );
        mu0=sqrt( min(min(mud))*max(max(mud)) );
        Km=lamdad+ 2/3*mud;
        K0=sqrt( min(min(Km))*max(max(Km)) );
        lamda0=K0-2/3*mu0;
        
        lamdam=lamdad+lamda0*ones(N);
        mum=mud+mu0*ones(N);
        
        INV_C=cell(N(1),N(2));
        for i=1:N(1)
            for j=1:N(2)
                C=[lamdam(i,j)+2*mum(i,j)       lamdam(i,j)       0
                    lamdam(i,j)     lamdam(i,j)+2*mum(i,j)    0
                    0                     0         2*mum(i,j)];
                INV_C{i,j}= inv(C);
            end
        end
        %% INPUT average strain
        %prescribed strain in the x1 direction
        eps00=[0;  0; 0.005] ;  % eps00(1)=e_11     eps00(2)=e_22      eps00(3)=e_12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eps0=cell(3,1);
        for i=1: 3
            eps0{i}=eps00(i)*ones(N);  % eps{1}=e_11     eps{2}=e_22      eps{3}=e_12
        end
        eps=eps0;
        %% iteration
        epsv=voitrans( eps, Nall );
        dlt=10^(-8);
        num=0;
        n_max=P2;
        Z=zeros(4*Nall,0);
        H=Z'*eval_derivative( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,Z,dlt );
        I=eye(size(H,2));
        inv_I_H=(I-H)\I;
        nepsv=EM( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,eps0,INV_C  );
        %%%%%%%%%%%%%%%%%%%%
        ite_cont=0;       %iteration count
        conv_test=1;    %convergence test parameter
        tic                       %time
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
            nepsv=EM( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,eps0,INV_C  );
            ite_cont=ite_cont+1;
            num=num+1;
            
            if num>n_max
                D=[q3-q2,q2-q1];
                Z= incr_basis_size( Z,D,N );
                H=Z'*eval_derivative( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,Z,dlt,eps0,INV_C );
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
            conv_rec(ite_cont)=conv_test;
        end
        it_time(b,a)=toc;
        it_num(b,a)=ite_cont;
        
    end
end
%% Plot

% figure
% plot(xd1,eps{3}(:,N(2)/2+1))
% eps3_0=eps{3}(:,N(2)/2+1);