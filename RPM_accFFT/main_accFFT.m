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
p1=[100]; %contrast K  5*10^4 10^5  5*10^5 ,10^(-5),10^(-6),10^(-1),10^(-2),10^(-3),10^(-4),10^(-5),10^(-6),10^(-7),10^(-8)
it_time=zeros(size(p1,2),1);
it_num=zeros(size(p1,2),1);
for b=1:size(p1,2)
        P1=p1(b);
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sig=epstosig( eps,lamdad,mud );
        fsig=cell(3,1);
        feb=cell(3,1);
        s=cell(3,1);
        
        ite_cont=0;       %iteration count
        conv_test=1;    %convergence test parameter
        tic                       %time
        while conv_test>10^(-4)
            
            
            sig_o=epstosig_o( eps,lamda0,mu0,N );
            sa=cell(3,1);
            sb=cell(3,1);
            for i=1: 3
                sa{i}=sig{i} - sig_o{i};
                sb{i}=2*sa{i};
            end
            
            fsb=cell(3,1);
            for i=1: 3
                fsb{i}=fft2(sb{i});
            end
            
            dfeb=greenop(  ks1,ks2,N,lamda0,mu0,fsb );
            for i=1: 3
                feb{i}=-dfeb{i};
            end
            
            eb=cell(3,1);
            for i=1: 3
                eb{i}=ifft2(feb{i}); % invese FFT of feps
                eb{i}=real(eb{i})+2*eps0{i};
            end
            
            s_o=epstosig_o( eb,lamda0,mu0,N );
            for i=1:3
                s{i}=sa{i} +s_o{i};
            end
            
            eps=sigtoe( s,N,INV_C );
            sig=epstosig( eps,lamdad,mud );
            
            %%%%%%%%%%%%%%%%%%%%%%
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
            ite_cont=ite_cont+1;
            conv_rec(ite_cont)=conv_test;
        end
        it_time(b)=toc;
        it_num(b)=ite_cont;
        
end
%% plot

% figure
% plot(xd1,(eps{3}(:,N(2)/2)))

