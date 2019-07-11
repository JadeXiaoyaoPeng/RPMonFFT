function [ dfeps ] = greenop( ks1,ks2,N,lamda0,mu0,fsig )
    dfeps=cell(3,1);    
    for i=1:N(1)
        for j=1:N(2)
            cks=[ks1(i);ks2(j)];
            cur_fsig=[fsig{1}(i,j);  fsig{2}(i,j);  fsig{3}(i,j)];
            cur_dfeps=zeros(3,1);
            
            A=1/(4*mu0*norm(cks)^2);
            if isinf(A) == 1
                A=0;
            end
            B=-(lamda0+mu0) / (mu0*(lamda0+2*mu0)*norm(cks)^4);
            if isinf(B) == 1
                B=0;
            end
            C=B*( cks(1)^2*cur_fsig(1) + 2*cks(1)*cks(2)*cur_fsig(3) + cks(2)^2*cur_fsig(2) );
            
            cur_dfeps(1)=4*A*( cks(1)*cks(1)*cur_fsig(1)+cks(1)*cks(2)*cur_fsig(3)) +C*cks(1)*cks(1);
            cur_dfeps(2)=4*A*( cks(1)*cks(2)*cur_fsig(3)+cks(2)*cks(2)*cur_fsig(2)) +C*cks(2)*cks(2);
            cur_dfeps(3)=2*A*( cks(1)*cks(2)*(cur_fsig(1)+cur_fsig(2))  +  (cks(1)*cks(1)+cks(2)*cks(2))*cur_fsig(3)  ) +C*cks(1)*cks(2);
  
            dfeps{1}(i,j)=cur_dfeps(1);
            dfeps{2}(i,j)=cur_dfeps(2);
            dfeps{3}(i,j)=cur_dfeps(3);
        end
    end


end

