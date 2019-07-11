function [ eps ] = invoitrans( epsv, N )
    eps=cell(3,1);
    for i=1:3
        eps{i}=zeros(N);
    end
    for i=1:N(1)*N(2)
        eps{1}(i)=epsv(4*i-3);
        eps{3}(i)=epsv(4*i-2);
        eps{2}(i)=epsv(4*i);
    end
end