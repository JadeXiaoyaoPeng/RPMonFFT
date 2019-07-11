function [ epsv ] = voitrans( eps, Nall )
    epsv=zeros(4*Nall, 1);
    for i=1:Nall
        epsv(4*i-3)=eps{1}(i);
        epsv(4*i-2)=eps{3}(i);
        epsv(4*i-1)=eps{3}(i);
        epsv(4*i)=eps{2}(i);
    end

end

