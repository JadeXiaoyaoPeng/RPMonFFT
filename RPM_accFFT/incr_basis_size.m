function [ Z ] = incr_basis_size( Z,D,N )
if size(Z,2)>=4*N(1)*N(2)
    return
end
[Q, R] = gs_m(D);
if abs(R(1))>1000*abs(R(4))
    Z=cat(2,Z,Q(:,1));
else
    Z=cat(2,Z,Q);
end
end

