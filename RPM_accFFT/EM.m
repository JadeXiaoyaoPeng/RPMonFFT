function [epsv ] = EM( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud,eps0,INV_C  )
Nall=N(1)*N(2);
eps= invoitrans( epsv, N );

sig=epstosig( eps,lamdad,mud );
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

feb=cell(3,1);
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
s=cell(3,1);
for i=1:3
    s{i}=sa{i} +s_o{i};
end

eps=sigtoe( s,N,INV_C );
epsv  = voitrans( eps, Nall );
end

