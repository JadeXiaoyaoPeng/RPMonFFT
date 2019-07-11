function [epsv ] = MS( ks1,ks2,lamda0,mu0,epsv,N,lamdad,mud )
Nall=N(1)*N(2);
eps= invoitrans( epsv, N );
feps=cell(3,1);
for i=1: 3
    feps{i}=fft2(eps{i});
end

sig=epstosig( eps,lamdad,mud );
fsig=cell(3,1);      % FFT of sig
for i=1: 3
    fsig{i}=fft2(sig{i});
end

dfeps=greenop(  ks1,ks2,N,lamda0,mu0,fsig );
    for i=1: 3
        feps{i}=fft2(eps{i})-dfeps{i};
    end
    
    for i=1: 3
    eps{i}=ifft2(feps{i}); % invese FFT of feps 
    eps{i}=real(eps{i});
    end 
    epsv  = voitrans( eps, Nall );
end

