function [ sig_o ] = epstosig_o( eps,lamda0,mu0, N )
sig_o=cell(3,1);
lamda0m= lamda0 * ones(N);
mu0m= mu0 * ones(N);
l2m=lamda0m+2*mu0m;
            sig_o{1}=l2m.*eps{1}+lamda0m.*eps{2};
            sig_o{2}=l2m.*eps{2}+lamda0m.*eps{1};
            sig_o{3}=2*mu0m.*eps{3}; 
end

