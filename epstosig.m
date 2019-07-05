% calculate stress for isotropic materials
function [ sig ] = epstosig( eps,lamdad,mud )
    sig=cell(3,1);
    l2m=lamdad+2*mud;
    sig{1}=l2m.*eps{1}+lamdad.*eps{2};
    sig{2}=l2m.*eps{2}+lamdad.*eps{1};
    sig{3}=2*mud.*eps{3};
end

