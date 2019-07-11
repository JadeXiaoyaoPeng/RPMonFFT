function [ e ] = sigtoe( s,N,INV_C )
e=cell(3,1);
for i=1:N(1)
    for j=1:N(2)
        
        ee=INV_C{i,j}*[s{1}(i,j);s{2}(i,j);s{3}(i,j)];
        e{1}(i,j)= ee(1);
        e{2}(i,j)= ee(2);
        e{3}(i,j)= ee(3);
    end
end
end
