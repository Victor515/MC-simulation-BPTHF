function [Mn,Mw,PDI] = result()
global chain;
for i = 1:length(chain)
    a(i) = chain(i).inserted_THF;%THF number of every chain
    b(i) = length(chain(i).inserted_chain_pos);%Number of Branching point of every chains
end
figure(1);
hist(a);
title('THF num distribution');
figure(2);
hist(b,0:10);
title('Branching point num distribution');
global weight
[Mn,Mw,PDI,weight] = calculate();
figure(3);
weight = weight(weight ~= 0);
hist(weight,10);
title('weight distribution');
length(weight)

