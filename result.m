function result()
global chain;
for i = 1:length(chain)
    a(i) = chain(i).inserted_THF;%THF num average over chains
    b(i) = length(chain(i).inserted_chain_pos);%Branching point average over chains
end
figure(1);
hist(a);
figure(2);
hist(b,0:10);
global weight
[Mn,Mw,PDI,weight] = calculate();
Mn
Mw
PDI
