function result()
global chain;
for i = 1:length(chain)
    a(i) = chain(i).inserted_THF;
    b(i) = length(chain(i).inserted_chain_pos);
end
figure(1);
hist(a);
figure(2);
hist(b,0:10);