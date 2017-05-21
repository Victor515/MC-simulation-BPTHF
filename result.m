function [Mn,Mw,PDI,weight,T_unit,DB,dist_to_core] = result()
global chain;
global conversion_record;
%% 
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
%% 
[Mn,Mw,PDI,weight,T_unit,DB,dist_to_core] = calculate();
%% 
figure(3);
weight = weight(weight ~= 0);
hist(weight,10);
title('weight distribution');
% figure(4);
% global MW_record;
% plot(conversion_record,MW_record);%plot Molecular weight(weight average) versus conversion
figure(5);
global DB_record;
plot(conversion_record,DB_record);
title('支化度随转化率变化');
%% 


