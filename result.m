function [Mn,Mw,PDI,weight,T_unit,DB,dist_to_core] = result()
global chain;
global conversion_record;
global record_instant_value;
%%
% %% 
% for i = 1:length(chain)
%     a(i) = chain(i).inserted_THF;%THF number of every chain
%     b(i) = length(chain(i).inserted_chain_pos);%Number of Branching point of every chains
% end
% figure(1);
% hist(a);
% title('THF num distribution');
% figure(2);
% hist(b,0:10);
% title('Branching point num distribution');
%%
[Mn,Mw,PDI,weight,T_unit,DB,dist_to_core] = calculate();
%% 
% figure(3);
% weight = weight(weight ~= 0);
% hist(weight,10);
% title('weight distribution');
%%
if record_instant_value
    figure(4);
    global Mw_record;
    plot(conversion_record,Mw_record);%plot Molecular weight(weight average) versus conversion
    xlabel('转化率');
    ylabel('重均分子量');
    figure(5);
    global Mn_record;
    plot(conversion_record,Mn_record);%plot Molecular weight(number average) versus conversion
    xlabel('转化率');
    ylabel('数均分子量');
    figure(6);
    global PDI_record;
    plot(conversion_record,PDI_record);%plot Molecular weight(weight average) versus conversion
    xlabel('转化率');
    ylabel('PDI');
    figure(7);
    global DB_record;
    plot(conversion_record,DB_record);
    xlabel('转化率');
    ylabel('支化度');
    figure(8);
    global T_unit_record;
    plot(conversion_record,T_unit_record);
    xlabel('转化率');
    ylabel('末端数');
    figure(9);
    global dist_to_core_record;
    plot(conversion_record,dist_to_core_record);
    xlabel('转化率');
    ylabel('末端到核心距离');
end
%% 
figure(10);
gpcplot(weight,1);
%%


