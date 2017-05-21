test = [1/0.05,1/0.1,1/0.12,1/0.14,1/0.15,1/0.17,1/0.19,1/0.2];
count = 1;
for i = test
     [Mn(count),Mw(count),PDI(count),avg_T(count), avg_DB(count)] = main(i);
     count = count + 1;
end