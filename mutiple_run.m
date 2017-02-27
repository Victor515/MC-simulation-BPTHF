Mn = zeros(1,10);
Mw = zeros(1,10);
PDI = zeros(1,10);
for i = 1:10
    [Mn(i),Mw(i),PDI(i)] = main();
end
fprintf('Mn = %d, Mw = %d, PDI = %d',mean(Mn),mean(Mw),mean(PDI));