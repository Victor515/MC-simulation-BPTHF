function [Mn,Mw,PDI,weight] = calculate()
%calculate molecular weight
global POLYMER;
[Mn,Mw,PDI,weight] = mw_calculate(POLYMER);
end

function [Mn,Mw,PDI,weight] = mw_calculate(polymer)
global chain;
polymer_temp = polymer;
weight = zeros(1,2000);
%delete all empty cells in polymer_temp, the length of remaining cell
%arrays is the number of polymer chains
polymer_temp(cellfun(@isempty,polymer_temp)) = [];

%calculate molecular weight of each chain, assume chains are terminated by
%water
for i = 1:length(polymer_temp)
    for j = 1:length(polymer_temp{i})
        weight(i) = weight(i) + 160 + chain(j).inserted_THF * 72;
    end
end

%calculate Mn, Mw and PDI
Mn = sum(weight) / length(weight);
Mw = weight * weight'/ sum(weight);
PDI = Mw / Mn;
end

