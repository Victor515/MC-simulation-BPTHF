function [Mn,Mw,PDI,weight,T_unit,DB] = calculate()
%calculate molecular weight
global POLYMER;
global chain;
[Mn,Mw,PDI,weight] = mw_calculate(POLYMER);
[T_unit,D_unit,L_unit,DB] = db_calculate(chain,POLYMER);
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

function [T_unit,D_unit,L_unit,DB] = db_calculate(chain,polymer)
polymer_temp = polymer;

%delete all cyclized polymer in polymer_temp
global CYCL_POLYMER;
for i = CYCL_POLYMER
    polymer_temp(i) = [];
end

%delete all empty cells in polymer_temp, the length of remaining cell
%arrays is the number of polymer chains
polymer_temp(cellfun(@isempty,polymer_temp)) = [];

DB = zeros(1,length(polymer_temp));
T_unit = zeros(1,length(polymer_temp));
D_unit = zeros(1,length(polymer_temp));
L_unit = zeros(1,length(polymer_temp));

% calculate number of terminal units, dendric units, linear units and
% degree of branching of every polymer
for i = 1:length(polymer_temp)
    T_unit(i) = length(polymer_temp{i});
    D_unit(i) = T_unit(i) - 1; % number of dendric units are always smaller than number of terminal units by one 
    for j = 1:length(polymer_temp{i})
        L_unit(i) = chain(polymer_temp{i}(j)).inserted_THF + L_unit(i);
    end
    DB(i) = (T_unit(i) + D_unit(i)) / (T_unit(i) + D_unit(i) + L_unit(i));
end
    
end

