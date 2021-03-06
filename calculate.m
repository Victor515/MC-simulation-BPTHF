function [Mn,Mw,PDI,weight,T_unit,DB,dist_to_core,D_unit,L_unit] = calculate()
%calculate molecular weight
global POLYMER;
global chain;


[Mn,Mw,PDI,weight] = mw_calculate(POLYMER);

%calculate DB parameters
[T_unit,D_unit,L_unit,DB,dist_to_core] = db_calculate(chain,POLYMER);
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
        weight(i) = weight(i) + 160 + chain(polymer_temp{i}(j)).inserted_THF * 72;
    end
end

%calculate Mn, Mw and PDI
weight = weight(weight ~= 0);
Mn = sum(weight) / length(weight);
Mw = weight * weight'/ sum(weight);
PDI = Mw / Mn;
end

% calculate parameters related to branching structure
function [T_unit,D_unit,L_unit,DB,dist_to_core] = db_calculate(chain,polymer)
polymer_temp = polymer;

%delete all empty cells in polymer_temp, the length of remaining cell
%array is the number of polymer chains
polymer_temp(cellfun(@isempty,polymer_temp)) = [];

% %delete all cyclized polymer in polymer_temp
polymer_temp = cycle_rm(polymer_temp); %create a new data structure, only used for calculating dist_to_core

DB = zeros(1,length(polymer_temp));
T_unit = zeros(1,length(polymer_temp));
D_unit = zeros(1,length(polymer_temp));
L_unit = zeros(1,length(polymer_temp));
dist_to_core = zeros(1,length(polymer_temp));

% calculate number of terminal units, dendric units, linear units and
% degree of branching of every polymer
for i = 1:length(polymer_temp)
    T_unit(i) = length(polymer_temp{i});
    D_unit(i) = T_unit(i) - 1; % number of dendric units are always smaller than number of terminal units by one 
    for j = 1:length(polymer_temp{i})
        L_unit(i) = chain(polymer_temp{i}(j)).inserted_THF + L_unit(i);
    end
    DB(i) = (D_unit(i) + T_unit(i)) / (D_unit(i) + T_unit(i) + L_unit(i));
    
    
    % calculate max distance to core for every polymer chain
    for k = polymer_temp{i}
        index = k;
        n = 1;
        count(n) = 0; % count number of branching point traversed
        while chain(index).chain_inserting ~= 0
            count(n) = count(n) + 1;
            index = chain(index).chain_inserting;
        end
        n = n + 1;
    end
    dist_to_core(i) = max(count);
            
end
%     dist_to_core = dist_to_core(dist_to_core ~= 0);
end

% remove cycled polymers
function polymer = cycle_rm(polymer_temp)
cycle_polymer = []; %create a list to record all serial numbers of cycled polymers
for i = 1:length(polymer_temp)
    if iscycle(polymer_temp{i})
        cycle_polymer = [cycle_polymer, i];
    end
end
for k = cycle_polymer
    polymer_temp{k} = [];
end
polymer_temp(cellfun(@isempty,polymer_temp)) = [];
polymer = polymer_temp;
end

% determine whether a polymer forms a circle
function boolcycle = iscycle(polymer)
global chain;
boolcycle = false;
index = polymer(1);
while chain(index).chain_inserting ~= 0
    index = chain(index).chain_inserting;
    if index == polymer(1)
        boolcycle = true;
        break;
    end
end
end