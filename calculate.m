function calculate()
%calculate molecular weight
global POLYMER;
[Mn,Mw,PDI,weight] = mw_calculate(POLYMER);
end

function [Mn,Mw,PDI,weight] = mw_calculate(polymer)
polymer_num = 1; %record the number of every single molecule
weight = zeros(2000); %pre-allocate memory for weight vector
for i = 1:2000
    if chain(i).chain_inserting == 0 %start from the chain that did not insert into any other chain
        weight(polymer_num) = weight_calculate(i);
        polymer_num = polymer_num + 1;
    end
end
end

function weight = weight_calculate(serial)
global chain;
weight = chain(serial).inserted_THF * THF_weight;
if ~isempty(chain(serial).inserted_chain_serial) %determine whether there is any branching point
    len = length(chain(serial).inserted_chain_serial);
    for j = 1:len
         weight = weight + weight_calculate(chain(serial).inserted_chain_serial(j));%recursion
    end
end
end