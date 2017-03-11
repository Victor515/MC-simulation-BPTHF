function reaction()
global chain;
global THF_NUM;
global EGDE_NUM;
initial_THF = THF_NUM;
initial_EGDE = length(chain);
conversion = (initial_THF - THF_NUM) / initial_THF;

%set up arrays to record instant values
global conversion_record;%record instant conversion rate
conversion_record = [];
global MW_record;%record instant Molecular weight
MW_record = [];
global T_unit_record;
T_unit_record = [];
global DB_record;
DB_record = [];
global dist_to_core_record;
dist_to_core_record = [];

% set up a cell array POLYMER,to record branched structures and calculate
% molecular weight
global POLYMER
POLYMER = cell(1,2000);

% % set up an array to record all cyclized polymer
% global CYCL_POLYMER
% CYCL_POLYMER = [];

polynum = 1; %a variable to count POLYMER chains generated
while conversion < 0.4
%choose the chain to be propagated
chain_serial = choose_chain(length(chain)); 

%decide the reaction type(THF or EGDE insertion) for the chosen chain, modify relevant parameters based on
%the reaction
reaction_type = choose_reaction(length(chain));
switch reaction_type
    case 0 % A THF is inserted
        chain(chain_serial).inserted_THF = chain(chain_serial).inserted_THF + 1;
        THF_NUM = THF_NUM -1;
    case 1 % A EGDE-PO is inserted
        %check whether ring clevage of EGDE has happened, if not, insert
        %one THF first
        if chain(chain_serial).inserted_THF == 0
            chain(chain_serial).inserted_THF = 1;
            THF_NUM = THF_NUM - 1;
        else
            %select the inserting chain
            chain_serial_ins = select_chain(length(chain));
            
%             if (chain_serial == chain_serial_ins && chain(chain_serial).polymer_num ~= 0)
%                 CYCL_POLYMER = [CYCL_POLYMER, chain(chain_serial).polymer_num];
%             end
            %modify POLYMER
            if (chain(chain_serial).polymer_num == 0 && chain(chain_serial_ins).polymer_num == 0) %if neither of inserted and inserting chain belongs to any polymer
                if chain_serial ~= chain_serial_ins
                    POLYMER{polynum} = [chain_serial,chain_serial_ins];
                    chain(chain_serial).polymer_num = polynum;
                    chain(chain_serial_ins).polymer_num = polynum;
                else
                    POLYMER{polynum} = chain_serial;
                    chain(chain_serial).polymer_num = polynum;
                    chain(chain_serial_ins).polymer_num = polynum;
%                     CYCL_POLYMER = [CYCL_POLYMER, chain(chain_serial).polymer_num];
                end
                polynum = polynum + 1;
            elseif (chain(chain_serial).polymer_num == 0 && chain(chain_serial_ins).polymer_num ~= 0) %if inserting chain belongs to any polymer, while inserted chain does not
                    POLYMER{chain(chain_serial_ins).polymer_num} = [chain_serial, POLYMER{chain(chain_serial_ins).polymer_num}]; % add inserted chain to the beginning
                    chain(chain_serial).polymer_num = chain(chain_serial_ins).polymer_num;
            elseif (chain(chain_serial_ins).polymer_num == 0 && chain(chain_serial).polymer_num ~= 0) %if inserted chain belongs to any polymer, while inserting chain does not
                    POLYMER{chain(chain_serial).polymer_num} = [POLYMER{chain(chain_serial).polymer_num},chain_serial_ins];
                    chain(chain_serial_ins).polymer_num = chain(chain_serial).polymer_num;
            else    %if both chains belong to any polymer, union all records and store in polymer of inserted chain, delete polymer of inserting chain
                    if chain(chain_serial).polymer_num ~= chain(chain_serial_ins).polymer_num %make sure two chains do not belong to the same polymer
                        POLYMER{chain(chain_serial).polymer_num} = [POLYMER{chain(chain_serial).polymer_num}, POLYMER{chain(chain_serial_ins).polymer_num}];
                        temp_poly_num = chain(chain_serial_ins).polymer_num;
                        for i = POLYMER{temp_poly_num}
                            chain(i).polymer_num = chain(chain_serial).polymer_num;
                        end
                        POLYMER{temp_poly_num} = [];
%                     else % if a chain inserts into itself, a cycle will form in the polymer the chain belong to
%                         CYCL_POLYMER = [CYCL_POLYMER, chain(chain_serial).polymer_num];
                    end
            end
            
            
            %modify inserted_chain_pos and inserted_chain_serial, inserted_THF is the position the chain
            %inserts into
            chain(chain_serial).inserted_chain_pos = [chain(chain_serial).inserted_chain_pos, chain(chain_serial).inserted_THF];
            chain(chain_serial).inserted_chain_serial = [chain(chain_serial).inserted_chain_serial, chain_serial_ins];
            
            %modify inserting_chain
            chain(chain_serial_ins).chain_inserting = chain_serial;
            
            %Decrease EGDE_NUM by 1
            EGDE_NUM = EGDE_NUM -1;
            
        end
        
end

%calculate characterization parameters
conversion = (initial_THF - THF_NUM) / initial_THF;
conversion
% % record instant conversion and parameters
% if ismember(conversion,[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
% conversion_record = [conversion_record,conversion]; 
% [Mn,Mw,PDI,weight,T_unit,DB,dist_to_core] = calculate();
% T_unit_record = [T_unit_record,mean(T_unit)];
% DB_record = [DB_record, mean(DB)];
% dist_to_core_record = [dist_to_core_record, mean(dist_to_core)];
% end
end

end


function chain_serial = choose_chain(maxium_num)
chain_serial = unidrnd(maxium_num);
end


function chain_serial = select_chain(maxium_num)
%When a PO inserts into another chain, we need to determine which chain it
%belongs to, based on random num. This function is similar to
%choose_chain()
global chain;
while 1
    chain_serial = unidrnd(maxium_num);
    if chain(chain_serial).chain_inserting == 0 %determine whether this chain already inserts into other chains
        break;
    end
end
end

function reaction_type = choose_reaction(chain_num)
global EGDE_NUM;
global THF_NUM;
global RATE_RATIO;
%generate a random number
r = rand;

%determine which reaction to proceed
if r > EGDE_NUM * chain_num / (EGDE_NUM*chain_num + RATE_RATIO * THF_NUM * chain_num)
    reaction_type = 0; %A THF is inserted
else
    reaction_type = 1; %A EGDE-PO is inserted
end
end



        

