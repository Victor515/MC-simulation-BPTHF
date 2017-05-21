function [Mn, Mw, PDI,avg_T, avg_DB, DB] = main(rate_ratio, feed_ratio,conversion)
% set reaction rate constant ratio, kTHF: kPO
global RATE_RATIO
RATE_RATIO = rate_ratio;

%set reactant numbers, EGDE and THF
global FEED_RATIO %EGDE:THF
FEED_RATIO = feed_ratio;
global EGDE_NUM
EGDE_NUM = 2000;
global THF_NUM
THF_NUM = EGDE_NUM / FEED_RATIO;

% represent EGDE chains using structure arrays
% inserted_THF -- how many THF are in chain i
% inserted_chain_pos -- the position where other chains inserts,array is
% used because muptile chains could be inserted
% inserted_chain_serial -- the serial number of above chains
% chain_inserting -- the serial number of the chain which chain i
% inserts into
% polymer_num -- records the number of the polymer the chain belongs to, one chain
% could only belong to one polymer
global chain;
for i = 1:EGDE_NUM
    chain(i).inserted_THF = 0;
    chain(i).inserted_chain_pos = [];
    chain(i).inserted_chain_serial = [];
    chain(i).chain_inserting = 0;
    chain(i).polymer_num = 0;
end


%conduct reactions, record parameters
rng('default');
rng('shuffle'); %generate different random sequence
reaction(conversion);

%demonstrate final results
[Mn,Mw,PDI,weight,T_unit,DB,dist_to_core] = result();
avg_T = mean(T_unit); % average number of terminal units per chain
avg_DB =  mean(DB); % average degree of branching
% mean(dist_to_core) % average distance to the core


