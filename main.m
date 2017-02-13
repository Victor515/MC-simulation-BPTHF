function [Mn,Mw,PDI] = main()
% set reaction rate constant ratio, kTHF: kPO
global RATE_RATIO
RATE_RATIO = 1/2;

%set reactant numbers, EGDE and THF
global EGDE_NUM
EGDE_NUM = 2000;
global THF_NUM
THF_NUM = 40*EGDE_NUM;

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
reaction();

%demonstrate final results
[Mn, Mw, PDI] = result();


