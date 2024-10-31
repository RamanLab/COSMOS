%load models of aerobic organisms
Soneidensis=readCbModel('iSO783_shewan.mat');
Kpneumoniae = readCbModel('iYL1228_klebsi.mat');

%define list of species names and biomass reactions
speciesList = {Soneidensis,Kpneumoniae};
biomassList = {'Biomass','Biomass'};
spNameList = {'Soneidensis','Kpneumoniae'};

%create community models
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);

%define medium
options.mediumMets = {'glyc_e[u]','k_e[u]','pi_e[u]','nh4_e[u]','mg2_e[u]','so4_e[u]','ca2_e[u]','fe2_e[u]','cl_e[u]','zn2_e[u]','man_e[u]','cobalt2_e[u]','cu2_e[u]','ni2_e[u]','na1_e[u]','mobd_e[u]','hco3_e[u]'}';
options.initMedium = [543.02,23.31,23.31,15.15,7.89,15.488,40.018,0.018,0.0448,0.001,0.001,0.002,0.0002,0.0002,0.002,0.0003,9.55]';

%define products and other options
options.Products = {'EX_13ppd_e'}';
options.ProductName = {'Propane-1,3-diol'}';
options.orgCount = length(speciesList);

%define Vmax, Km and delt
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;
options.delt = 0.1;
options.carbonSource = 'glyc_e[u]';
options.maxTime = 100;
biomassRatios = [0.066 0.133; 0.1 0.1; 0.133 0.066; 0.15 0.05];

% Store results for Pareto front computation
CommunityBMResults = struct('biomassRatio', {}, 'totalBiomass', {}, 'totalProduct', {});

allResults={};
% Loop through all biomass ratios
for r = 1:size(biomassRatios, 1)
    options.initBiomass = biomassRatios(r, :);
    
    % Set medium for each community
    for i = 1:length(Communities)
        Community = setMediumCom(Communities{1, i}, options.mediumMets, options.initMedium);
    end
    
    % Run dFBA simulations for all communities
    allResults{1,r} = sokp_fedBatch(Community, options);
    [CommunityResults{1,r}] = compileSingleCom(Community,options,allResults{1,r});
    
    % Store the results for this biomass ratio
    CommunityBMResults(r).biomassRatio = options.initBiomass;
    CommunityBMResults(r).totalBiomass = sum(CommunityResults{1,r}.biomass);
    CommunityBMResults(r).totalProduct = CommunityResults{1,r}.totPrdtConc;
    
end

%monoculture simulation-sokp
optionsMono = options;
glycerolConc = [434.42,543.02,760.23];
optionsMono.initBiomass = [0.1,0.1]; %resetting biomass ratio

for i=1:length(glycerolConc)
    %change glycerol concentration
    optionsMono.initMedium(find(ismember(options.mediumMets,'glyc_e[u]')),1) = glycerolConc(i);
    
    % Set medium for each community
    Communities{1,1} = setMediumCom(Communities{1, 1}, optionsMono.mediumMets, optionsMono.initMedium);
    [MonoResults{i},~] = monocultureGrowth(Communities,pairedModelInfo,optionsMono);
    
end

