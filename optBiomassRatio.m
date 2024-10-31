%% run optimal biomass analysis
%load models of aerobic organisms
Soneidensis=readCbModel('iSO783_shewan.mat');
Kpneumoniae = readCbModel('iYL1228_klebsi.mat');

speciesList = {Soneidensis,Kpneumoniae};
biomassList = {'Biomass','Biomass'};
spNameList = {'Soneidensis','Kpneumoniae'};

%create community models
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);
Community = Communities{1,1};

%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};

%setup anaerobic conditions
options.mediumMets = setdiff(options.mediumMets, 'o2_e[u]','stable'); %remove oxygen

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;

%define products and other options
options.Products = {'EX_13ppd_e'}';
options.ProductName = {'1,3-propanediol'}';
options.orgCount = length(speciesList);
options.carbonSource = 'glc__D_e[u]';

%define Vmax, Km and delt
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;
%define size of timestep and total time
options.delt = 0.1;
options.maxTime = 12;

%define the minimum %improvement in biomass and abundance required for
%communityto be considered viable
options.abdCutoff = 10;

% Define a range of initial biomass ratios to explore
biomassRatios = [0.1 0.9; 0.2 0.8; 0.3 0.7;0.4 0.6; 0.5 0.5; 0.6 0.4; 0.7 0.3; 0.8 0.2; 0.9 0.1];
biomassRatios = biomassRatios/5; %so that sum of initial biomass is 0.2 as other analysis
options.maxTime = 12;
% Store results for all biomass ratios
CommunityBMResults = struct('biomassRatio', {}, 'totalBiomass', {}, 'totalProduct', {});

tic
allResults={};
% Loop through all biomass ratios
for r = 1:size(biomassRatios, 1)
    options.initBiomass = biomassRatios(r, :);
    
    % Set medium for each community
    Community = setMediumCom(Community, options.mediumMets, options.initMedium);
    
    % Run dFBA simulations for all communities
    allResults{1,r} = dFBACom(Community, options);
    CommunityResults{1,r} = compileSingleCom(Community,options,allResults{1,r});
    
    % Store the results for this biomass ratio
    CommunityBMResults(r).biomassRatio = options.initBiomass;
    
    CommunityBMResults(r).biomassA = CommunityResults{1,r}.biomass(1);
    CommunityBMResults(r).biomassB = CommunityResults{1,r}.biomass(2);
    
    CommunityBMResults(r).totalBiomass = sum(CommunityResults{1,r}.biomass);
    CommunityBMResults(r).totalProduct = CommunityResults{1,r}.totPrdtConc;
    
    CommunityBMResults(r).abundance1 = CommunityResults{1,r}.abundance(1);
    CommunityBMResults(r).abundance2 = CommunityResults{1,r}.abundance(2);
    
end
toc
