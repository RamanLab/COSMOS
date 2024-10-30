%% perform carbon source analysis
%aerobic rich environment

%% initialize models and parameters
%load models of aerobic organisms
Ecoli = readCbModel('iML1515.mat'); Scerevisiae = readCbModel('iMM904_Scerevisiae.mat');
Bsubtilis = readCbModel('iYO844_subtilis.mat'); Pputida = readCbModel('iJN1463_putida.mat');
Llactis = readCbModel('iNF517_lactis.mat'); Scystis = readCbModel('iJN678_sycystis.mat');
Paeruginosa = readCbModel('iMO1056_pseudo.mat');Soneidensis=readCbModel('iSO783_shewan.mat');
Kpneumoniae = readCbModel('iYL1228_klebsi.mat');Tmaritima = readCbModel('iLJ478_thermo.mat');
Mbarkeri = readCbModel('iAF692_barkerii.mat');

%define list of species names and biomass reactions
speciesList = {Ecoli,Scerevisiae,Pputida,Llactis,Bsubtilis,Scystis,Paeruginosa, Soneidensis, Kpneumoniae};
biomassList = {'BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_SC5_notrace','BIOMASS_KT2440_WT3','BIOMASS_LLA','BIOMASS_BS_10','BIOMASS_Ec_SynHetero','PA_Biomass6_DM','Biomass','Biomass'};
spNameList = {'Ecoli','Scerevisiae','Pputida','Llactis','Bsubtilis','Scystis','Paeruginosa','Soneidensis','Kpneumoniae'};

%create community models
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);

%define rich medium
options.mediumMets_ini = {'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};

%define medium concentration
options.initMedium_ini(1:length(options.mediumMets_ini),1) = 10;

%define products and other options
options.Products = {'EX_succ_e','EX_pyr_e','EX_for_e','EX_lac__L_e','EX_lac__D_e','EX_ac_e','EX_fum_e','EX_glcn_e','EX_ppa_e','EX_adpac_e','EX_sbt__D_e','EX_xylt_e','EX_etoh_e','EX_meoh_e','EX_btoh_e','EX_12ppd__R_e','EX_13ppd_e','EX_btd_RR_e','EX_glyc_e','EX_h2_e','EX_but_e','EX_spmd_e','EX_ptrc_e','EX_xan_e','EX_gthrd_e'}';
options.ProductName = {'Succinate','Pyruvate','Formate','L-Lactate','D-Lactate','Acetate','Fumarate','Gluconate','Propionate','Adipic acid','Sorbitol','Xylitol','Ethanol','Methanol','Butanol','Propane-1,2-diol','Propane-1,3-diol','2,3 Butanediol','Glycerol','Hydrogen','Butyrate','Spermidine','Putrescine','Xanthine','Glutathione'}';
options.orgCount = length(speciesList);

%define initial biomass concentration
options.initBiomass = [0.1,0.1];

%define size of timestep and total time
options.delt = 0.1;
options.maxTime = 12;

%define the minimum %improvement in biomass and abundance required for
%communityto be considered viable
options.abdCutoff = 10;

%define carbon source list
carbonSourceList = {'glc__D_e[u]';'fru_e[u]';'sucr_e[u]';'malt_e[u]';'lac__L_e[u]';'xyl__D_e[u]';'gal_e[u]';};
csMedium(1:length(carbonSourceList),1) = 10;

%% cs analysis
environment = getEnvironment();
parfor cs = 1:length(carbonSourceList)
    restoreEnvironment(environment);
    solverOK = changeCobraSolver('ibm_cplex', 'all');
    [Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);
    
    optionsCS = options;
    optionsCS.mediumMets = optionsCS.mediumMets_ini;
    optionsCS.initMedium = optionsCS.initMedium_ini;
    optionsCS.carbonSource = carbonSourceList(cs);
    
    %to accomodate presence of carbon source in initial media
    if ~ismember(optionsCS.mediumMets, optionsCS.carbonSource)
        optionsCS.mediumMets = vertcat(optionsCS.mediumMets, optionsCS.carbonSource);
        optionsCS.initMedium(end+1) = csMedium(cs);
        optionsCS.csMedium = csMedium(cs);
    else
        optionsCS.initMedium(find(ismember(optionsCS.mediumMets, optionsCS.carbonSource)), 1) = ...
        optionsCS.initMedium(find(ismember(optionsCS.mediumMets, optionsCS.carbonSource)), 1) + csMedium(cs);
        optionsCS.csMedium = optionsCS.initMedium(find(ismember(optionsCS.mediumMets, optionsCS.carbonSource)), 1);
    end
    
    % Define Vmax, Km
    [optionsCS.Vmax, optionsCS.Km] = deal(zeros(length(optionsCS.mediumMets), 1));
    optionsCS.Vmax(:) = 20;
    optionsCS.Km(:) = 0.05;
    
    % Set medium for the communities
    for i = 1:length(Communities)
        Communities{1, i} = setMediumCom(Communities{1, i}, optionsCS.mediumMets, optionsCS.initMedium);
    end
    
    % Run dFBA for all the pairwise communities
    allResults = cell(1, size(Communities, 2));
    for i = 1:size(Communities, 2)
        allResults{i} = dFBACom(Communities{1, i}, optionsCS);
    end
    
    % Initialize result containers
    CommunityResults = {};
    CommunityModelInfo = {};
    MonoResults = {};
    singleModelInfo = {};
    pairwiseInteractions = {};
    PrdtAnalysis = {};
    interactionGeneric = {};
    interactionSpecific = {};
    
    %compile results for communities using abundance cutoff
    [CommunityResults, CommunityModelInfo] = compileCom(Communities, pairedModelInfo, optionsCS, allResults);
    
    % Perform monoculture analysis
    [MonoResults, singleModelInfo] = monocultureGrowth(Communities, pairedModelInfo, optionsCS);
    
    % Perform interaction analysis if results exist
    if ~isempty(CommunityResults) && ~isempty(MonoResults)
        %computes interaction-type for all communities
        pairwiseInteractions = pairInteractions(CommunityResults,CommunityModelInfo,MonoResults,singleModelInfo);
        %comparison of productivity - alternatively use YieldAnalysis function
        [PrdtAnalysis,BestSystem] = ProductivityAnalysis(CommunityResults,CommunityModelInfo,MonoResults,singleModelInfo,options);
        %find influence of interaction on productivity
        [interactionGeneric,interactionSpecific] = InteractionAnalysis(PrdtAnalysis,pairwiseInteractions);
        
        % Store results
        CSResults{cs}.pairwiseInteractions = pairwiseInteractions;
        CSResults{cs}.PrdtAnalysis = PrdtAnalysis;
        CSResults{cs}.BestSystem = BestSystem;
        CSResults{cs}.interactionGeneric = interactionGeneric;
        CSResults{cs}.interactionSpecific= interactionSpecific;
    else
        disp('No viable communities found, so comparison analysis was not performed');
    end
    CSResults{cs}.carbonSource = optionsCS.carbonSource;
    CSResults{cs}.options = optionsCS;
    CSResults{cs}.CommunityResults = CommunityResults;
    CSResults{cs}.CommunityModelInfo = CommunityModelInfo;
    CSResults{cs}.MonoResults = MonoResults;
    CSResults{cs}.singleModelInfo = singleModelInfo;
    CSResults{cs}.rawResults = allResults;
end
