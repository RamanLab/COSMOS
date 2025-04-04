% To analyze the sensitivity of community growth to Vmax and Km parameters
%% LHS sampling
% Define parameter bounds
Vmax_range = [1, 100];   % mmol/gDW/hr
Km_range = [0.01, 1];    % mM

% Number of samples 
N = 10;  

% Generate Latin Hypercube Samples
lhs_samples = lhsdesign(N, 2);  

% Scale samples to parameter ranges
Vmax_samples = Vmax_range(1) + lhs_samples(:,1) * (Vmax_range(2) - Vmax_range(1));
Km_samples = Km_range(1) + lhs_samples(:,2) * (Km_range(2) - Km_range(1));

% Combine into parameter sets
parameter_sets = [Vmax_samples, Km_samples];

%% run aerobic rich analysis
%load models of aerobic organisms
Ecoli = readCbModel('iML1515.mat'); Scerevisiae = readCbModel('iMM904_Scerevisiae.mat');
Bsubtilis = readCbModel('iYO844_subtilis.mat'); Pputida = readCbModel('iJN1463_putida.mat');
Llactis = readCbModel('iNF517_lactis.mat'); Scystis = readCbModel('iJN678_sycystis.mat');
Paeruginosa = readCbModel('iMO1056_pseudo.mat');Soneidensis=readCbModel('iSO783_shewan.mat');
Kpneumoniae = readCbModel('iYL1228_klebsi.mat');

%define list of species names and biomass reactions
speciesList = {Ecoli,Scerevisiae,Pputida,Llactis,Bsubtilis,Scystis,Paeruginosa, Soneidensis, Kpneumoniae};
biomassList = {'BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_SC5_notrace','BIOMASS_KT2440_WT3','BIOMASS_LLA','BIOMASS_BS_10','BIOMASS_Ec_SynHetero','PA_Biomass6_DM','Biomass','Biomass'};
spNameList = {'Ecoli','Scerevisiae','Pputida','Llactis','Bsubtilis','Scystis','Paeruginosa','Soneidensis','Kpneumoniae'};

%create community models
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);

%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;

%define products and other options
options.Products = {'EX_succ_e','EX_pyr_e','EX_for_e','EX_lac__L_e','EX_lac__D_e','EX_ac_e','EX_fum_e','EX_glcn_e','EX_ppa_e','EX_adpac_e','EX_sbt__D_e','EX_xylt_e','EX_etoh_e','EX_meoh_e','EX_btoh_e','EX_12ppd__R_e','EX_13ppd_e','EX_btd_RR_e','EX_glyc_e','EX_h2_e','EX_but_e','EX_spmd_e','EX_ptrc_e','EX_catechol_e','EX_gthrd_e'}';
options.ProductName = {'Succinate','Pyruvate','Formate','L-Lactate','D-Lactate','Acetate','Fumarate','Gluconate','Propionate','Adipic acid','Sorbitol','Xylitol','Ethanol','Methanol','Butanol','Propane-1,2-diol','Propane-1,3-diol','2,3 Butanediol','Glycerol','Hydrogen','Butyrate','Spermidine','Putrescine','Catechol','Glutathione'}';
options.orgCount = length(speciesList);

%define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
%define size of timestep and total time
options.delt = 0.1;
options.maxTime = 12;
options.solver = 'ibm_cplex';

%define initial biomass concentration
options.initBiomass = [0.1,0.1];

%define the minimum %improvement in biomass and abundance required for
%communityto be considered viable
options.abdCutoff = 10;
%set medium for the communities
    for i = 1: length(Communities)
        Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
    end

%% sensitivity analysis
% to evaluate parameter sensitivity
for p = 1: length(parameter_sets)
    options.Vmax(:) = parameter_sets(p,1);
    options.Km(:) = parameter_sets(p,2);
    
    %run dfba for all the pairwise communities
    environment = getEnvironment();
    allResults = cell(1, size(Communities, 2));

    for i= 1 : size(Communities,2)
        restoreEnvironment(environment);
        solverOK = changeCobraSolver(options.solver,'all'); 
        allResults{i} = dFBACom(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
    end

    %compile community results using abundance cutoff
    [CommunityResults{p},CommunityModelInfo{p}] = compileCom(Communities,pairedModelInfo,options,allResults);
end

CommGrowthA = zeros(length(parameter_sets),length(pairedModelInfo));
CommGrowthB = zeros(length(parameter_sets),length(pairedModelInfo));
CommGrowthTotal = zeros(length(parameter_sets),length(pairedModelInfo));


for p=1:length(parameter_sets)
    for i=1:length(CommunityModelInfo{1,p})
        [LocA,LocB] = ismember(CommunityModelInfo{1,p}{i,1},pairedModelInfo(:,1));
        if LocB>0
            CommGrowthA(p,LocB) = CommunityResults{1,p}{1,i}.biomass(1);
            CommGrowthB(p,LocB) = CommunityResults{1,p}{1,i}.biomass(2);
            CommGrowthTotal(p,LocB) = CommunityResults{1,p}{1,i}.biomass(1)+CommunityResults{1,p}{1,i}.biomass(2);
        end
    end
end
    
Vmax_values = repmat(Vmax_samples,1,length(pairedModelInfo));
Km_values = repmat(Km_samples,1,length(pairedModelInfo));

%flatten data
Vmax_flat = Vmax_values(:);
Km_flat = Km_values(:);
CommGrowthTotal_flat = CommGrowthTotal(:);

% Calculate correlation coefficients
corr_Vmax = corr(Vmax_flat, CommGrowthTotal_flat, 'Type', 'Pearson');
corr_Km = corr(Km_flat, CommGrowthTotal_flat, 'Type', 'Pearson');

% Construct the design matrix
X = [Vmax_flat, Km_flat];
X = [ones(size(X, 1), 1), X]; % Add intercept term
y = CommGrowthTotal_flat;

% Perform regression
[b,~,~,~,stats] = regress(y, X);

% Display regression results
fprintf('Regression coefficients: \n');
fprintf('Intercept: %.2f\n', b(1));
fprintf('Vmax: %.2f\n', b(2));
fprintf('Km: %.2f\n', b(3));
fprintf('R-squared: %.2f\n', stats(1));
