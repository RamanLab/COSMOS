%to find cross-fed metabolites in a co-culture
%load models of organisms
Soneidensis=readCbModel('iSO783_shewan.mat');
Kpneumoniae = readCbModel('iYL1228_klebsi.mat');

%define list of species names and biomass reactions
speciesList = {Soneidensis, Kpneumoniae};
biomassList = {'Biomass','Biomass'};
spNameList = {'Soneidensis','Kpneumoniae'};

%create community models
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);

%define size of timestep and total time
options.delt = 0.1;
options.maxTime = 12;
options.solver = 'ibm_cplex';

%define initial biomass concentration
options.initBiomass = [0.1,0.1];
options.abdCutoff = 10;

%% aerobic rich medium
%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults1 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    OptionsAll{1,i} = options;
    CommExchangeMets{i} = options.Products;
    allResults1{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{1},CommunityModelInfo{1}] = compileCom(Communities,pairedModelInfo,options,allResults1);

%% aerobic minimal medium
%define minimal medium
options.mediumMets = {'glc__D_e[u]';'4abz_e[u]'; 'btn_e[u]'; 'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'pnto__R_e[u]';'ribflv_e[u]';'ura_e[u]';'h2s_e[u]';'so4_e[u]';'zn2_e[u]';'sel_e[u]';'xan_e[u]';'thm_e[u]';'glu__L_e[u]';'leu__L_e[u]';'thr__L_e[u]';'val__L_e[u]';'ile__L_e[u]';'arg__L_e[u]';'ser__L_e[u]';'nac_e[u]';'so3_e[u]';'ala__L_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults2 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    allResults2{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{2},CommunityModelInfo{2}] = compileCom(Communities,pairedModelInfo,options,allResults2);

%% anaerobic rich medium
%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%setup anaerobic conditions
options.mediumMets = setdiff(options.mediumMets, 'o2_e[u]','stable'); %remove oxygen

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults3 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    allResults3{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{3},CommunityModelInfo{3}] = compileCom(Communities,pairedModelInfo,options,allResults3);


%% anaerobic minimal medium
%define minimal medium
options.mediumMets = {'glc__D_e[u]';'4abz_e[u]'; 'btn_e[u]'; 'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'pnto__R_e[u]';'ribflv_e[u]';'ura_e[u]';'h2s_e[u]';'so4_e[u]';'zn2_e[u]';'sel_e[u]';'xan_e[u]';'thm_e[u]';'glu__L_e[u]';'leu__L_e[u]';'thr__L_e[u]';'val__L_e[u]';'ile__L_e[u]';'arg__L_e[u]';'ser__L_e[u]';'nac_e[u]';'so3_e[u]';'ala__L_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults4 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    allResults4{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{4},CommunityModelInfo{4}] = compileCom(Communities,pairedModelInfo,options,allResults4);


%compute cross-fed metabolites in all four environments
CrossFedMets_Comm1 = []; CrossFedMets_Comm2 =[];
CrossFedMets_Comm3 = []; CrossFedMets_Comm4=[];
threshold = 0.001; %threshold for metabolite to be considered cross-fed
eumoniae = readCbModel('iYL1228_klebsi.mat');

%define list of species names and biomass reactions
speciesList = {Soneidensis, Kpneumoniae};
biomassList = {'Biomass','Biomass'};
spNameList = {'Soneidensis','Kpneumoniae'};

%create community models
[Communities,pairedModelInfo] = createPairwiseCommunity(speciesList,biomassList,spNameList);

%define size of timestep and total time
options.delt = 0.1;
options.maxTime = 12;
options.solver = 'ibm_cplex';

%define initial biomass concentration
options.initBiomass = [0.1,0.1];
options.abdCutoff = 10;

%% aerobic rich medium
%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults1 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    OptionsAll{1,i} = options;
    CommExchangeMets{i} = options.Products;
    allResults1{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{1},CommunityModelInfo{1}] = compileCom(Communities,pairedModelInfo,options,allResults1);

%% aerobic minimal medium
%define minimal medium
options.mediumMets = {'glc__D_e[u]';'4abz_e[u]'; 'btn_e[u]'; 'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'pnto__R_e[u]';'ribflv_e[u]';'ura_e[u]';'h2s_e[u]';'so4_e[u]';'zn2_e[u]';'sel_e[u]';'xan_e[u]';'thm_e[u]';'glu__L_e[u]';'leu__L_e[u]';'thr__L_e[u]';'val__L_e[u]';'ile__L_e[u]';'arg__L_e[u]';'ser__L_e[u]';'nac_e[u]';'so3_e[u]';'ala__L_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults2 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    allResults2{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{2},CommunityModelInfo{2}] = compileCom(Communities,pairedModelInfo,options,allResults2);

%% anaerobic rich medium
%define rich medium
options.mediumMets = {'glc__D_e[u]';'ala__L_e[u]';'arg__L_e[u]';'asn__L_e[u]';'asp__L_e[u]';'cys__L_e[u]';'gln__L_e[u]';'glu__L_e[u]';'gly_e[u]';'his__L_e[u]';'ile__L_e[u]';'leu__L_e[u]';'lys__L_e[u]';'met__L_e[u]';'phe__L_e[u]';'pro__L_e[u]';'ser__L_e[u]';'thr__L_e[u]';'trp__L_e[u]';'tyr__L_e[u]';'val__L_e[u]';'4abz_e[u]';'adocbl_e[u]';'btn_e[u]';'ca2_e[u]';'cbl1_e[u]';'chol_e[u]';'cl_e[u]';'co2_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'gua_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'nac_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'orot_e[u]';'pb2_e[u]';'h2s_e[u]';'xan_e[u]';'pheme_e[u]';'pi_e[u]';'pime_e[u]';'pnto__R_e[u]';'pydx_e[u]';'ribflv_e[u]';'sel_e[u]';'so3_e[u]';'so4_e[u]';'thm_e[u]';'thymd_e[u]';'ura_e[u]';'zn2_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%setup anaerobic conditions
options.mediumMets = setdiff(options.mediumMets, 'o2_e[u]','stable'); %remove oxygen

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults3 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    allResults3{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{3},CommunityModelInfo{3}] = compileCom(Communities,pairedModelInfo,options,allResults3);


%% anaerobic minimal medium
%define minimal medium
options.mediumMets = {'glc__D_e[u]';'4abz_e[u]'; 'btn_e[u]'; 'ca2_e[u]';'cbl1_e[u]';'cl_e[u]';'cobalt2_e[u]';'cu2_e[u]';'fe2_e[u]';'fe3_e[u]';'fol_e[u]';'h2_e[u]';'h2o_e[u]';'k_e[u]';'mg2_e[u]';'mn2_e[u]';'mobd_e[u]';'na1_e[u]';'ncam_e[u]';'nh4_e[u]';'ni2_e[u]';'no3_e[u]';'o2_e[u]';'pi_e[u]';'pnto__R_e[u]';'ribflv_e[u]';'ura_e[u]';'h2s_e[u]';'so4_e[u]';'zn2_e[u]';'sel_e[u]';'xan_e[u]';'thm_e[u]';'glu__L_e[u]';'leu__L_e[u]';'thr__L_e[u]';'val__L_e[u]';'ile__L_e[u]';'arg__L_e[u]';'ser__L_e[u]';'nac_e[u]';'so3_e[u]';'ala__L_e[u]'};
options.carbonSource = 'glc__D_e[u]';

%define medium concentration
options.initMedium(1:length(options.mediumMets),1) = 10;
% define Vmax, Km
[options.Vmax,options.Km] = deal(zeros(length(options.mediumMets),1));
options.Vmax(:) = 20;
options.Km(:) = 0.05;

%set medium for the communities
for i = 1: length(Communities)
    Communities{1,i} = setMediumCom(Communities{1,i}, options.mediumMets,options.initMedium);
end

allResults4 = cell(1, size(Communities, 2));

for i= 1 : size(Communities,2)
    options.Products={}; options.ProductName = {};
    for k=1:length(Communities{1,i}.infoCom.EXsp)
        if (Communities{1,i}.indCom.EXsp(k,1)>0) && (Communities{1,i}.indCom.EXsp(k,2)>0)
            options.ProductName(end+1) = Communities{1,i}.infoCom.Mcom(k);
        end
    end
    for j=1:length(options.ProductName)
        options.Products{j} = erase(options.ProductName{j},'[u]');
        options.Products{j} = strcat('EX_',options.Products{j});
    end
    allResults4{i} = dFBAComCross(Communities{1,i},options); %uses fastFVA with ibmcplex and fluxVariability for other solvers
end

%compile community results using abundance cutoff
[CommunityResults{4},CommunityModelInfo{4}] = compileCom(Communities,pairedModelInfo,options,allResults4);


%compute cross-fed metabolites in all four environments
CrossMets = []; 
threshold = 0.001; %threshold for metabolite to be considered cross-fed

for i=1:length(Communities)
    for j = 1:4
        for k= 1:length(CommExchangeMets{1,i})
            Exc1 = []; Exc2 = [];
            Exc1 = CommunityResults{1,j}{1,i}.prdtFVAConc1(k,:);
            Exc2 = CommunityResults{1,j}{1,i}.prdtFVAConc2(k,:);
            for m=1:length(Exc1)
                if ((Exc1(m)>=threshold) && (Exc2(m)<=-(threshold)) || ((Exc1(m)<=-(threshold)) && (Exc2(m)>=threshold)))
                    CrossMets{1,i}{j,1}{k,m}=1;
                else
                    CrossMets{1,i}{j,1}{k,m}=0;
                end
            end           
        end
        CrossFedMetsComm{1,i}{j,1} = [(CommExchangeMets{1,i})' (CrossMets{1,i}{j,1})];
    end
end
