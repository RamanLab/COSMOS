%find monoculture growth for the organisms that constitute the communities
function [MonoResults,singleModelInfo] = monocultureGrowth(Communities,pairedModelInfo, options)

% modelName: Community model name
% model: Community model structure
% biomass: Final biomass concentration
% abundance: Final abundance of the community
% cs_cons: Carbon source consumed

%% initiate arrays
MonoResults = {}; singleModelInfo ={};
noOfPoints = 5; %feasible solutions across which growth and production rates are averaged

%% Extract the monoculture inputs
for i=1:options.orgCount
    if i==1
        singleModelInfo{i,1} = pairedModelInfo{i,2};
        model = Communities{1,i};
        Monoculture{i} = changeRxnBounds(model,model.rxns(strmatch(strcat(pairedModelInfo{i,4},'_'),model.rxns)), 0,'b');
    elseif i>1
        singleModelInfo{i,1} = pairedModelInfo{i-1,4};
        model = Communities{1,i-1};
        Monoculture{i} = changeRxnBounds(model,model.rxns(strmatch(strcat(pairedModelInfo{i-1,2},'_'),model.rxns)), 0,'b');
    end
    %% run dfba for the monocultures
    MonoResults{i} = dFBACom(Monoculture{i},options);
    
    [prdtID,prdtConc] = deal(zeros(length(options.Products),1));
    biomass = 0;cs_cons =0;
    result = MonoResults{i};
    startIndex = max(find(result.solnstat==1))-noOfPoints; %to find the final concentrations
    csID = find(ismember(options.mediumMets,options.carbonSource));
    
    %% compute final biomass, prdt, cs concentrations
    if startIndex > 0
        index = startIndex + find(result.solnstat(startIndex:startIndex-1+noOfPoints)==1);
        if i==1
            biomass = mean(result.biomassarr(1,index),2);
            prdtConc = mean(result.prdtFVAConc1(:,index),2);
            cs_cons = mean(options.initMedium(find(ismember(options.mediumMets,options.carbonSource)),1)-result.medium_nutrient(csID,index));
            
        elseif i>1
            biomass = mean(result.biomassarr(2,index),2);
            prdtConc = mean(result.prdtFVAConc2(:,index),2);
            cs_cons = mean(options.initMedium(find(ismember(options.mediumMets,options.carbonSource)),1)-result.medium_nutrient(csID,index));
        end
    end
    
    %% update results
    MonoResults{i}.modelName = singleModelInfo{i};
    MonoResults{i}.model = Monoculture{i};
    MonoResults{i}.biomass = biomass;
    MonoResults{i}.prdtConc = prdtConc;
    MonoResults{i}.cs_cons = cs_cons;
    
end
end