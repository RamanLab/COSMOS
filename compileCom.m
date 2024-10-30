%% Compiles the community results for list of communities using abundance cutoff
function [CommunityResults,CommunityModelInfo] = compileCom(Communities,pairedModelInfo,options,allResults)
% medium_nutrient: change in medium at each time step
% biomassarr: change in biomass at each time step
% mu_est:growth rate at each time step
% solnstat: array of solution stat (returns 1 if feasible and 0 otherwise)
% bioAdded: biomass added due to growth at each time step
% timearr: array of timeSteps
% prdtFVAConc1: product concentration for organism A in community at each time step
% prdtFVAConc2: product concentration for organism B in community at each time step
% totPrdtFVAConc: total product concentration of community at each time step
% prdtID: Product reaction numbers
% prdtRxns: Product reaction names
% modelName: Community model name
% model: Community model structure
% biomass: Final biomass concentration
% abundance: Final abundance of the community
% cs_cons: Carbon source consumed
% productConc1: Final product concentration secreted by organism A
% productConc2: Final product concentration secreted by organism B
% totPrdtConc: Final product concentration secreted by the community

%% initialize variables
temp=0; CommunityResults = {}; CommunityModelInfo = {};
% mean of solutions across a timeframe
noOfPoints = 5; %feasible solutions across which growth and production rates are averaged

%% compilation of final values of biomass, abundance, etc
for i = 1:length(allResults)
    %Compile each community result
    result = allResults{1,i}; Community = Communities{1,i};
    biomass = [0;0]; abundance =[0;0]; cs_cons=[];
    %Products to be analyzed
    [prdtConc1, prdtConc2, totPrdtConc] = deal(zeros(length(options.Products),1));
    
    %Analyze data for a set of datapoints at stationary phase
    startIndex = max(find(result.solnstat==1))-noOfPoints; %to find the final concentrations
    if startIndex > 0 %to ensure enough datapoints
        index = startIndex + find(result.solnstat(startIndex:startIndex-1+noOfPoints)==1);
        %mean biomass and abundance across datapoints
        biomass = mean(result.biomassarr(:,index),2);
        abundance = biomass/sum(biomass,1);
        %for minimum growth
        if min(abundance)>(options.abdCutoff/100) && min(biomass)> (min(options.initBiomass)+(options.abdCutoff/100*min(options.initBiomass))) %more than x percent of maximum values
            temp = temp+1;
            %Carbon source consumed
            csID = find(ismember(options.mediumMets,options.carbonSource));
            cs_cons = mean(options.initMedium(find(ismember(options.mediumMets,options.carbonSource)),1)-result.medium_nutrient(csID,index));
            
            %Product concentration
            prdtConc1 = mean(result.prdtFVAConc1(:,index),2);
            prdtConc2 = mean(result.prdtFVAConc2(:,index),2);
            totPrdtConc = mean(result.totPrdtFVAConc(:,index),2);
            %% update results
            result.modelName = pairedModelInfo{i,1};
            result.model = Community;
            
            result.biomass = biomass;
            result.abundance = abundance;
            result.cs_cons = cs_cons;
            
            result.productConc1 = prdtConc1;
            result.productConc2 = prdtConc2;
            result.totPrdtConc = totPrdtConc;
            
            if temp>0
                CommunityResults{1,temp} = result;
                CommunityModelInfo(temp,:) = pairedModelInfo(i,:);
            end
        end
    end
    
end
end
