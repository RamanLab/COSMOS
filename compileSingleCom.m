%% Compile the community results for a single community
%% gives biomass and abundance without using the abundance cutoff
function [result] = compileSingleCom(Community,options,result)

% mean of solutions across a timeframe
noOfPoints = 5; %time interval across which growth and production rates are averaged

%Compile each community result
biomass = [0;0]; abundance =[0;0]; cs_cons =[];
%Products to be analyzed
[prdtConc1, prdtConc2, totPrdtConc] = deal(zeros(length(options.Products),1));


%Analyze data for a set of datapoints at stationary phase
startIndex = max(find(result.solnstat==1))-noOfPoints;
if startIndex > 0 %to ensure enough datapoints
    index = startIndex + find(result.solnstat(startIndex:startIndex-1+noOfPoints)==1);
    
    %mean biomass and abundance across datapoints
    biomass = mean(result.biomassarr(:,index),2);
    abundance = biomass/sum(biomass,1);
    
    csID = find(ismember(options.mediumMets,options.carbonSource));
    cs_cons = mean(options.initMedium(find(ismember(options.mediumMets,options.carbonSource)),1)-result.medium_nutrient(csID,index));
    
    %Product concentration
    prdtConc1 = mean(result.prdtFVAConc1(:,index),2);
    prdtConc2 = mean(result.prdtFVAConc2(:,index),2);
    totPrdtConc = mean(result.totPrdtFVAConc(:,index),2);
end

result.biomass = biomass;
result.abundance = abundance;
result.cs_cons = cs_cons;

result.productConc1 = prdtConc1;
result.productConc2 = prdtConc2;
result.totPrdtConc = totPrdtConc;
end

