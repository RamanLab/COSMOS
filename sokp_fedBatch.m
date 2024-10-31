%% to simulate fed batch in the S. oneidensis- K.pneumoniae community
function [result] = sokp_fedBatch(model,options)

%% to initiate necessary arrays
result = struct;
minConc = 0.001; %the minimal detectable concentration in medium
temp = 0; %iteration number through timeSteps

biomass = zeros(length(model.infoCom.spName),1);
prdtID = zeros(length(options.Products),length(options.initBiomass));
prdtRxns= cell(length(options.Products),length(options.initBiomass));
%to monitor how the medium and biomas etc change over time
%array of nutrient concentrations across timeSteps
medium_nutrient = zeros(length(options.initMedium),round(options.maxTime/options.delt,0));
medium = options.initMedium; %to set initial state of medium
medium_nutrient(:,1) = medium;
%array of biomass, added biomass and growth rate across timeSteps
[biomassarr,biomass_added,mu_est] = deal(zeros(length(options.initBiomass),round(options.maxTime/options.delt,0)));
biomassarr(:,1) = biomass;
%other arrays
[timearr,statarr] = deal(zeros(1,round(options.maxTime/options.delt,0)));
infeasibleCount=0;
%initiate biomass
biomass(:) = options.initBiomass;
%to estimate concentration of products
[prdtFVAConc1,prdtFVAConc2]= deal(zeros(length(options.Products),round(options.maxTime/options.delt,0)));

%% dynamic FBA
for timeStep = 1: round(options.maxTime/options.delt,0)
    tcurr = timeStep*options.delt;
    temp = temp+1; %iteration number
    timearr(temp) = tcurr;
    %exchange rxns of each community
    for i=1:length(biomass)
        mediumExcRxns{:,i} = intersect(model.indCom.EXsp(:,i),find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets))));
    end
    %estimating medium consumption
    for i = 1:length(options.mediumMets)
        uptake = computeUptake(options.mediumMets(i), medium_nutrient(i), mediumExcRxns(1,:), model, biomass, options.Vmax(i), options.Km(i), options.delt);
        if temp>1
            uptake = computeUptake(options.mediumMets(i), medium_nutrient(i,temp-1), mediumExcRxns(1,:), model, biomass, options.Vmax(i), options.Km(i), options.delt);
        end
        for j = 1:length(biomass)
            mediumExchange = intersect(model.indCom.EXsp(:,j),find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets(i)))));
            %update model bounds with computed uptake
            if ~isempty(mediumExchange) && ~isempty(uptake) && length(uptake)==length(biomass)
                model.lb(mediumExchange) = - uptake(j)/(options.delt*biomass(j)); %uptake fluxes are negative
            end
        end
    end
    
    %standard FBA of model
    solution = optimizeCbModel(model,'max','one');
    %update biomass and growth rates
    for i=1:length(biomass)
        statarr(temp) = solution.stat;
        if biomass(i)==0 && solution.stat==1
            biomass_added(i,temp) = 0;
        else
            if solution.stat ==1 && infeasibleCount<=5
                biomass_added(i,temp) = solution.x(model.indCom.spBm(i))*biomass(i)*options.delt;
                mu_est(i,temp) = biomass_added(i,temp)/(biomass(i)*options.delt);
                
                %to perform FVA for products
                for k=1:length(options.Products)
                    maxFlux1=[]; maxFlux2=[];
                    if ~isempty(find(contains(model.infoCom.EXsp(:,i),options.Products(k))==1,1))
                        prdtID(k,i) = (find(contains(model.infoCom.EXsp(:,i),options.Products(k))==1));
                    end
                    if i==1 && prdtID(k,i)~=0
                        prdtRxns(k,1) = model.infoCom.EXsp(prdtID(k,1),1);
                        [~, maxFlux1] = fastFVA(model, 90, 'max','ibm_cplex',prdtRxns(k,1));
                        if ~isempty(maxFlux1) && ~isnan(maxFlux1)
                            prdtFVAConc1(k,temp) = maxFlux1*biomass(i)*options.delt;
                        end
                    end
                    if i==2 && prdtID(k,i)~=0
                        prdtRxns(k,2) = model.infoCom.EXsp(prdtID(k,2),2);
                        [~, maxFlux2] = fastFVA(model, 90, 'max','ibm_cplex',prdtRxns(k,2));
                        if ~isempty(maxFlux2) && ~isnan(maxFlux2)
                            prdtFVAConc2(k,temp) = maxFlux2*biomass(i)*options.delt;
                        end
                    end
                end
                
            else
                biomass_added(i,temp)=0;
                infeasibleCount = infeasibleCount + 1;
                mu_est(i,temp) = NaN;
            end
        end
        %update all arrays
        biomass(i) = biomass(i) + biomass_added(i,temp);
        biomassarr(i,temp) = biomass(i);
    end
    
    %to update the concentration of all media components
    for m =1:length(options.mediumMets)
        for i=1:length(biomass)
            %only update media components for species that are actually present
            if biomass(i)>0  && solution.stat ==1
                %in case a specific model does not have metabolite m
                metIndex = intersect(model.indCom.EXcom,find(ismember(model.rxns,findRxnsFromMets(model,options.mediumMets(m)))));
                if solution.stat == 1 && ~isempty(metIndex)
                    medium(m) = medium(m) + solution.x(metIndex)*biomass(i)*options.delt;
                end
            end
        end
    end
    
    %set all medium components that have fallen below a threshold value to zero, to avoid solver problems with very small values
    for m = 1: length(options.mediumMets)
        if medium(m) <= minConc
            medium(m)=0;
        end
        %% nutrient feed in fed-batch
        if medium(1)< 108.6
            medium(1) = 325.7;
        end
        
    end
    
    %change in media components and biomass with time
    medium_nutrient(:,temp) = medium;
    
end

%% update arrays
result.medium_nutrient = medium_nutrient;
result.biomassarr = biomassarr;
result.mu_est = mu_est;  result.solnstat = statarr;
result.bioAdded = biomass_added;
result.timearr = timearr;
result.prdtFVAConc1 = prdtFVAConc1; result.prdtFVAConc2 = prdtFVAConc2;
result.totPrdtFVAConc = prdtFVAConc1 + prdtFVAConc2;
result.prdtID = prdtID; result.prdtRxns = prdtRxns;
result.modelName = model.description;

end
