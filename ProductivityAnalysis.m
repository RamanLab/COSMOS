%% Productivity of Community vs Monoculture analysis
function [PrdtAnalysis,BestSystem] = ProductivityAnalysis(CommunityResults,CommunityModelInfo,MonoResults,singleModelInfo,options)

%commprdt: productivity in community
%commprdt1: productivity in organism A in community
%commprdt2: productivity in organism B in community
%monoprdt1: productivity in organism A in monoculture
%monoprdt2: productivity in organism B in monoculture
%monoprdt: productivity in monoculture

%prdtChange :boolean change in product concentration between community and monocultures
%prdtChange1 :boolean change in product concentration between community and first monoculture
%prdtChange2 :boolean change in product concentration between community and second monoculture
%prdtRatio :ratio of change in product concentration between community and monocultures

%maxPrdt: maximum productivity acheived
%OptSys: Microbial system with the maximum productivity
%abdSys: Abundance of the optimal system

noOfPoints=5; %number of feasible solns approximated
ratio = 0.1; %minimum ratio of increase in production
PrdtAnalysis={};OptSys={};abdSys ={};

for i=1:size(CommunityResults,2)
    [prdtChange,prdtRatio,prdtChange1,prdtChange2] = deal(zeros(length(options.Products),1));
    PrdtAnalysis{i}.modelName = CommunityResults{1,i}.modelName;
    PrdtAnalysis{i}.Products = options.Products;
    model1 = find(ismember(singleModelInfo,CommunityModelInfo{i,2}));
    model2 = find(ismember(singleModelInfo,CommunityModelInfo{i,4}));
    
    startIndex1 = max(find(CommunityResults{1,i}.solnstat==1))-noOfPoints;
    if startIndex1 > 0 %to ensure enough datapoints
        index1 = startIndex1 + find(CommunityResults{1,i}.solnstat(startIndex1:startIndex1-1+noOfPoints)==1);
        CommTime = mean(CommunityResults{1,i}.timearr(:,index1),2);
    end
    startIndex2 = max(find(MonoResults{1,model1}.solnstat==1))-noOfPoints;
    if startIndex2 > 0 %to ensure enough datapoints
        index2 = startIndex2 + find(MonoResults{1,model1}.solnstat(startIndex2:startIndex2-1+noOfPoints)==1);
        Mono1Time = mean(MonoResults{1,model1}.timearr(:,index2),2);
    end
    startIndex3 = max(find(MonoResults{1,model2}.solnstat==1))-noOfPoints;
    if startIndex3 > 0 %to ensure enough datapoints
        index3 = startIndex3 + find(MonoResults{1,model2}.solnstat(startIndex3:startIndex3-1+noOfPoints)==1);
        Mono2Time = mean(MonoResults{1,model2}.timearr(:,index3),2);
    end
    
    for j=1:length(options.Products)
        if Mono1Time > 0.001
            monoprdt1(j) = MonoResults{1,model1}.prdtConc(j)/ Mono1Time;
        else
            monoprdt1(j) = NaN;
        end
        if Mono2Time > 0.001
            monoprdt2(j) = MonoResults{1,model2}.prdtConc(j)/ Mono2Time;
        else
            monoprdt2(j) = NaN;
        end
        if CommTime > 0.001
            commprdt(j) = CommunityResults{1,i}.totPrdtConc(j)/CommTime;
            commprdt1(j) = CommunityResults{1,i}.productConc1(j)/CommTime;
            commprdt2(j) = CommunityResults{1,i}.productConc2(j)/CommTime;
        else
            commprdt(j) = NaN; commprdt1(j) = NaN; commprdt2(j) = NaN;
        end
        
        %boolean change and ratio of change in product concentration between
        %community and respective monocultures
        if abs(monoprdt1(j)+monoprdt2(j)) < 0.0001 && abs(commprdt(j)) > 0.0001
            prdtRatio(j) = 1000;
            prdtChange(j) = 1;
        elseif abs(monoprdt1(j)+monoprdt2(j)) > 0.0001 && abs(commprdt(j)) < 0.0001
            prdtRatio(j) = ((commprdt(j)-max(monoprdt1(j),monoprdt2(j)))/abs(max(monoprdt1(j),monoprdt2(j))));
            prdtChange(j) = -1;
        elseif abs(monoprdt1(j)+monoprdt2(j)) > 0.0001 && abs(commprdt(j)) > 0.0001
            prdtRatio(j) = ((commprdt(j)-max(monoprdt1(j),monoprdt2(j)))/abs(max(monoprdt1(j),monoprdt2(j))));
            if commprdt(j) > (1+ratio)*max(monoprdt1(j),monoprdt2(j))
                prdtChange(j) = 1;
            elseif commprdt(j) < (1-ratio)*max(monoprdt1(j),monoprdt2(j))
                prdtChange(j) = -1;
            end
        elseif abs(monoprdt1(j)+monoprdt2(j)) < 0.0001 && abs(commprdt(j)) < 0.0001
            prdtChange(j) = NaN;
        end
        
        %boolean change and ratio of change in product concentration between
        %each organism in the community and respective monoculture
        if commprdt1(j) > (1+ratio)*monoprdt1(j)
            prdtChange1(j) = 1;
        elseif commprdt1(j) < (1-ratio)*monoprdt1(j)
            prdtChange1(j) = -1;
        elseif abs(commprdt1(j))<0.0001 && abs(monoprdt1(j))<0.0001
            prdtChange1(j) = NaN;
        end
        
        if commprdt2(j) > (1+ratio)*monoprdt2(j)
            prdtChange2(j) = 1;
        elseif commprdt2(j) < (1-ratio)*monoprdt2(j)
            prdtChange2(j) = -1;
        elseif abs(commprdt2(j))<0.0001 && abs(monoprdt2(j))<0.0001
            prdtChange2(j) = NaN;
        end
        
        %to find the system that produces maximum product
        [maxPrdt(j),index1] = max([monoprdt1(j),monoprdt2(j),commprdt(j)]);
        if index1==1
            OptSys{j} = MonoResults{1,model1}.modelName;
            abdSys{j} = 1;
        elseif index1==2
            OptSys{j} = MonoResults{1,model2}.modelName;
            abdSys{j} = 1;
        elseif index1==3
            OptSys{j} = CommunityResults{1,i}.modelName;
            abdSys{j} = CommunityResults{1,i}.abundance;
        end
        
    end
      
    PrdtAnalysis{i}.commprdt = commprdt';
    PrdtAnalysis{i}.commprdt1 = commprdt1';
    PrdtAnalysis{i}.commprdt2 = commprdt2';
    
    PrdtAnalysis{i}.monoprdt1 = monoprdt1';
    PrdtAnalysis{i}.monoprdt2 = monoprdt2';
    PrdtAnalysis{i}.monoprdt = max(monoprdt1',monoprdt2');
    
    PrdtAnalysis{i}.prdtChange = prdtChange;
    PrdtAnalysis{i}.prdtChange1 = prdtChange1;
    PrdtAnalysis{i}.prdtChange2 = prdtChange2;
    PrdtAnalysis{i}.prdtRatio = prdtRatio;
    PrdtAnalysis{i}.maxPrdt = maxPrdt';
    PrdtAnalysis{i}.OptSys = OptSys';
    PrdtAnalysis{i}.abdSys = abdSys';
end

for i=1:length(PrdtAnalysis)
    for j=1:length(options.Products)
        maxPrdtList(i,j)= PrdtAnalysis{1,i}.maxPrdt(j);
        OptSysListOld{i,j} = PrdtAnalysis{1,i}.OptSys{j};
        abundanceList{i,j} = PrdtAnalysis{1,i}.abdSys{j};
    end
end

for j=1:length(options.Products)
    [sorted_maxPrdtConcList(:,j),idx1] = sort(maxPrdtList(:,j),'descend');
    OptSysList(:,j) = OptSysListOld(idx1,j);
    abundanceListNew(:,j) = abundanceList(idx1,j);
end

BestSystem.AllPrdtsConc = vertcat(options.ProductName',num2cell(sorted_maxPrdtConcList));
BestSystem.AllPrdtsSysName = vertcat(options.ProductName',OptSysList);

%to give a unique value table
for j=1:length(options.Products)
    [ProductwiseResults{1,j},index2] = unique(sorted_maxPrdtConcList(:,j),'stable');
    for i=1:length(index2)
        if ProductwiseResults{1,j}(i,1) <= 0.0001
            [ProductwiseResults{1,j},index2] = unique(sorted_maxPrdtConcList(:,j),'stable');
            index2(i)=[];
            BestSystem.ProductwiseResults{:,j} = horzcat(OptSysList(index2,j),num2cell(ProductwiseResults{1,j}(1:end-1,1)),abundanceListNew(index2,j));
        else
            BestSystem.ProductwiseResults{:,j} = horzcat(OptSysList(index2,j),num2cell(ProductwiseResults{1,j}),abundanceListNew(index2,j));
        end
    end
end

%to find top performing systems/organisms across products
TotSysList={};

for i=1:length(BestSystem.ProductwiseResults)
    TotSysList = vertcat(TotSysList,BestSystem.ProductwiseResults{1,i}(:,1));
end

c = categorical(TotSysList);
TopOrgName = categories(c);
TopOrgCount = countcats(c);

[sorted_TopOrgCount,idx2] = sort(TopOrgCount,'descend');
sorted_TopOrgName = TopOrgName(idx2);
BestSystem.TopOrgs = horzcat(sorted_TopOrgName,num2cell(sorted_TopOrgCount));
end
