%% Yield of Community vs Monoculture analysis
function [PrdtAnalysis,BestSystem] = YieldAnalysis(CommunityResults,CommunityModelInfo,MonoResults,singleModelInfo,options)

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
    
    for j=1:length(options.Products)
        if MonoResults{1,model1}.cs_cons > 0.001
            monoprdt1(j) = MonoResults{1,model1}.prdtConc(j)/ MonoResults{1,model1}.cs_cons;
        else
            monoprdt1(j) = NaN;
        end
        if MonoResults{1,model2}.cs_cons > 0.001
            monoprdt2(j) = MonoResults{1,model2}.prdtConc(j)/ MonoResults{1,model2}.cs_cons;
        else
            monoprdt2(j) = NaN;
        end
        if CommunityResults{1,i}.cs_cons > 0.001
            commprdt(j) = CommunityResults{1,i}.totPrdtConc(j)/CommunityResults{1,i}.cs_cons;
            commprdt1(j) = CommunityResults{1,i}.productConc1(j)/CommunityResults{1,i}.cs_cons;
            commprdt2(j) = CommunityResults{1,i}.productConc2(j)/CommunityResults{1,i}.cs_cons;
        else
            commprdt(j) = NaN; commprdt1(j) = NaN; commprdt2(j) = NaN;
        end
        
        %boolean change and ratio of change in product concentration between
        %community and respective monocultures
        if abs(monoprdt1(j)+monoprdt2(j)) < 0.0001 && abs(commprdt(j)) > 0.0001
            prdtRatio(j) = 1000;
            prdtChange(j) = 1;
        elseif abs(monoprdt1(j)+monoprdt2(j)) > 0.0001 && abs(commprdt(j)) < 0.0001
            prdtRatio(j) = ((commprdt(j)-max(monoprdt1(j),monoprdt2(j)))/abs(max(monoprdt1(j),monoprdt2(j))));                      prdtChange(j) = -1;
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
        [maxPrdt(j),index] = max([monoprdt1(j),monoprdt2(j),commprdt(j)]);
        if index==1
            OptSys{j} = MonoResults{1,model1}.modelName;
        elseif index==2
            OptSys{j} = MonoResults{1,model2}.modelName;
        elseif index==3
            OptSys{j} = CommunityResults{1,i}.modelName;
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
end

for i=1:length(PrdtAnalysis)
    for j=1:length(options.Products)
        maxPrdtList(i,j)= PrdtAnalysis{1,i}.maxPrdt(j);
        OptSysListOld{i,j} = PrdtAnalysis{1,i}.OptSys{j};
    end
end

for j=1:length(options.Products)
    [sorted_maxPrdtConcList(:,j),idx] = sort(maxPrdtList(:,j),'descend');
    OptSysList(:,j) = OptSysListOld(idx,j);
end

BestSystem.AllPrdtsConc = vertcat(options.ProductName',num2cell(sorted_maxPrdtConcList));
BestSystem.AllPrdtsSysName = vertcat(options.ProductName',OptSysList);

%to give a unique value table
for j=1:length(options.Products)
    [ProductwiseResults{1,j},index] = unique(sorted_maxPrdtConcList(:,j),'stable');
    for i=1:length(ProductwiseResults{1,j})
        if ProductwiseResults{1,j}(i,1) ==0
            index(i)=[];
            BestSystem.ProductwiseResults{:,j} = horzcat(OptSysList(index,j),num2cell(ProductwiseResults{1,j}(1:end-1,1)));
        else
            BestSystem.ProductwiseResults{:,j} = horzcat(OptSysList(index,j),num2cell(ProductwiseResults{1,j}));
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

[sorted_TopOrgCount,idx] = sort(TopOrgCount,'descend');
sorted_TopOrgName = TopOrgName(idx);
BestSystem.TopOrgs = horzcat(sorted_TopOrgName,num2cell(sorted_TopOrgCount));

end

