%to collate the data from all four environments
%collate PrdtAnalysis of all the environments under a nested cell array PrdtAnalysisEnv
load('aerRich.mat','PrdtAnalysis','pairwiseInteractions');
PrdtAnalysisEnv{1} = PrdtAnalysis;
interactionEnv{1} = pairwiseInteractions;
load('aerMin.mat','PrdtAnalysis','pairwiseInteractions');
PrdtAnalysisEnv{2} = PrdtAnalysis;
interactionEnv{2} = pairwiseInteractions;
load('anaerRich.mat','PrdtAnalysis','pairwiseInteractions');
PrdtAnalysisEnv{3} = PrdtAnalysis;
interactionEnv{3} = pairwiseInteractions;
load('anaerMin.mat','PrdtAnalysis','pairwiseInteractions','options');
PrdtAnalysisEnv{4} = PrdtAnalysis;
interactionEnv{4} = pairwiseInteractions;
% Headers for all sheets
Header = {'Environment','Community','InteractionA','InteractionB','Product','Productivity_Ratio','ProductivityA','ProductivityB'};
Environments = {'Aer_Rich','Aer_Min','Anaer_Rich','Anaer_Min'};
Data={}; tempA=0; 
%collate data according to header
for i = 1:length(Environments)
    PrdtAnalysis = PrdtAnalysisEnv{1,i}; 
    interactions = interactionEnv{1,i};
    for j = 1:length(PrdtAnalysis)
        
        Data(tempA+1 : tempA+length(options.ProductName),6) = num2cell(PrdtAnalysis{1,j}.prdtRatio);
        Data(tempA+1 : tempA+length(options.ProductName),5) = options.ProductName;
        Data(tempA+1 : tempA+length(options.ProductName),2) = {PrdtAnalysis{1,j}.modelName};
        Data(tempA+1 : tempA+length(options.ProductName),1) = Environments(i); 
        Data(tempA+1 : tempA+length(options.ProductName),7) = num2cell(PrdtAnalysis{1,j}.commprdt1);
        Data(tempA+1 : tempA+length(options.ProductName),8) = num2cell(PrdtAnalysis{1,j}.commprdt2);
        
        for k = tempA+1:tempA+length(options.ProductName)
            [Lia,Locb] = ismember(interactions(:,1),Data(k,2));
            Data(k,3) = interactions(find(Locb==1),8);
            Data(k,4) = interactions(find(Locb==1),9);
        end
        tempA = tempA+length(options.ProductName);
    end
end
% to write the data
 writecell([Header;Data], 'TotalData.xls', 'Sheet','Data_AllEnv');
