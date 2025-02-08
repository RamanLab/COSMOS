%to collate the data from all four environments
%collate carbon source analysis results of all the environments under a nested cell array CSEnv
load('csAerRich.mat','CSResults');
CSEnv{1} = CSResults;
load('csAerMin.mat','CSResults');
CSEnv{2} = CSResults;
load('csAnaerRich.mat','CSResults');
CSEnv{3} = CSResults;
load('csAnaerMin.mat','CSResults');
CSEnv{4} = CSResults;
% Headers for all sheets
Header = {'Environment','CarbonSource','Community','Product','Productivity_Ratio'};
Environments = {'Aer_Rich','Aer_Min','Anaer_Rich','Anaer_Min'};
Data={}; tempA=0;
%collate data according to header
for i = 1:length(Environments)
    CSResults = CSEnv{1,i};
    for j=1:length(CSResults)
        PrdtAnalysis = CSResults{1,j}.PrdtAnalysis;
        options = CSResults{1,j}.options;
        for k = 1:length(PrdtAnalysis)
            Data(tempA+1 : tempA+length(options.ProductName),5) = num2cell(PrdtAnalysis{1,k}.prdtRatio);
            Data(tempA+1 : tempA+length(options.ProductName),4) = options.ProductName;
            Data(tempA+1 : tempA+length(options.ProductName),3) = {PrdtAnalysis{1,k}.modelName};
            Data(tempA+1 : tempA+length(options.ProductName),2) = (options.carbonSource);
            Data(tempA+1 : tempA+length(options.ProductName),1) = Environments(i);
            
            tempA = tempA+length(options.ProductName);
        end
    end
end
for i=1:length(Data(:,2))   
    Data{i, 2} = strrep(Data{i, 2}, 'glc__D_e[u]', 'Glucose');
    Data{i, 2} = strrep(Data{i, 2}, 'fru_e[u]', 'Fructose');
    Data{i, 2} = strrep(Data{i, 2}, 'sucr_e[u]', 'Sucrose');
    Data{i, 2} = strrep(Data{i, 2}, 'malt_e[u]', 'Maltose');
    Data{i, 2} = strrep(Data{i, 2}, 'lac__L_e[u]', 'Lactose');
    Data{i, 2} = strrep(Data{i, 2}, 'xyl__D_e[u]', 'Xylose');
    Data{i, 2} = strrep(Data{i, 2}, 'gal_e[u]', 'Galactose');
end
%to write the data
 xlswrite('CSData.xls', [Header;Data]);
