%% to write the results into excel

%Headers for all sheets

TopOrgs_header = {'Microbial system', 'Number of products for which it has maximum productivity'};
Productivity_header = {'Community Name','Productivity','Abundance A','Abundance B'};

%write the data
xlswrite('results.xls', pairwiseInteractions, 'Communities');

xlswrite('results.xls', InteractionGeneric, 'InteractionGeneric');

xlswrite('results.xls', InteractionSpecific, 'InteractionSpecific');

xlswrite('results.xls', [TopOrgs_header; BestSystem.TopOrgs], 'TopOrganisms');

%split abundance to two columns
for i=1:length(options.ProductName)
    BestSystemMod = BestSystem.ProductwiseResults;
    for j=1:size(BestSystem.ProductwiseResults{1,i},1)
        abd_combined = BestSystem.ProductwiseResults{1,i}{j,3};
        for k = 1:length(abd_combined)
            BestSystemMod{1,i}{j,2+k} = abd_combined(k);
            if k==1
                BestSystemMod{1,i}{j,3+k} = 0; 
            end
        end
    end

    xlswrite('results.xls', [Productivity_header; BestSystemMod{1,i}], options.ProductName{i});  

end