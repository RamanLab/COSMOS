function writeResultsToExcel(filename,data_filename)
    % Function to write microbial system analysis results to an Excel file
    % Inputs:
    %   filename - Name of the output Excel file (e.g., 'results.xls')

    % Load the required variables from aerRich.mat
    results = load(data_filename); 
    
    % Headers for different sheets
    TopOrgs_header = {'Microbial system', 'Number of products for which it has maximum productivity'};
    Productivity_header = {'Community Name', 'Productivity', 'Abundance A', 'Abundance B'};

    % Write the data using writecell (recommended over xlswrite)
    writecell(results.pairwiseInteractions, filename, 'Sheet', 'Communities');
    writecell(results.interactionGeneric, filename, 'Sheet', 'InteractionGeneric');
    writecell(results.interactionSpecific, filename, 'Sheet', 'InteractionSpecific');
    writecell([TopOrgs_header; results.BestSystem.TopOrgs], filename, 'Sheet', 'TopOrganisms');

    % Split abundance into two columns and write per product
    for i = 1:length(results.options.ProductName)
        BestSystemMod = results.BestSystem.ProductwiseResults;
        for j = 1:size(results.BestSystem.ProductwiseResults{1, i}, 1)
            abd_combined = results.BestSystem.ProductwiseResults{1, i}{j, 3};
            for k = 1:length(abd_combined)
                BestSystemMod{1, i}{j, 2 + k} = abd_combined(k);
                if k == 1
                    BestSystemMod{1, i}{j, 3 + k} = 0;
                end
            end
        end
        writecell([Productivity_header; BestSystemMod{1, i}], filename, 'Sheet', results.options.ProductName{i});
    end
end
