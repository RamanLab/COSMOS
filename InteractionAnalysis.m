%% effect of interaction on production
function [interactionAnalysisA,interactionAnalysisB] = InteractionAnalysis(PrdtAnalysis,pairwiseInteractions)
%interactionAnalysisA - The count of generic or common interaction types
%interactionAnalysisB - The count of specific interaction types (giver and
%taker)

interactionAnalysisA = cell(7,5);
interactionAnalysisA(2:end,1) = {'Amensalism';'Commensalism';'Competition';'Mutualism';'Neutralism';'Parasitism'};
interactionAnalysisA(1,1:end) = {'Interaction_Type','Increase','Decrease','Nochange','Not produced'};
interactionAnalysisA(2:end,2:end) = {0};

interactionAnalysisB = cell(10,5);
interactionAnalysisB(2:end,1) = {'Amensalism_Giver';'Amensalism_Unaffected';'Commensalism_Unaffected';'Commensalism_Taker';'Competition';'Mutualism';'Neutralism';'Parasitism_Giver';'Parasitism_Taker'};
interactionAnalysisB(1,1:end) = {'Interaction_Type','Increase','Decrease','Nochange','Not produced'};
interactionAnalysisB(2:end,2:end) = {0};

if ~isempty(PrdtAnalysis)
    for i = 1:length(PrdtAnalysis{1,1}.Products)
        for j=1:size(pairwiseInteractions,1)-1
            [Lia,Locb] = ismember(interactionAnalysisA(2:end,1),pairwiseInteractions(j+1,10));
            if nnz(Locb)~=0
                k = find(Locb==1)+1;
                if PrdtAnalysis{1,j}.prdtChange(i) == 1
                    interactionAnalysisA{k,2} = interactionAnalysisA{k,2}+1;
                elseif PrdtAnalysis{1,j}.prdtChange(i) == -1
                    interactionAnalysisA{k,3} = interactionAnalysisA{k,3}+1;
                elseif PrdtAnalysis{1,j}.prdtChange(i) == 0
                    interactionAnalysisA{k,4} = interactionAnalysisA{k,4}+1;
                elseif isnan(PrdtAnalysis{1,j}.prdtChange(i))
                    interactionAnalysisA{k,5} = interactionAnalysisA{k,5}+1;
                end
            end
        end
    end
    
    for i = 1:length(PrdtAnalysis{1,1}.Products)
        for j=1:size(pairwiseInteractions,1)-1
            [Lia,Locb] = ismember(interactionAnalysisB(2:end,1),pairwiseInteractions(j+1,8));
            if nnz(Locb)~=0
                k = find(Locb==1)+1;
                if PrdtAnalysis{1,j}.prdtChange1(i) == 1
                    interactionAnalysisB{k,2} = interactionAnalysisB{k,2}+1;
                elseif PrdtAnalysis{1,j}.prdtChange1(i) == -1
                    interactionAnalysisB{k,3} = interactionAnalysisB{k,3}+1;
                elseif PrdtAnalysis{1,j}.prdtChange1(i) == 0
                    interactionAnalysisB{k,4} = interactionAnalysisB{k,4}+1;
                elseif isnan(PrdtAnalysis{1,j}.prdtChange1(i))
                    interactionAnalysisB{k,5} = interactionAnalysisB{k,5}+1;
                end
            end
        end
    end
    
    for i = 1:length(PrdtAnalysis{1,1}.Products)
        for j=1:size(pairwiseInteractions,1)-1
            [Lia,Locb] = ismember(interactionAnalysisB(2:end,1),pairwiseInteractions(j+1,9));
            if nnz(Locb)~=0
                k = find(Locb==1)+1;
                if PrdtAnalysis{1,j}.prdtChange2(i) == 1
                    interactionAnalysisB{k,2} = interactionAnalysisB{k,2}+1;
                elseif PrdtAnalysis{1,j}.prdtChange2(i) == -1
                    interactionAnalysisB{k,3} = interactionAnalysisB{k,3}+1;
                elseif PrdtAnalysis{1,j}.prdtChange2(i) == 0
                    interactionAnalysisB{k,4} = interactionAnalysisB{k,4}+1;
                elseif isnan(PrdtAnalysis{1,j}.prdtChange2(i))
                    interactionAnalysisB{k,5} = interactionAnalysisB{k,5}+1;
                end
            end
        end
    end
end
end
