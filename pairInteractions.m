%% finds the pairwise interactions for a two-member community
function [pairwiseInteractions] = pairInteractions(CommunityResults,CommunityModelInfo,MonoResults,singleModelInfo)
pairwiseInteractions={};
pairwiseInteractions{1, 1} = 'pairedModelID';
pairwiseInteractions{1, 2} = 'ModelID1';
pairwiseInteractions{1, 3} = 'ModelID2';
pairwiseInteractions{1, 4} = 'pairedGrowth_Model1';
pairwiseInteractions{1, 5} = 'pairedGrowth_Model2';
pairwiseInteractions{1, 6} = 'singleGrowth_Model1';
pairwiseInteractions{1, 7} = 'singleGrowth_Model2';
pairwiseInteractions{1, 8} = 'Interaction_Model1';
pairwiseInteractions{1, 9} = 'Interaction_Model2';
pairwiseInteractions{1, 10} = 'Interaction';
ratio=0.1; %percentage difference to categorize interaction

for i= 1:size(CommunityModelInfo,1)
    pairwiseInteractions{i+1, 1} = CommunityModelInfo{i,1};
    pairwiseInteractions{i+1, 2} = CommunityModelInfo{i,2};
    pairwiseInteractions{i+1, 3} = CommunityModelInfo{i,4};
    model1 = find(ismember(singleModelInfo,CommunityModelInfo{i,2}));
    model2 = find(ismember(singleModelInfo,CommunityModelInfo{i,4}));
    % combined growth
    pairwiseInteractions{i + 1, 4} = CommunityResults{1,i}.biomass(1,1);
    pairwiseInteractions{i + 1, 5} = CommunityResults{1,i}.biomass(2,1);
    
    pairwiseInteractions{i+1, 6} = MonoResults{1,model1}.biomass;
    pairwiseInteractions{i+1, 7} = MonoResults{1,model2}.biomass;
end

for i= 2:size(pairwiseInteractions,1)
    outcome1=0; outcome2=0;
    if pairwiseInteractions{i,4}<0.01 || pairwiseInteractions{i,5}<0.01 || pairwiseInteractions{i,6}<0.01 || pairwiseInteractions{i,7}<0.01
        pairwiseInteractions{i , 8} = 'Infeasible';
        pairwiseInteractions{i , 9} = 'Infeasible';
        pairwiseInteractions{i , 10} = 'Infeasible';
        outcome1=NaN; outcome2=NaN;
    end
    if pairwiseInteractions{i,4} > (1+ratio)*(pairwiseInteractions{i,6}/2)
        outcome1 = 1;
    elseif pairwiseInteractions{i,4} < (1-ratio)*(pairwiseInteractions{i,6}/2)
        outcome1 = -1;
    end
    if pairwiseInteractions{i,5} > (1+ratio)*(pairwiseInteractions{i,7}/2)
        outcome2 = 1;
    elseif pairwiseInteractions{i,5} < (1-ratio)*(pairwiseInteractions{i,7}/2)
        outcome2 = -1;
    end
    
    if outcome1 == 1 && outcome2 == 1
        pairwiseInteractions{i,8} = 'Mutualism';
        pairwiseInteractions{i,9} = 'Mutualism';
        pairwiseInteractions{i,10} = 'Mutualism';
    elseif outcome1 == 1 && outcome2 == -1
        pairwiseInteractions{i,8} = 'Parasitism_Taker';
        pairwiseInteractions{i,9} = 'Parasitism_Giver';
        pairwiseInteractions{i,10} = 'Parasitism';
    elseif outcome1 == -1 && outcome2 == 1
        pairwiseInteractions{i,8} = 'Parasitism_Giver';
        pairwiseInteractions{i,9} = 'Parasitism_Taker';
        pairwiseInteractions{i,10} = 'Parasitism';
    elseif outcome1 == -1 && outcome2 == -1
        pairwiseInteractions{i,8} = 'Competition';
        pairwiseInteractions{i,9} = 'Competition';
        pairwiseInteractions{i,10} = 'Competition';
    elseif outcome1 == 0 && outcome2 == 1
        pairwiseInteractions{i,8} = 'Commensalism_Unaffected';
        pairwiseInteractions{i,9} = 'Commensalism_Taker';
        pairwiseInteractions{i,10} = 'Commensalism';
    elseif outcome1 == 1 && outcome2 == 0
        pairwiseInteractions{i,8} = 'Commensalism_Taker';
        pairwiseInteractions{i,9} = 'Commensalism_Unaffected';
        pairwiseInteractions{i,10} = 'Commensalism';
    elseif outcome1 == 0 && outcome2 == 0
        pairwiseInteractions{i,8} = 'Neutralism';
        pairwiseInteractions{i,9} = 'Neutralism';
        pairwiseInteractions{i,10} = 'Neutralism';
    elseif outcome1 == 0 && outcome2 == -1
        pairwiseInteractions{i,8} = 'Amensalism_Unaffected';
        pairwiseInteractions{i,9} = 'Amensalism_Giver';
        pairwiseInteractions{i,10} = 'Amensalism';
    elseif outcome1 == -1 && outcome2 == 0
        pairwiseInteractions{i,8} = 'Amensalism_Giver';
        pairwiseInteractions{i,9} = 'Amensalism_Unaffected';
        pairwiseInteractions{i,10} = 'Amensalism';
    end
end
end
