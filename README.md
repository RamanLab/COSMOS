# MonoVsCocultures
Compares monocultures and cocultres in an identical environment to identify the optimal microbial system. The comparison is done across diverse environments anf the effect of various factors such as medium composition, community interaction, carbon source etc are analysed. 

Masterscripts for the four environments – aerobicRich.m, aerobicMinimal.m, anaerobicRich.m, anaerobicMinimal.m

Effect of interactions – InteractionAnalysis.m

Effect of carbon source – CSanalysis.m

Dynamin FBA – dFBACom.m

Comparison based on productivity – ProductivityAnalysis.m

Compare based on yield – YieldAnalysis.m

Optimize biomass ratio – optBiomassRatio.m

Validation with experimental study – validation_SoKp.m

All models used in this analysis were sourced from the BiGG Models and BioModels databases. Some annotations were adjusted to maintain consistency across the dataset.

Prerequisites

All simulations were performed in MATLAB R2018a (MathWorks Inc., USA) using:

    COBRA Toolbox v3.0
    IBM ILOG CPLEX 12.8

Research Article Link: 

All the analysis were performed on MATLAB R2018a and COBRAToolbox v.3.0 and IBM Cplex 12.8 solver.
