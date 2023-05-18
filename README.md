# SharedEvolutionSocialMicrobiomes
This project contains the code to simulate evolution under mutation-selection-migration and genetic drift. The model was used to predict the effect of bacterial transmission on the landscape of mutations shared among populations of E. coli colonizing the gut of mice (Frazao, N. and Gordo, I. 2023).
There are two models for the effects of mutations: 
one where mutations have a fixed effect MRSocialSingleSnew.cpp; 
another where mutations have exponentially distributed effects: 
SocialExpST400freq.cpp (that runs for 400 generations) or SocialExpST1500freq.cpp (that runs for 1500 generations)
Examples of scripts to run the codes are in scriptSocialExpS400.txt and scriptSocialExpS1500.txt
