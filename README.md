# Estimating Distance-to-Default with a Sector-Specific Liability Adjustment via Sequential Monte Carlo

Name of QuantLet: DTD_SMC

Published in: Applied Quantitative Finance third Edition

Description: Distance-to-Default (DTD), a widely adopted corporate default predictor, arises from the classical structural credit risk model of Merton (1974). The modern way of estimating DTD applies the model on an observed time series of equity values along with the default point definition made popular by the commercial KMV model. It is meant to be a default trigger level one year from the evaluation time, and is assumed to be the short-term debt plus 50% of the long-term debt. This default point assumption, however, leaves out other corporate liabilities, which can be substantial and particularly so for financial firms. Duan, et al (2012) rectified it by adding other liabilities after applying an unknown but estimable haircut. We assume a common haircut for all firms in a sector and devise a novel density-tempered expanding-data sequential Monte Carlo method to jointly estimate this common and other firm-specific parameters. Joint estimation is challenging due to a large number of parameters, but the benefits are manifold, for example, rigorous statistical inference on the common parameter becomes possible and estimates for asset correlations are a by-product. Four industry groups of US firms in 2009 and 2014 are used to demonstrate this estimation method.

Keywords: Distance-to-Default, credit risk model, sequential Monte Carlo method, density-tempered, expanding-data

Author: Jin-Chuan Duan and Christine Wei-Ting Wang

Submitted:
Input files: Y2009G020018.mat, Y2014G020018.mat, airlines in 2009 and 2014, Y2009G020051.mat, Y2014G020051.mat, banks in 2009 and 2014, Y2009G020055.mat, Y2014G020055.mat, insurance in 2009 and 2014, Y2014G020082.mat, Y2009G020082.mat, engineering&construction in 2009 and 2014.
•	Date: 250 daily observations for each firm up to the end of the year
•	FirmCode: classified into 76 industry groups by Bloomberg Industry Classification System (BICS)
•	Equity: market capitalization
•	ShortD: short term debt
•	LongD: long term debt
•	otherL: other liability
•	TotalA: total asset
•	Rf: interest rate
RandFirm_20051.mat: 40 Banks selected in paper. 
RandFirm_20055.mat: 40 Insurances selected in paper.

Output files:
Y2009_Sec20018.mat, Y2014_Sec20018.mat, Y2009_Sec20082.mat, Y2014_Sec20082.mat, Y2009_Sec20051.mat, Y2014_Sec20051.mat, Y2009_Sec20055.mat, Y2014_Sec20055.mat OutY2009_Sec20018.xlsx, OutY2009_Sec20051.xlsx, OutY2009_Sec20055.xlsx, OutY2009_Sec20082.xlsx, OutY2014_Sec20018.xlsx, OutY2014_Sec20051.xlsx, OutY2014_Sec20055.xlsx, OutY2014_Sec20082.xlsx

Example: There are two code files. One is RunSMC_est.m, and the other is TableGenerate.m. In RunSMC_est.m, data from "Input" folder are loaded and all other functions are stored in "subfunciton" folder. The parameters for the initialization sampler are set in funciton SMCsettingbyFirm.m. Then, we sample 1024 particles using the initialization sampler and proceed with the function SMC_Step_MixW.m, which conducts the reweighting and resampling steps. After the resampling step, a function MoveSet_MixW.m runs Metropolis-Hastings (MH) moves. At this stage, the acceptant rate of each block will be displayed in the Command Window. The SMC_Step_MixW.m function will be repeatedly run until it finishes the likelihood density-tempering for the first 5 firms. Then, it will go back to RunSMC_est.m file. 5 firms are added each time, and re-initialization followed by calling SMC_Step_MixW.m will be conducted again and again until reaching N firms. The following files in "Examples" folder are 12 firms from Airline sector in 2014 and 40 random selected firms from Insurance sector in 2009. 
OutY2014_Sec20018.xlsx/ Y2014_Sec20018.mat (Time required to estimate: 0.73 hours) 
OutY2009_Sec20055.xlsx/ Y2009_Sec20055.mat (Time required to estimate: 13.5 hours)
In TableGenerate.m file, Table 1 to Table 4 can be generated from the "Output" data.

Computer system: 
Processor: Intel(R) Xeon(R) CPU W3530 @ 2.80GHz 
Installed memory: 8.00 GB 
Software: MATLAB R2013b
