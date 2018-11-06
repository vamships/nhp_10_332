# ******************************************************************************
# Author      : Srivamshi Pittala (srivamshi.pittala.gr@dartmouth.edu)
# Advisor     : Prof. Chris Bailey-Kellogg (cbk@cs.dartmouth.edu) & Prof. Margaret E. Ackerman
# Project     : NIH 10-332
# Description : Help file for running the scripts to analyze luminex and functional measurements
#				associated with the HIV vaccine trial 10-332
# Cite        : TBD
# ******************************************************************************

# Copyright (C) <2018>  <Srivamshi Pittala>

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

#***********************
Input Files (in the directory in/)
1. luminex_tp5.csv	: biophysical measurements
2. functions.csv 	: functional measurements
3. subjects.csv		: group and survival information
#***********************

#***********************
For survival analysis, use the following three scripts
#***********************

#-----------------------
1. coxPredRisk_biophysical.R
Takes about 10 min to complete
#-----------------------

===> This script performs the survival analysis using the biophysical measurements (luminex). 
		(1)	First, the features are clustered using hierarchical clustering.
		(2) Then repeated cross-validation is performed on the features returned from clustering. In this 	cross-validation, greedy backward elimination is used to select features for each fold using the training set. The most-frequent features appearing in the repeated cross-validation are used to build a final model, the results from which are used for the figures.
		(3) The features in the final model are used to discover other features that could be equally predictive of risk but were not considered since they are correlated with the final feature set.
		(4) Permutation tests are performed by shuffling the rows of the outcome labels, but keeping the rows of feature matrix the same. The permuted data are sent through the same pipeline as were the actual data (i.e. pre-filtering, repeated cross-validation, and final model evaluation). This is repeated multiple times, independently permuting the outcome labels every time.
		(5) Robustness is estimated by comparing the C-indices from using actual features to those of using permuted features.

#-----------------------
2. coxPredRisk_randomSel.R
Takes about 5 min to complete
#-----------------------

===> This script performs survival analysis using the biophysical measurements, but with robustness tests by size-matched random feature selection
		(1) Robustness tests are performed by randomly selecting features from the actual feature set. The features are size-matched with the final feature set. The randomly selected features are sent through the same pipeline as were the final features (i.e. final model evaluation). This is repeated multiple times, independently selecting random features every time.
		(2) Robustness is estimated by comparing the C-indices from using actual features to those of using random features.

#-----------------------
3. coxPredRisk_functional.R
Takes about 5 min to complete
#-----------------------
===> This script performs survival analysis using the two functional measurements ADCP and ADNP
		(1) Repeated cross-validation is performed using the two features.
		(2) Permutation tests are performed by shuffling the rows of the outcome labels, but keeping the rows of feature matrix the same. The permuted data are sent through the same pipeline as were the actual data. This is repeated multiple times, independently permuting the outcome labels every time.
		(3) Robustness is estimated by comparing the C-indices from using actual features to those of using permuted features.

#***********************
For multinomial logistic classification, use the following script
#***********************

#-----------------------
1. predClassM_biophysical.R
Takes about 5 min to complete
#-----------------------
===> This script performs multinomial logistic classification to identify vaccine groups using the biophysical measurements
		(1) Classification is done using the lasso-regularized multinomial logistic regression on the features. The best regularization parameter (lambda) is chosen to be the one with lowest classification error. This is repeated multiple times. A final model is trained and evaluated by using a fixed seed to determine folds.
		(2) Permutation tests are performed by shuffling the rows of the class labels, but keeping the rows of feature matrix the same. The permuted data are sent through the same pipeline as were the actual data. This is repeated multiple times, independently permuting the class labels every time.
		(3) Robustness tests are also performed by randomly selecting features from the actual feature set. The features are size-matched with the final feature set. The randomly selected features are sent through the same pipeline as were the actual features. This is repeated multiple times, independently selecting random features every time.
		(4) Robustness is estimated by comparing the accuracies from using actual features to those of using permuted and random features.

#***********************
For generating the figures as in the manuscript, use the following two scripts
#***********************

#-----------------------
1. generate_networkPlot.R
Takes about 2 min to complete
#-----------------------
===> This script generates the network plot to visualize the co-correlates of the biophysical measurements identified by the final model. Run this after 'coxPredRisk_biophysical.R' has successfully finished without errors.
	(1) The nodes of the network are the features. Edges are used to representation correlations that are above 0.75
	(2) The output of the file is an xml file, that can be used by a network visualization software like Cytoscape to finally generate the figure.
	(3) I used Cytoscape software (version 3.5.1). The steps follow:
		(i)		Open the Cytoscape software application
		(ii)	File->Import->Network->File...
		(iii)	Choose the nodes.xml file
		(iv)	Layout->Attribute Circle Layout->group
		(v)		File->Export as Image...
		(vi)	Save file as Features.png in the same directory as nodes.xml

#-----------------------
2. generate_figures.R
#-----------------------
===> This script generates the figures as shown in the manuscript. Run this after after all the above scripts have finished successfully (i.e. coxPredRisk_biophysical.R, coxPredRisk_randomSel.R, coxPredRisk_functional.R, predClassM_biophysical.R, and generate_networkPlot.R). The directory ‘results_figures_reference/’ can be used as a reference to what the final figure outputs should look like.

# ******************************************************************************

#-----------------------
System configuration 
#-----------------------
OS 	: Ubuntu 14.04.5 LTS
CPU : i7 8 cores @ 3.6 GHz
RAM : 24 GB

#-----------------------
Software and packages (version)
#-----------------------
01. R (3.4.3) and RStudio (1.0.143)
02. glmnet (2.0-13)
03. gplots (3.0.1)
04. ggplot2 (2.2.1)
05. RColorBrewer (1.1-2)
06. effsize (0.7.1)
07. survival (2.41-3)
08. corrplot (0.84)
09. caret (6.0-77)	(Install with dependencies=TRUE)
10. survcomp (1.24.0)	(Installed via bioconductor)
11. polycor (0.7-9)
12. elasticnet (1.1)
13. e1071(1.6-8)

#-----------------------
Functions used by the scripts (in the directory funcs/)
#-----------------------
01. clusterFeats.R
02. convertKMtoEvents.R
03. coxSurvivalFinalModel.R
04. coxSurvivalLookBack.R
05. coxSurvivalWithBackSearch.R
06. createColumnColors.R
07. createSubjectColors.R
08. doCoxBackSearch.R
09. doFullSurvivalAnalysis.R
10. extractProbabilityFromKM.R
11. glmnetMultiClass.R
12. heatmap4.R
13. plotConfusion.R
14. selectFromClusters.R
15. takeOffOneFeat.R
