README - Figure3_Machine_learning
#################################

*************************
*** Package versions ***
caret v.6.0-80
randomForest v.4.6-14
gbm v.2.1.3
kernlab v.0.9-27
*************************
***   Other software   ***
DiNUP v1.3
DANPOS v2.2
homer	v4.10.3
*************************

* Feature_extraction:
	- Get_features_FC_v7.R: script to extract features for Machine Learning
	- Features_extraction_functions_FC_v6.R: functions called in above script
	
	The output of feature extraction is provided in:
		"Figure3_Machine_learning/Feature_extraction/Output_RData/All_ML_input_matrices_up_to_date.RData"
		
 * Build_Models:
 	- Both
 		* Learn_aging_change_NEW_oldie_SVM_vNewCaret.R: run RF machine learning across tissues
 		* SVM_helper_functions_v2.R: functions called in above script
 		
 		* Learn_aging_change_NEW_oldie_NNET_vNewCaret.R: run RF machine learning across tissues
 		* NNET_helper_functions_v1.R: functions called in above script
 		
 		* Learn_aging_change_NEW_oldie_RF_vNewCaret.R: run RF machine learning across tissues
 		* RF_helper_functions_v3.R: functions called in above script
 		
 		* Learn_aging_change_NEW_oldie_GBM_vNewCaret.R: run RF machine learning across tissues
 		* GBM_helper_functions_v3.R: functions called in above script
 		
 		* RF_GBM_summary_feature_importance: 
 			summary_importance_GBM_vnewCaret.R: parse feature importance in GBM models
			summary_importance_RF_vnewCaret.R: parse feature importance in RF models

	- Dynamic (helper functions identical to 'Both' case)
		* Learn_aging_change_NEW_oldie_GBM_vDYNAMIC.R
		* Learn_aging_change_NEW_oldie_NNET_vDYNAMIC.R
		* Learn_aging_change_NEW_oldie_RF_vDYNAMIC.R
		* Learn_aging_change_NEW_oldie_SVM_vDYNAMIC.R
		
	- Static (helper functions identical to 'Both' case)
		* Learn_aging_change_NEW_oldie_GBM_vSTATIC.R
		* Learn_aging_change_NEW_oldie_NNET_vSTATIC.R
		* Learn_aging_change_NEW_oldie_RF_vSTATIC.R
		* Learn_aging_change_NEW_oldie_SVM_vSTATIC.R
		
		All_ML_input_matrices_up_to_date.RData: RData file with all extracted features used in Machine Learning

* Accuracy_parsing:
	- GBM: (Both, Dynamic and Static)
		Always a parsing code paired with a function code
		
	- RF: (Both, Dynamic and Static)
		Always a parsing code paired with a function code

* Accuracies_comparison;
	- Compare_accuracies.R: plot heatmap of accuracies from 'static', 'dynamic' and 'both' models.
	
* Nucleosomes: used in Machine learning
	- dinup_scg4.sh: DINUP calling
	- nucleosome_calling_batch_local_v2.sh: DANPOS calling
	- consensus, DANPOS and DiNUP folders contain the final differential nucleosome calls
