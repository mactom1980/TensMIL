function [mean_CV_TEST_accuracy, std_val_test_acc, b_PCA_p, quant_p, bins_p, CV_test_accuracy_array, CV_confusMatrix] = ...
    f_bayesParameterLearningFS(X, Y, imID, k, TensMIL_Alg, nr_of_bins, acc_bacc)

% X = instances matrix (n x m)
% Y = instances label matrix (n x m)
% imID = instances ID (1 x m)
% k = nr. of folds for crovalidation


% TensMIL_Alg = 1 means TensMIL (optimizing PCA varince retained recentage
% and number of bins - no intances selection - quant is set to 0)
% nr_of_bins is not used in this frame

% TensMIL_Alg = 2 means TensMIL2 (optimizing PCA variance retained and
% quant - instance selection)
% nr_of_bins is the fixed number of bins used in TENSMIL2

% FOLDS = number of folds for bayes optimization, i.e. the bojective
% function optimizes the 1-accuracy of the mean two fold cross validation
% of the validation set

[X_train_cell, Y_train_cell, patientID_train_cell, X_val_cell, Y_val_cell, patientID_val_cell]  =...
    f_produce_stratified_kFoldValidationSets(X, Y, imID, k);

FOLDS = 2;

% optimizing the hyperparameters

if TensMIL_Alg == 1
    PCA_perc = optimizableVariable('PCA_p',[.05,.6],'Type','real');
    bins_nr = optimizableVariable('bins',[3, 120],'Type','integer');
    bins_p = nan(k,1);
    b_PCA_p = nan(k,1);
    
    for fold = 1:k

        disp(['fold ' num2str(fold)])

        [X_train_cell_val, Y_train_cell_val, patientID_train_cell_val, X_val_cell_val, Y_val_cell_val,...
        patientID_val_cell_val]  = f_produce_stratified_kFoldValidationSets(X_train_cell{fold}, Y_train_cell{fold}, patientID_train_cell{fold}, FOLDS);

        fun = @(x)f_kfold_Accuracy_regr_sub_class_linear_v3_OBJfun...
            (X_train_cell_val, Y_train_cell_val, patientID_train_cell_val, X_val_cell_val, Y_val_cell_val,...
            patientID_val_cell_val, FOLDS, 1, x.PCA_p, x.bins, 0, acc_bacc);

        results = bayesopt(fun,[PCA_perc bins_nr ], 'AcquisitionFunctionName','expected-improvement-plus','Verbose',1,'PlotFcn',[], 'MaxObjectiveEvaluations', 30,'UseParallel',false );%, 'UseParallel',true);
        zbest = bestPoint(results);
       
        b_PCA_p(fold) = zbest.PCA_p;
        bins_p(fold) = zbest.bins;
    end
elseif TensMIL_Alg == 2
    PCA_perc = optimizableVariable('PCA_p',[.05,.8],'Type','real');
    quant_per = optimizableVariable('quant',[.05,.5],'Type','real');
        
    b_PCA_p = nan(k,1);
    quant_p = nan(k,1);
    
    for fold = 1:k

        disp(['fold ' num2str(fold)])

        [X_train_cell_val, Y_train_cell_val, patientID_train_cell_val, X_val_cell_val, Y_val_cell_val,...
        patientID_val_cell_val]  = f_produce_stratified_kFoldValidationSets(X_train_cell{fold}, Y_train_cell{fold}, patientID_train_cell{fold}, FOLDS);

        fun = @(x)f_kfold_Accuracy_regr_sub_class_linear_v3_OBJfun...
            (X_train_cell_val, Y_train_cell_val, patientID_train_cell_val, X_val_cell_val, Y_val_cell_val,...
            patientID_val_cell_val, FOLDS, 1, x.PCA_p, nr_of_bins, x.quant, acc_bacc);

        results = bayesopt(fun,[PCA_perc quant_per], 'AcquisitionFunctionName','expected-improvement-plus','Verbose',1,'PlotFcn',[], 'MaxObjectiveEvaluations', 30,'UseParallel',false );%, 'UseParallel',true);

        zbest = bestPoint(results);
        
        b_PCA_p(fold) = zbest.PCA_p;
        quant_p(fold) = zbest.quant;
       
    end
else
    disp('use TensMIL_Alg = 1 for TensMIL and TensMIL_Alg = 2 for TensMIL2')
    mean_CV_TEST_accuracy = nan; 
    std_val_test_acc= nan;
    b_PCA_p = nan;
    quant_p = nan;
    bins_p = nan;
    CV_test_accuracy_array = nan;
    CV_confusMatrix = nan;
    return
end


% training and testing the model using the learned hyper-parameters
%%
if TensMIL_Alg == 1
    [CV_test_accuracy_array, mean_CV_TEST_accuracy, mean_CV_train_accuracy, T_train, T_test, CV_confusMatrix] = f_kfold_Accuracy_regr_sub_class_linear_v3...
    (X_train_cell, Y_train_cell, patientID_train_cell, X_val_cell, Y_val_cell,...
    patientID_val_cell, k, 1,b_PCA_p, bins_p, zeros(k,1));

    quant_p = 0;
    std_val_test_acc = std(CV_test_accuracy_array);
    disp(['Mean Train Accuracy: ' num2str(mean_CV_train_accuracy)])
    disp(['Mean Test Accuracy: ' num2str(mean_CV_TEST_accuracy) ' (std ' num2str(std_val_test_acc) ')'])
    disp(['Mean Train Time: ' num2str(mean(T_train))])
    disp(['Mean Test Time: ' num2str(mean(T_test))])
else
    [CV_test_accuracy_array, mean_CV_TEST_accuracy, mean_CV_train_accuracy, T_train, T_test, CV_confusMatrix] = f_kfold_Accuracy_regr_sub_class_linear_v3...
    (X_train_cell, Y_train_cell, patientID_train_cell, X_val_cell, Y_val_cell,...
    patientID_val_cell, k, 1,b_PCA_p,  nr_of_bins*ones(k,1), quant_p);

    bins_p = nr_of_bins;
    std_val_test_acc = std(CV_test_accuracy_array);
    disp(['Mean Train Accuracy: ' num2str(mean_CV_train_accuracy)])
    disp(['Mean Test Accuracy: ' num2str(mean_CV_TEST_accuracy) ' (std ' num2str(std_val_test_acc) ')'])
    disp(['Mean Train Time: ' num2str(mean(T_train))])
    disp(['Mean Test Time: ' num2str(mean(T_test))])
end
%%
end