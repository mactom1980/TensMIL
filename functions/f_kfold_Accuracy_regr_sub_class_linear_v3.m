function [CV_val_accuracy_array, mean_CV_val_accuracy, mean_CV_train_accuracy, T_train, T_test, CV_confusMatrix, pred_score_per_fold] = f_kfold_Accuracy_regr_sub_class_linear_v3...
    (X_train_cell, Y_train_cell, patientID_train_cell, X_val_cell, Y_val_cell,...
    patientID_val_cell, k, SVD_on,SVD_per, nrOfBinsPerSub, quant)

rng(0,'twister')

CV_train_error = nan(k,1);
CV_val_error = nan(k,1);
CV_confusMatrix = cell(k,1);
pred_score_per_fold = cell(k,1);

T_train = zeros(k,1);
T_test = zeros(k,1);
for fold = 1:k
    tic;
    if SVD_on == 1
        A_train_val = [X_train_cell{fold};X_val_cell{fold}];
        [U, S, ~] = svd(A_train_val, "econ");
        sigValues = diag(S);
        cSum = cumsum(sigValues)./sum(sigValues);
        ind_of_least_informativeSigVals = find(cSum>=SVD_per(fold));
        A_train_val = U*S(:,1:ind_of_least_informativeSigVals(1));
        A_train = A_train_val(1:size(X_train_cell{fold},1),:);
        A_val = A_train_val(size(X_train_cell{fold},1)+1:end,:);
    else
        A_train = X_train_cell{fold};
        A_val = X_val_cell{fold};
    end

    [CV_train_error(fold), linmodel, LinDiscAnalysis, cuting_points] =...
    f_regr_subClass_train_linear_v3(A_train, Y_train_cell{fold}, patientID_train_cell{fold}, nrOfBinsPerSub(fold),quant(fold));
    
    T_train(fold) = toc;

    tic;
    
    [CV_val_error(fold), ~, CV_confusMatrix{fold}, pred_score_per_fold{fold}] = ...
    f_regr_subClass_val_linear_v3(A_val, Y_val_cell{fold}, patientID_val_cell{fold}, linmodel, LinDiscAnalysis, nrOfBinsPerSub(fold), cuting_points,quant(fold));
    
    T_test(fold) = toc;

    disp(['fold ' num2str(fold) ': Tr. Ac. = ' num2str(1-CV_train_error(fold)) ', Test acc. = ' num2str(1-CV_val_error(fold))])
end

mean_CV_train_accuracy = 1 - mean(CV_train_error);
mean_CV_val_accuracy = 1 - mean(CV_val_error);
CV_val_accuracy_array = 1-CV_val_error;

end
