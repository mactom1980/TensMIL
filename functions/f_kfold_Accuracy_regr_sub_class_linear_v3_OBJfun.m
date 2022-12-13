function [mean_CV_val_error] = f_kfold_Accuracy_regr_sub_class_linear_v3_OBJfun...
    (X_train_cell, Y_train_cell, patientID_train_cell, X_val_cell, Y_val_cell,...
    patientID_val_cell, k, SVD_on,SVD_per, nrOfBinsPerSub, quant, acc_bacc)

if SVD_on==1
    SVD_disp = 'SVD on';
else
    SVD_disp = 'SVD off';
end

rng(0,'twister')

CV_train_error = nan(k,1);
CV_val_error = nan(k,1);
CV_bal_error = nan(k,1);
conMatrix = cell(k,1);
for fold = 1:k

    if SVD_on == 1
        A_train_val = [X_train_cell{fold};X_val_cell{fold}];
        [U, S, ~] = svd(A_train_val, "econ"); % use the econ
        sigValues = diag(S);
        cSum = cumsum(sigValues)./sum(sigValues);
        ind_of_least_informativeSigVals = find(cSum>=SVD_per);
        A_train_val = U*S(:,1:ind_of_least_informativeSigVals(1));
        A_train = A_train_val(1:size(X_train_cell{fold},1),:);
        A_val = A_train_val(size(X_train_cell{fold},1)+1:end,:);
    else
        A_train = X_train_cell{fold};
        A_val = X_val_cell{fold};
    end

    try
        [CV_train_error(fold), linmodel, LinDiscAnalysis, cuting_points] =...
        f_regr_subClass_train_linear_v3(A_train, Y_train_cell{fold}, patientID_train_cell{fold}, nrOfBinsPerSub,quant);
        [CV_val_error(fold), CV_bal_error(fold), conMatrix{fold}] = ...
        f_regr_subClass_val_linear_v3(A_val, Y_val_cell{fold}, patientID_val_cell{fold}, linmodel, LinDiscAnalysis, nrOfBinsPerSub, cuting_points,quant);
    catch
        disp('Warning! Something went worng with the classifier. Test and Validation Error is set manually to one (1).')
        CV_train_error(fold) = 1;
        CV_val_error(fold) = 1;
        CV_bal_error(fold) = 1;
    end
end

if acc_bacc == 1 % we optimize the accuracy
    mean_CV_val_error = mean(CV_val_error);
    
elseif acc_bacc == 2 % we optimize balanced accuracy
    mean_CV_val_error = mean(CV_bal_error);
else
    diplay(['use acc_bacc = 1 for accurcy and acc_bacc = 2 for balanced accuracy'])
    return
end
