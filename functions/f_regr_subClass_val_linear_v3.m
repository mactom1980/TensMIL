function [Val_error, Val_bal_error, conMatrix, pred_score] = f_regr_subClass_val_linear_v3...
    (A_val, Y_Val, tWind_ID_Val, linmodel, LinDiscAnalysis, nrOfBinsPerSub, cuting_points, quant)

[y_pred_by_tWindow_val, y_confidence_val] = predict(linmodel,A_val);

conf_range = y_confidence_val(:,2)-y_confidence_val(:,1);

threshold = min(conf_range(conf_range>=quantile(conf_range, quant)));

y_pred_by_tWindow_val = y_pred_by_tWindow_val(conf_range>= threshold);
Y_Val = Y_Val(conf_range>=threshold);
tWind_ID_Val = tWind_ID_Val(conf_range>=threshold);

patients_val = unique(tWind_ID_Val);
nrOfPatients_val = length(patients_val);
Y_sub_val = zeros(nrOfPatients_val,1);

H_subjs_val = zeros(nrOfPatients_val,nrOfBinsPerSub);

for i=1:nrOfPatients_val
    
    pred_sub_val = y_pred_by_tWindow_val(tWind_ID_Val==patients_val(i));
    Y_sub_val(i) = unique(Y_Val(tWind_ID_Val==patients_val(i)));
%     hist_subj_val = hist(pred_sub_val,cuting_points);
    hist_subj_val = histcounts(pred_sub_val,nrOfBinsPerSub);
    H_subjs_val(i,:) = cumsum(hist_subj_val/sum(hist_subj_val));
end


[y_pred_by_sub_val, pred_score, ~] = predict(LinDiscAnalysis,H_subjs_val);

Val_error = sum(Y_sub_val~=y_pred_by_sub_val)/length(Y_sub_val);
conMatrix = confusionmat(Y_sub_val, y_pred_by_sub_val, 'Order', unique(Y_sub_val));

% balanced error
NrOfClasses = length(unique(Y_Val));
% confMat = confusionmat(Y_sub_val, y_pred_by_sub_val, 'Order', unique(Y_sub_val));
nrOfsubPerClass = nan(NrOfClasses,1);
for i=0:NrOfClasses-1

    nrOfsubPerClass(i+1) = sum(Y_sub_val==i);
    
end

Val_bal_error = 1-sum(diag(conMatrix)./nrOfsubPerClass)/NrOfClasses;
 
end