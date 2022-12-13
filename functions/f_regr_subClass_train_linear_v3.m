function [Tr_error, linmodel, LinDiscAnalysis, cuting_points] = ...
    f_regr_subClass_train_linear_v3(A_train, Y_Train, tWind_ID_Train, nrOfBinsPerSub, quant)

robust_vals.RobustWgtFun = 'logistic';

f_turn_robustfittingWarning_off(true)

linmodel = fitlm(A_train,Y_Train,'quadratic','RobustOpts',robust_vals);

[y_pred_by_tWindow_train, y_confidence_range_train] = predict(linmodel,A_train);

conf_range = y_confidence_range_train(:,2)- y_confidence_range_train(:,1);
threshold = min(conf_range(conf_range>=quantile(conf_range, quant)));

y_pred_by_tWindow_train = y_pred_by_tWindow_train(conf_range>=threshold);
tWind_ID_Train = tWind_ID_Train(conf_range>=threshold);
Y_Train = Y_Train(conf_range>=threshold);

patients_train = unique(tWind_ID_Train);
nrOfPatients_train = length(patients_train);

H_subjs = zeros(nrOfPatients_train,nrOfBinsPerSub);
Y_sub = zeros(nrOfPatients_train,1);

%find quantiles
cuting_points = zeros(nrOfBinsPerSub+1, 1);
cuting_points(2:end) = quantile(y_pred_by_tWindow_train,nrOfBinsPerSub);

% cuting_points = quantile(y_pred_by_tWindow_train,nrOfBinsPerSub);

for i=1:nrOfPatients_train
    
    pred_sub = y_pred_by_tWindow_train(tWind_ID_Train==patients_train(i));
    Y_sub(i) = unique(Y_Train(tWind_ID_Train==patients_train(i)));
%     hist_subj = hist(pred_sub,cuting_points);
    hist_subj = histcounts(pred_sub,nrOfBinsPerSub);
    
    H_subjs(i,:) = cumsum(hist_subj/sum(hist_subj));
end

LinDiscAnalysis = fitcdiscr(H_subjs,Y_sub,'DiscrimType','pseudoquadratic');
y_pred_by_sub = predict(LinDiscAnalysis,H_subjs);

Tr_error = sum(Y_sub~=y_pred_by_sub)/length(Y_sub);

end
