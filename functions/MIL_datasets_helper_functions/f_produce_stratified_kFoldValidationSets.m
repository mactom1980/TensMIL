function [X_train_cell, Y_train_cell, patientID_train_cell, X_val_cell, Y_val_cell,...
    patientID_val_cell]  = f_produce_stratified_kFoldValidationSets(X, Y, patientID, k)

rng(0,'twister')

X_train_cell = cell(k,1);
Y_train_cell = cell(k,1);
patientID_train_cell = cell(k,1);
X_val_cell = cell(k,1);
Y_val_cell = cell(k,1);
patientID_val_cell = cell(k,1);

[Y_sub, X_sub, patientID_sub] = f_ativity2subject(patientID,X,Y);
nrOfPatients = length(Y_sub);

Y_sub_array = -ones(nrOfPatients,1);

for i=1:nrOfPatients
        if length(unique(Y_sub{i}))==1
            Y_sub_array(i) = unique(Y_sub{i});
        else
            disp(['A bag contains not homogeneous lebels: ' num2str(patientID_sub{i} )])
        end
end

c = cvpartition(Y_sub_array,'KFold',k);

for i=1:k
    X_train_sub = [];
    Y_train_sub = [];
    patient_ID_tainsub = [];
    
    X_val_sub = [];
    Y_val_sub = [];
    patientID_val_sub = [];
    
    X_train_sub = X_sub(c.training(i));
    Y_train_sub = Y_sub(c.training(i));
    patient_ID_tainsub = patientID_sub(c.training(i));
    
    [X_train_cell{i},Y_train_cell{i}, patientID_train_cell{i}] = f_subject2activity(X_train_sub,Y_train_sub,patient_ID_tainsub);

    X_val_sub = X_sub(c.test(i));
    Y_val_sub = Y_sub(c.test(i));
    patientID_val_sub = patientID_sub(c.test(i));
    
    [X_val_cell{i},Y_val_cell{i}, patientID_val_cell{i}] = f_subject2activity(X_val_sub,Y_val_sub,patientID_val_sub);
end

end