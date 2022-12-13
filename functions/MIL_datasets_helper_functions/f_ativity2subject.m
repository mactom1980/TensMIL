function [Y_sub, X_sub, patient_ID] = f_ativity2subject(patientIDperActivity,X,Y)

%this function organizes the time windows of the initial 3D tensor X and
%the corresponding label vector Y into two cell arrays X_sub and Y_sub that
%contain the slice of X and Y that correspond to a particular subject.

%In other words f_activity2subject maps the time-windows to subjects.
%
%
%Input: patientID -> an mx1 array containing the labels of each subject for
%each time window
% Y -> an mx1 array conaining the labels of the classes for each time
% window
%X -> an mxnxk tensor containing the data where m -> NrOf time windows,
%n->NrOf bins, k->NrOf Channels (in this framework).

%Actually this function does the maping according the first dimension



patients = unique(patientIDperActivity);
nrOfPatients = length(patients);

Y_sub = cell(nrOfPatients,1);
X_sub = cell(nrOfPatients,1);
patient_ID = cell(nrOfPatients,1);

for i = 1:nrOfPatients
    
    Y_sub{i} = Y(patientIDperActivity==patients(i));
    X_sub{i} = X(patientIDperActivity==patients(i),:,:);
    patient_ID{i} = patients(i);
    
end



end