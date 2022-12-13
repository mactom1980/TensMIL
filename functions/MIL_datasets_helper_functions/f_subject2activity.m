function [X,Y, tWindow_ID] = f_subject2activity(Xsub_cell,Ysub_cell,sub_ID)

%Xsub_cell and Ysub_cell are the cell containing the activities of each subject
%Xsub_cell and Ysub_cell musts have the same number of rows
%X,Y are a matrix and a vector that contain the activities

X = [];
Y = [];
tWindow_ID = [];
for i = 1:length(Ysub_cell)
    
    X = [X;Xsub_cell{i}];
    Y = [Y;Ysub_cell{i}];
    
    tWindow_ID = [tWindow_ID;repmat(sub_ID{i},length(Ysub_cell{i}),1)];
end

end