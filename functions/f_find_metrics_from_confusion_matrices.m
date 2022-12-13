function [mean_precision, mean_recall, mean_f1_score] = f_find_metrics_from_confusion_matrices(CV_confusMatrix, nrOfFolds)

precision = cell(3,1);
recall = cell(3,1);
f1_score = cell(3,1);



tmp = size(CV_confusMatrix{1});
nrOfClasses = tmp(1);
for i = 1:3
    precision{i} =  nan(nrOfFolds,1);
    recall{i} = nan(nrOfFolds,1);
    f1_score{i} = nan(nrOfFolds,1);
end




for i=1:nrOfClasses
    for j=1:nrOfFolds
%         if j==4 & i==3
%             disp('huston we got a problem')
%         end
        tmp_mat = CV_confusMatrix{j};
        precision{i}(j) = tmp_mat(i,i)/sum(tmp_mat(:,i));
        recall{i}(j) = tmp_mat(i,i)/sum(tmp_mat(i,:));
        f1_score{i}(j) = 2*(precision{i}(j)*recall{i}(j))/(precision{i}(j) + recall{i}(j));
    end
end
mean_precision = nan(nrOfClasses,1);
mean_recall = nan(nrOfClasses,1);
mean_f1_score = nan(nrOfClasses,1);
for i =1:nrOfClasses
    
    mean_precision(i) = mean(precision{i});
    mean_recall(i) = mean(recall{i});
    mean_f1_score(i) = mean(f1_score{i});
end