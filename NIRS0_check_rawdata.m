clc; clear;
sub = [7:16 18:19 21:28];
cond = 'N';
for isub = 1:length(sub)
    currsub = sub(isub);
    disp(['Now checking subject ' num2str(currsub)]);
    eval(strcat('load(''./nirs-transfer/sub',num2str(currsub),'_',cond,'.mat'');'));
    means(isub,:) = mean(nirs_data.oxyData);
    stds(isub,:) = std(nirs_data.oxyData);
    maxs(isub,:) = max(nirs_data.oxyData);
    mins(isub,:) = min(nirs_data.oxyData);
    badstd = find(stds(isub,:)>0.5);
    badscale = union(find(maxs(isub,:)>5),find(mins(isub,:)<-5));
    badCH = intersect(badstd,badscale);
    if isempty(badCH)
        disp('No bad channel.');
    else
        disp('Bad channels:');
        disp(badCH);
    end
    disp(' ');
end

clear currsub isub nirs_data

%% std check
clc; clear;
sub = [7:16 18:19 21:28];
cond = 'N';
for isub = 1:length(sub)
    currsub = sub(isub);
    disp(['Now checking subject ' num2str(currsub)]);
    eval(strcat('load(''./nirs-transfer/sub',num2str(currsub),'_',cond,'.mat'');'));
    means(isub,:) = mean(nirs_data.oxyData);
    stds(isub,:) = std(nirs_data.oxyData);
    maxs(isub,:) = max(nirs_data.oxyData);
    mins(isub,:) = min(nirs_data.oxyData);
    badCH = union(find(means(isub,:)>(mean(means(isub,:))+mean(stds(isub,:))*3)),find(means(isub,:)<(mean(means(isub,:))-mean(stds(isub,:))*3)));
    if isempty(badCH)
        disp('No bad channel.');
    else
        disp('Bad channels:');
        disp(badCH);
    end
    disp(' ');
end

clear currsub isub nirs_data