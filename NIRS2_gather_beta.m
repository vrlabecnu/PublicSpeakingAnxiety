clc; clear;
sub = [7:16 18:19 21:28];
stage = {'pre','sp','rest','bl'};
cond = {'N'};

%% gather beta
for icond = 1:length(cond)
    currcond = cond{icond};
    for isub = 1:length(sub)
        currsub = sub(isub);
        eval(strcat('load(''./nirs-spm/sub',num2str(currsub),'_',currcond,'/SPM_indiv_HbO.mat'');'));
        for istage = 1:length(stage)
            eval([stage{istage} '_beta = SPM_nirs.nirs.beta(' num2str(istage) ',1:44);']);
            eval([currcond '_' stage{istage} '_beta(' num2str(isub) ',1:44)=' stage{istage} '_beta;']);
        end
    end
    dlmwrite(strcat('./nirs-beta/',currcond,'_pre_beta.csv'),eval(strcat(currcond,'_pre_beta')),',');
    dlmwrite(strcat('./nirs-beta/',currcond,'_sp_beta.csv'),eval(strcat(currcond,'_sp_beta')),',');
    dlmwrite(strcat('./nirs-beta/',currcond,'_rest_beta.csv'),eval(strcat(currcond,'_rest_beta')),',');
    dlmwrite(strcat('./nirs-beta/',currcond,'_bl_beta.csv'),eval(strcat(currcond,'_bl_beta')),',');
end
clear pre_beta sp_beta rest_beta bl_beta
clear cond currsub isub SPM_nirs istage

%% plot activation
N_pre_beta = dlmread('./nirs-beta/N_pre_beta.csv',',',1,1);
N_sp_beta = dlmread('./nirs-beta/N_sp_beta.csv',',',1,1);
N_rest_beta = dlmread('./nirs-beta/N_rest_beta.csv',',',1,1);
N_bl_beta = dlmread('./nirs-beta/N_bl_beta.csv',',',1,1);
%%%% bad channels: sub8 - CH18; sub12 - CH18, CH22; sub26 - CH40, CH44.
N_pre_beta(2,18) = NaN; N_pre_beta(6,18) = NaN; N_pre_beta(6,22) = NaN; N_pre_beta(18,40) = NaN; N_pre_beta(18,44) = NaN;
N_sp_beta(2,18) = NaN; N_sp_beta(6,18) = NaN; N_sp_beta(6,22) = NaN; N_sp_beta(18,40) = NaN; N_sp_beta(18,44) = NaN;
N_rest_beta(2,18) = NaN; N_rest_beta(6,18) = NaN; N_rest_beta(6,22) = NaN; N_rest_beta(18,40) = NaN; N_rest_beta(18,44) = NaN;
N_bl_beta(2,18) = NaN; N_bl_beta(6,18) = NaN; N_bl_beta(6,22) = NaN; N_bl_beta(18,40) = NaN; N_bl_beta(18,44) = NaN;

currcond = 'N';
for istage = 1:length(stage)
    eval(['avg_' currcond '_' stage{istage} ' = nanmean(' currcond '_' stage{istage} '_beta);']);
    eval(['plotTopoMap(avg_' currcond '_' stage{istage} '(1:22),''3x5'',[-0.1,0.2]);']);
    colormap hot;
    saveas(gcf, ['./nirs-activation/' currcond '_' stage{istage} '_RH'], 'fig');
    eval(['plotTopoMap(avg_' currcond '_' stage{istage} '(23:44),''3x5'',[-0.1,0.2]);']);
    colormap hot;
    saveas(gcf, ['./nirs-activation/' currcond '_' stage{istage} '_LH'], 'fig');
end
% clear -regexp _beta$
clear istage

%% plot contrast
for iCH = 1:44
    [H,P,CI,STATS] = ttest(N_pre_beta(:,iCH),N_rest_beta(:,iCH),0.05,'both');
    contrast1_N(iCH) = STATS.tstat;
    contrast1_N_p(iCH) = P;
    [H,P,CI,STATS] = ttest(N_sp_beta(:,iCH),N_bl_beta(:,iCH),0.05,'both');
    contrast2_N(iCH) = STATS.tstat;
    contrast2_N_p(iCH) = P;
end
plotTopoMap(contrast1_N(1:22),'3x5',[-1,5]); colormap hot; saveas(gcf, './nirs-activation/contrast1_N_RH', 'fig');
plotTopoMap(contrast1_N(23:44),'3x5',[-1,5]); colormap hot; saveas(gcf, './nirs-activation/contrast1_N_LH', 'fig');
plotTopoMap(contrast2_N(1:22),'3x5',[-1,5]); colormap hot; saveas(gcf, './nirs-activation/contrast2_N_RH', 'fig');
plotTopoMap(contrast2_N(23:44),'3x5',[-1,5]); colormap hot; saveas(gcf, './nirs-activation/contrast2_N_LH', 'fig');

% clear -regexp _beta$
clear istage

%% transactivation
currcond = 'N';
N_dev1 = N_sp_beta - N_pre_beta;
N_dev2 = N_bl_beta - N_rest_beta;
avg_N_dev1 = mean(N_dev1);
avg_N_dev2 = mean(N_dev2);

plotTopoMap(avg_N_dev1(:,1:22),'3x5',[-0.1,0.2]);
colormap hot;
saveas(gcf,'./nirs-activation/N_dev1_RH', 'fig');
plotTopoMap(avg_N_dev1(:,23:44),'3x5',[-0.1,0.2]);
colormap hot;
saveas(gcf,'./nirs-activation/N_dev1_LH', 'fig');

plotTopoMap(avg_N_dev2(:,1:22),'3x5',[-0.1,0.2]);
colormap hot;
saveas(gcf,'./nirs-activation/N_dev2_RH', 'fig');
plotTopoMap(avg_N_dev2(:,23:44),'3x5',[-0.1,0.2]);
colormap hot;
saveas(gcf,'./nirs-activation/N_dev2_LH', 'fig');