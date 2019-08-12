% fNIRS preprocessing and GLM
clc; clear;

study_path = '.';
subjectnum = [7:16 18:19 21:28];
condition = {'N'};

for isub = subjectnum
    for condnum = 1:length(condition)
        icond = condition{condnum};
        disp(['Now processing sub' num2str(isub) '_' icond '...']);
        %         %% step 1: data conversion (Hitachi ETG-7100)
        %         fname_nirs = '..\rawdata\*.csv';
        %         [nirs_data] = data_conversion_batch(fname_nirs, system, fs, dist, wavelength, DPF, flag_DPF_correction, [], [], []);
        %         save('C:\NIRS_SPM_v4\Sample_data\NIRS_data_file\converted_NIRS.mat', 'nirs_data');
        
        %% step 2: model specification
        fname_nirs = [study_path filesep 'nirs-transfer\sub' num2str(isub) '_' icond '.mat'];
        hb = 'hbo';
        HPF = 'wavelet';
        LPF = 'hrf';
        method_cor = 0;
        dir_save = [study_path filesep 'nirs-spm\sub' num2str(isub) '_' icond];
        mkdir(dir_save);
        flag_window = 0;
        hrf_type = 0;
        units = 1;
        names{1} = 'prepare';
        onsets{1} = 20;
        durations{1} = 20;
        names{2} = 'speaking';
        onsets{2} = 60;
        durations{2} = 180;
        names{3} = 'rest';
        onsets{3} = 280;
        durations{3} = 20;
        names{4} = 'baseline';
        onsets{4} = 320;
        durations{4} = 20;
        [SPM_nirs] = specification_batch(fname_nirs, hb, HPF, LPF, method_cor, dir_save, flag_window, hrf_type, units,  names, onsets, durations);
        
        %% step 3: model estimation
        fname_nirs = [study_path filesep 'nirs-transfer\sub' num2str(isub) '_' icond '.mat'];
        fname_SPM = [study_path filesep 'nirs-spm\sub' num2str(isub) '_' icond filesep 'SPM_indiv_HbO.mat'];
        [SPM_nirs] = estimation_batch(fname_SPM, fname_nirs);
        
        %% step 4: activation map (T-statistic, EC-corrected p-value < 0.05)
        fname_SPM = ['.\nirs-spm\sub' num2str(isub) '_' icond filesep 'SPM_indiv_HbO.mat'];
        fname_ch = '.\trans-registration\registration.mat';
        
        con_name = 'prepare';
        con_vec = [1 0 0 0 0]';
        STAT = 'T';
        spec_hemi = 'dorsal';
        p_value = 0.05;
        correct_p = 'EC';
        disp_fig = 0;
        [stat_brain, act_brain, threshold] = activation_map_batch(fname_SPM, fname_ch, con_name, con_vec, STAT, spec_hemi, p_value, correct_p, disp_fig);
        
        con_name = 'speaking';
        con_vec = [0 1 0 0 0]';
        STAT = 'T';
        spec_hemi = 'dorsal';
        p_value = 0.05;
        correct_p = 'EC';
        disp_fig = 0;
        [stat_brain, act_brain, threshold] = activation_map_batch(fname_SPM, fname_ch, con_name, con_vec, STAT, spec_hemi, p_value, correct_p, disp_fig);
        
        con_name = 'rest';
        con_vec = [0 0 1 0 0]';
        STAT = 'T';
        spec_hemi = 'dorsal';
        p_value = 0.05;
        correct_p = 'EC';
        disp_fig = 0;
        [stat_brain, act_brain, threshold] = activation_map_batch(fname_SPM, fname_ch, con_name, con_vec, STAT, spec_hemi, p_value, correct_p, disp_fig);
        
        con_name = 'baseline';
        con_vec = [0 0 0 1 0]';
        STAT = 'T';
        spec_hemi = 'dorsal';
        p_value = 0.05;
        correct_p = 'EC';
        disp_fig = 0;
        [stat_brain, act_brain, threshold] = activation_map_batch(fname_SPM, fname_ch, con_name, con_vec, STAT, spec_hemi, p_value, correct_p, disp_fig);
    end
end