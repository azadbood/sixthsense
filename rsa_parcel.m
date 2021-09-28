% pattern similarity analysis in cortical brain parcels (Schaefer et al 2018)
% written by Asieh Zadbood
function rsa_parcel(rr)
rand('state', sum(100*clock));
rng(sum(100*clock)+rr,'twister');

%% setting up
cond1=1; % spoiled
cond2=2; % twist
cond3=3; % no-twist

% which events
impscenes = [1 4 5 6 9 12 16];

% setting paths
basepath = 'Y:\Asieh\MovieTwistProject\'; addpath(basepath);
addpath('Y:\Asieh\isc_nifti_kit\NIFTI_tools');addpath('Y:\Asieh\isc_nifti_kit\');
addpath('Y:\Asieh\MovieTwistProject\funcs\');

% loading rois
roi_fnames = 'Schaefer_100_icm152';
roi_mask = load_nii(fullfile(basepath,'masks','mm',[roi_fnames '.nii']));
mask = reshape(roi_mask.img,[],1);
fprintf([' rois read from : ' roi_fnames '\n']);

% read event info
load(fullfile(basepath,'segments','final_recall_srm_cleanmotion.mat'));
load(fullfile(basepath,'segments','final_recall_cleanmotion.mat'));
load(fullfile(basepath,'segments','movie_events.mat'));
fid = fopen(fullfile(basepath,'bids','testnames_cleanmotion.txt'));
data = textscan(fid,'%s%s%s%s','HeaderLines',0,'CollectOutput',1);
data = data{:};
fclose(fid);
mevents = movie_events;
mevents(end,end)=2266; % cut movie end before the final scene for all groups

% create name lists for each condition
names1 = data(str2num(cell2mat(data(:,3)))==cond1,2);
names2 = data(str2num(cell2mat(data(:,3)))==cond2,2);
names3 = data(str2num(cell2mat(data(:,3)))==cond3,2);

%% read TR by TR and avg movie data, avg all for no-twist to get the G
for subj=1:length(names1) % spoiled condition
    load(fullfile(basepath,'data',roi_fnames,names1{subj},[names1{subj} '_' num2str(cond1) '_movie_' num2str(rr) '.mat'])); %loading roi data
    
    for e=1:length(impscenes)
        avgdata(:,e) = mean(subj_tcc(:,mevents(impscenes(e),1):mevents(impscenes(e),2)),2); % selecting events and average timecourse within events
    end
    mgroup1_tr(:,:,subj) = subj_tcc; % keep the entire timecourse data
    mgroup1_avg(:,:,subj) = avgdata; % keep the scene averaged data
end

for subj=1:length(names2) % twist condition
    load(fullfile(basepath,'data',roi_fnames,names2{subj},[names2{subj} '_' num2str(cond2) '_movie_' num2str(rr) '.mat'])); %loading roi data
    
    for e=1:length(impscenes)
        avgdata(:,e) = mean(subj_tcc(:,mevents(impscenes(e),1):mevents(impscenes(e),2)),2); % selecting events and average timecourse within events
    end
    mgroup2_tr(:,:,subj) = subj_tcc; % keep the entire timecourse data
    mgroup2_avg(:,:,subj) = avgdata; % keep the scene averaged data
end

for subj=1:length(names3) % no-twist condition
    load(fullfile(basepath,'data',roi_fnames,names3{subj},[names3{subj} '_' num2str(cond3) '_movie_' num2str(rr) '.mat'])); %loading roi data
    
    for e=1:length(impscenes)
        avgdata(:,e) = mean(subj_tcc(:,mevents(impscenes(e),1):mevents(impscenes(e),2)),2); % selecting events and average timecourse within events
    end
    mgroup3_tr(:,:,subj) = subj_tcc; % keep the entire timecourse data
    mgroup3_avg(:,:,subj) = avgdata; % keep the scene averaged data
    clear subj_tcc
end
fprintf(' movie data read \n');

%% read and average (recall) twist and no-twist movie
for subj=1:length(names1) % spoiled condition
    corr_scn = correct{1,cond1}(:,subj);
    
    subject_b = load_nii(fullfile(basepath,'bids','output','LSS',names1{subj},[names1{subj} '_recall_LSS.nii']));
    % add nan for single events in cut out scans
    if strcmp(names1{subj},'sub-08') || strcmp(names1{subj},'sub-35') || strcmp(names1{subj},'sub-44')
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),17);
        betas = zscore(betas,0,2);
        betas = cat(2,nan(size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),1),betas);
    elseif strcmp(names1{subj},'sub-59')
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),17);
        betas = zscore(betas,0,2);
        betas = cat(2,betas,nan(size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),1));
    else
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),18);
        betas = zscore(betas,0,2);
    end
    subj_r = betas(find(mask==rr),impscenes);
    rgroup1_r(:,:,subj) = subj_r; % recall data for all selected scenes
    corr_scn = corr_scn(impscenes);
    subj_r(:,find(corr_scn==0))=nan;
    rgroup1(:,:,subj) = subj_r; % recall data for all selected scenes with CORRECT recall
    
    clear subject_b subj_r betas
end

for subj=1:length(names2) % twist condition
    corr_scn = correct{1,cond2}(:,subj);
    subject_b = load_nii(fullfile(basepath,'bids','output','LSS',names2{subj},[names2{subj} '_recall_LSS.nii']));
    % add nan for single events in cut out scans
    if strcmp(names2{subj},'sub-08') || strcmp(names2{subj},'sub-35') || strcmp(names2{subj},'sub-44')
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),17);
        betas = zscore(betas,0,2);
        betas = cat(2,nan(size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),1),betas);
    elseif strcmp(names2{subj},'sub-59')
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),17);
        betas = zscore(betas,0,2);
        betas = cat(2,betas,nan(size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),1));
    else
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),18);
        betas = zscore(betas,0,2);
    end
    subj_r = betas(find(mask==rr),impscenes);
    rgroup2_r(:,:,subj) = subj_r; % recall data for all selected scenes
    corr_scn = corr_scn(impscenes);
    subj_r(:,find(corr_scn==0))=nan;
    rgroup2(:,:,subj) = subj_r;  % recall data for all selected scenes with CORRECT recall
    
    clear subject_b subj_r betas
end

for subj=1:length(names3) % no-twist condition
    corr_scn = correct{1,cond3}(:,subj);
    subject_b = load_nii(fullfile(basepath,'bids','output','LSS',names3{subj},[names3{subj} '_recall_LSS.nii']));
    % add nan for single events in cut out scans
    if strcmp(names3{subj},'sub-08') || strcmp(names3{subj},'sub-35') || strcmp(names3{subj},'sub-44')
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),17);
        betas = zscore(betas,0,2);
        betas = cat(2,nan(size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),1),betas);
    elseif strcmp(names3{subj},'sub-59')
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),17);
        betas = zscore(betas,0,2);
        betas = cat(2,betas,nan(size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),1));
    else
        betas = reshape(subject_b.img,size(subject_b.img,1)*size(subject_b.img,2)*size(subject_b.img,3),18);
        betas = zscore(betas,0,2);
    end
    subj_r = betas(find(mask==rr),impscenes);
    rgroup3_r(:,:,subj) = subj_r; % recall data for all selected scenes
    corr_scn = corr_scn(impscenes);
    subj_r(:,find(corr_scn==0))=nan;
    rgroup3(:,:,subj) = subj_r; % recall data for all selected scenes with CORRECT recall
    
    clear subject_b subj_r betas
end

fprintf(' recall data read \n');

%% correlations

% %%%%%%%%%%%%%% MOVIE %%%%%%%%%%%%%%%%%%
[all.mm.diag_sub_m2m3,all.mm.pval_m2m3] = corr_one2avgofothers(mgroup2_avg,mgroup3_avg,1,0,'m2m3'); % avg of all others
[all.mm.diag_sub_m2m1,all.mm.pval_m2m1] = corr_one2avgofothers(mgroup2_avg,mgroup1_avg,1,0,'m2m1'); % avg of all others
[all.mm.realdiff, all.mm.pvaldiff] = non_param_t(all.mm.diag_sub_m2m3,all.mm.diag_sub_m2m1,1); % based on vec1-vec2, type=1 paired,type2=unpaired
[all.mm.realdiff_z, all.mm.pvaldiff_z] = non_param_t(atanh(all.mm.diag_sub_m2m3),atanh(all.mm.diag_sub_m2m1),1); % based on vec1-vec2, type=1 paired,type2=unpaired

% %%%%%%%%%%%%%% RECALL %%%%%%%%%%%%%%%%%
[all.rr.diag_sub_r2r3,all.rr.pval_r2r3] = corr_one2avgofothers(rgroup2,rgroup3,1,0,'r2r3'); % 1 = avg of all others
[all.rr.diag_sub_r2r1,all.rr.pval_r2r1] = corr_one2avgofothers(rgroup2,rgroup1,1,0,'r2r1'); % 1 = avg of all others
[all.rr.realdiff, all.rr.pvaldiff] = non_param_t(all.rr.diag_sub_r2r1,all.rr.diag_sub_r2r3,1); % based on vec1-vec2, type=1 paired,type2=unpaired
[all.rr.realdiff_z, all.rr.pvaldiff_z] = non_param_t(atanh(all.rr.diag_sub_r2r1),atanh(all.rr.diag_sub_r2r3),1); % based on vec1-vec2, type=1 paired,type2=unpaired
% 
% %%%%%%% RECALL 2&3 to MOVIE 1 %%%%%%%%%
[all.r23m.diag_sub_r2m1,all.r23m.pval_r2m1] = corr_one2avgofothers(rgroup2,mgroup1_avg,1,0,'r2m1'); % 1 = avg of all others
[all.r23m.diag_sub_r3m1,all.r23m.pval_r3m1] = corr_one2avgofothers(rgroup3,mgroup1_avg,1,0,'r3m1'); % 1 = avg of all others
[all.r23m.realdiff, all.r23m.pvaldiff] = non_param_t(all.r23m.diag_sub_r2m1,all.r23m.diag_sub_r3m1,2); % based on vec1-vec2, type=1 paired,type2=unpaired
[all.r23m.realdiff_z, all.r23m.pvaldiff_z] = non_param_t(atanh(all.r23m.diag_sub_r2m1),atanh(all.r23m.diag_sub_r3m1),2); % based on vec1-vec2, type=1 paired,type2=unpaired
% 
% %%%%%%% RECALL 2 to MOVIE 1&2 %%%%%%%%%
[all.r2m.diag_sub_r2m1,all.r2m.pval_r2m1] = corr_one2avgofothers(rgroup2,mgroup1_avg,3,0,'r2m1'); % 3 = avg others between
[all.r2m.diag_sub_r2m2,all.r2m.pval_r2m2] = corr_one2avgofothers(rgroup2,mgroup2_avg,2,0,'r2m2'); % 2 = avg others within
[all.r2m.realdiff, all.r2m.pvaldiff] = non_param_t(all.r2m.diag_sub_r2m1,all.r2m.diag_sub_r2m2,1); % based on vec1-vec2, type=1 paired,type2=unpaired
[all.r2m.realdiff_z, all.r2m.pvaldiff_z] = non_param_t(atanh(all.r2m.diag_sub_r2m1),atanh(all.r2m.diag_sub_r2m2),1); % based on vec1-vec2, type=1 paired,type2=unpaired

save(fullfile(basepath,'results','rsa','rois',['pattern_corr_imp_' roi_fnames '_' num2str(rr) '.mat']),'all');



