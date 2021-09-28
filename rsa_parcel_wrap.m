%% code to collect pattern similarity results calculated by "rsa_parcel.m"
% cretaes maps (nii) with correlation (r) values in DMN regions, nan in
% non-DMN regions
% saves FDR values and other info

clear
rand('state', sum(100*clock));

savenames = {'mm','rr','r23m','r2m','interact'};

roi_fnames = 'Schaefer_100_icm152';
roi_num = 100;

basepath = 'Y:\Asieh\MovieTwistProject\'; addpath(basepath);
addpath('Y:\Asieh\isc_nifti_kit\NIFTI_tools');addpath('Y:\Asieh\isc_nifti_kit\');
addpath('Y:\Asieh\MovieTwistProject\funcs\');

fid = fopen(fullfile(basepath,'rois','Shaefer2018','Parcellations','MNI',['Schaefer2018_' num2str(roi_num) 'Parcels_7Networks_order.txt']));%testnames.txt')); % all or clean
data = textscan(fid,'%s%s%s%s%s%s','HeaderLines',0,'CollectOutput',1);
data = data{:};
fid = fclose(fid);
dm_parcels = find(cellfun('isempty',strfind(data(:,2),'Default'))==0);

for cond = 1:length(savenames)
    savename = savenames{cond};
    
    roi_mask = load_nii(fullfile(basepath,'MovieTwistProject','masks',[roi_fnames '.nii']));
    mask = reshape(roi_mask.img,[],1);

    mask_rdiff_dm_z = nan(size(reshape(roi_mask.img,[],1)));
    mask_pval_dm_z = nan(size(reshape(roi_mask.img,[],1)));
    
    test = [];
    for thisROI=1:roi_num
        
        roi_bin = find(mask==thisROI);
        
        load(fullfile(basepath,'MovieTwistProject','results','rsa','roiout',['pattern_corr_imp_' roi_fnames '_' num2str(thisROI) '.mat']));
        
        if strcmp(savename,'interact')
            [rdiff_z(thisROI), pval_z(thisROI)] = non_param_t(atanh(all.mm.diag_sub_m2m3-all.mm.diag_sub_m2m1),atanh(all.rr.diag_sub_r2r3-all.rr.diag_sub_r2r1),1); % based on vec1-vec2, type=1 paired,type2=unpaired
        else
            rdiff_z(thisROI) = all.(savename).realdiff_z;
            pval_z(thisROI) = all.(savename).pvaldiff_z;            
        end
                                
        if ismember(thisROI,dm_parcels)
            mask_rdiff_dm_z(roi_bin)= rdiff_z(thisROI);
            mask_pval_dm_z(roi_bin)= pval_z(thisROI);
        end
        
        clear roi_bin
        
    end
    
    savenames(cond);

% FDR
    rdiff_dm_z = rdiff_z(dm_parcels);
    pval_dm_z = pval_z(dm_parcels);
    
    allrealdiff_dm{cond} = rdiff_dm_z;
    allpvaldiff_dm{cond} = pval_dm_z;
    
    fdr_dm(cond) = fdr_BH(pval_dm_z, 0.05);
    fdr_effect_dm{1,cond} = rdiff_dm_z(pval_dm_z<=fdr_dm(cond));
    fdr_effect_dm{2,cond} = pval_dm_z(pval_dm_z<=fdr_dm(cond));
    effect_dm{cond} = rdiff_dm_z(pval_dm_z<0.05);
    
    clear rdiff_dm rdiff_dm_z pval_z pval_dm_z
%     
    %% saving fdr realdiff, everything else nan
    
    mask_vox_nan = mask_rdiff_dm_z; % take dmn filled matrix

    if isnan(fdr_dm(cond))
        mask_vox_nan(mask_pval_dm_z > 0.05) = nan; % if nothing passes fdr, put >0.05 = nan
    else 
        mask_vox_nan(mask_pval_dm_z > fdr_dm(cond)) = nan; % if >fdr put nans
    end
    nii = load_nii(fullfile(basepath,'MovieTwistProject','masks',[roi_fnames '.nii']));
    img = nii.img;
    cubenii = nii;
    cubenii.img = reshape(mask_vox_nan,size(img));
    save_nii(cubenii,fullfile(basepath,'MovieTwistProject','results','rsa','maps',[savename '_pattern_corr_imp_' num2str(thisROI) '_realdiff_dm.nii']));
    clear mask_vox_nan
end

save(fullfile(basepath,'MovieTwistProject','results','rsa','maps',['pattern_corr_imp_' roi_fnames '.mat']),...
    'allrealdiff_dm','allpvaldiff_dm','fdr_effect_dm','fdr_dm','savenames','dm_parcels','effect_dm');

