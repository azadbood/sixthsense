%% pattern similarity correlation
% condition defines whether it's avg of all = 1 (averages all subj in the
% other group), avg of others within = 2 (takes out that subj out and
% calculates the average), avg of others between = 3 (randomizes the other
% group and takes one subject out and calculates the average) subject
% numbers in groups
% plt = 0 no plot, 1 plot
% ttl gets title as string input
% by Asieh Zadbood
function [diag_sub,pval_scn] = corr_one2avgofothers(one,avgofothers,cond,plt,ttl)
subnum = size(one,3)
subnum_oth = size(avgofothers,3)
avgofothers_rand = avgofothers(:,:,randperm(subnum_oth));
%% run corr each subj to the avg of others
for sub=1:subnum
    subject = one(:,:,sub);
    if cond == 1
        avgof = nanmean(avgofothers,3);
    elseif cond == 2
        others = setdiff(1:subnum,sub);
        avgof = nanmean(avgofothers(:,:,others),3);
    elseif cond == 3
        others = setdiff(1:subnum_oth,sub);
        avgof = nanmean(avgofothers_rand(:,:,others),3);
    end
    corrval = corr(subject,avgof);
    rsall(:,:,sub) = corrval;
    scenes(sub,:)=diag(corrval);
    diag_sub(sub) = nanmean(diag(corrval));
    offdiag_sub(sub) = (nansum(nansum(corrval))-nansum(diag(corrval)))/(size(corrval,1)*size(corrval,2)-length(diag(corrval)));
end
%% stats
% keep the subjects, shuffle the scene labels
scnum = size(one,2);
for k=1:1000
    for sub=1:subnum
        corrval_shuff = rsall(randperm(scnum),randperm(scnum),sub);
        diag_sub_shuff(sub) = nanmean(diag(corrval_shuff));
    end
    diag_allsub_shuff(:,k) = diag_sub_shuff;
end
mean(diag_sub)
pval_scn =sum(mean(diag_allsub_shuff)>mean(diag_sub))/k % 1-tailed

if plt==1
    % plot
    rsa_plot(diag_sub,offdiag_sub,pval_scn,ttl)
    %   rsa_plot(diag_sub,mean(diag_allsub_scn_shuff,2),pval_scn,ttl)
end

end