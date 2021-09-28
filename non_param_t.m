%% nonparametric t test
% type = 1 --> paired
% type = 2 --> non-paired
function [rdiff, pvaldiff,realdiff] = non_param_t(vec1,vec2,type)

if type==1
    realdiff = vec1 - vec2;
    for k=1:10000
        sgn = randsample([-1,1],length(vec1),true);
        shuffdiff(k) = mean(sgn.*realdiff);
    end
    
elseif type==2
    realdiff = mean(vec1) - mean(vec2);
    for k=1:10000
        vec = [vec1 vec2];
        vec = vec(randperm(length(vec)));
        vec1_rand = vec(1:length(vec1));
        vec2_rand = vec(length(vec1)+1:end);
        shuffdiff(k) = mean(vec1_rand) - mean(vec2_rand);
    end
end

pvaldiff = sum(shuffdiff>mean(realdiff))/k; % 1-tailed
rdiff = mean(realdiff);
end