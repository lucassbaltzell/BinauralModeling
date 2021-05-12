function cf_fib = getCFfib(an_cfs,stim_cfs,bw_ext,stimtype,stimpar,eqflg)
%This function returns a matrix of AN cfs, with each row corresponding to a
%stimulus n. Each stimulus n receives an equal number of AN cfs
%OUPTUT
%cf_fib = matrix of desired AN cf fibers for each stimulus cf
%INPUTS
%an_cfs = vector of all possible AN cfs
%stim_cfs = vector of stimulus cfs
%bw_ext = amount to extend bandwidth (octave units)
%stimtype = type of stimulus
%stimpar = stimulus paramters
%eqflg = 0 or 1. if 1, ensure equal number of AN fibers for each stimulus cf

%created by Luke Baltzell 05/11/21

S = length(stim_cfs);
%get number of fibers (nfib) for each stimulus cf (stim_cf)
cf_fib = cell(1,S);
flims_ext = cell(1,S);
for s = 1:S
    if strcmp(stimtype,'nbNoise') == 1
        [~,flims] = genNBnoise(stimpar.dur,stimpar.fs,stim_cfs(s),...
            stimpar.bw,stimpar.tflg);
    elseif strcmp(stimtype,'rustleNoise') == 1
        [~,flims] = genNBrustle(stimpar.dur,stimpar.fs,stimpar.gpwidth,...
            stim_cfs(s),stimpar.bw,stimpar.tflg,stimpar.nord);
    elseif strcmp(stimtype,'pureTone') == 1
        flims = [stim_cfs(s) stim_cfs(s)];
    end
    flims_ext{1,s} = [flims(1)/2^(1/(2/bw_ext)) flims(2)*2^(1/(2/bw_ext))]; %extend bw
    [~,ind] = find(flims_ext{1,s}(1) < an_cfs & flims_ext{1,s}(2) > an_cfs);
    cf_fib{1,s} = an_cfs(ind);
    nfibs(s) = length(ind);
end

if sum(nfibs - max(nfibs)) ~= 0
    if eqflg == 1
       ind = find(nfibs < max(nfibs));
       for s = 1:length(ind)
           flg = 0;
           while flg == 0
               lc = flims_ext{1,s}(1);
               uc = flims_ext{1,s}(2);
               li = find(an_cfs == min(cf_fib{1,s}));
               ui = find(an_cfs == max(cf_fib{1,s}));
               if an_cfs(li)/lc > uc/an_cfs(ui)
                   %add lower cf
                   cf_fib{1,s} = cat(2,an_cfs(li-1),cf_fib{1,s});
                   lc = an_cfs(li-1);
               elseif an_cfs(li)/lc < uc/an_cfs(ui)
                   %add upper cf
                   cf_fib{1,s} = cat(2,cf_fib{1,s},an_cfs(ui+1));
                   uc = an_cfs(ui+1);
               elseif an_cfs(li)/lc == uc/an_cfs(ui)
                   %add upper cf
                   cf_fib{1,s} = cat(2,cf_fib{1,s},an_cfs(ui+1));
                   uc = an_cfs(ui+1);
               end
               if length(cf_fib{1,s}) == max(nfibs)
                   flg = 1;
               end
           end
       end
    end
end
                   
               