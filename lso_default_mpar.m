function  [lso] = lso_default_mpar()
% [lso] = lso_default_param() loads mpar for lso

lso = [];
lso.fs = 100000;
lso.max_rate = 600;

type = 'coc_mid';
lso.(type).type = type;
lso.(type).data_name = 'ashida_2016';

% structure parameters
lso.(type).inputfiberLists_ipsi_ex =[];
lso.(type).inputfiberLists_contra_in =[];
lso.(type).inputfiberLists_ipsi_in =[];
lso.(type).inputfiberLists_contra_ex =[];
lso.(type).mapping_method = @uniform_within_cf_mapping;
lso.(type).cf = [];
lso.(type).neurons_per_cf = 1;

% parameters (ashida 2016)
lso.(type).Tref = 1.6e-3; % refractory time
lso.(type).ThEx = 8e-3;   % threshold
lso.(type).WiEx = 0.8e-3; % length coincidence window
lso.(type).gIn = 2;   % threshold increase by inhibition
lso.(type).WiIn = 1.6e-3; % length inhibition window
lso.(type).fibersPerNeuron_ipsi_ex = 20;
lso.(type).fibersPerNeuron_contra_in = 8;
lso.(type).fibersPerNeuron_ipsi_in = 0;
lso.(type).fibersPerNeuron_contra_ex = 0;
lso.(type).best_ipd = 0;


default_set = 'default';

if ~isfolder('param_store')
    mkdir('param_store')
end
save(['param_store/',default_set,'_lso.mat'],'lso')

end

