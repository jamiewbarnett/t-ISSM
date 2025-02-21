%% Project Description
% some shit...
%%
%% Toggles

smb_scale = 0;
melt_rate = 40; %Default of 40 m/yr
calving_threshold = 200e3; %Default 200e3

number_of_years = 100; %Length of simulation from 2300.

run_name = 'ice_go_byebye';


load(Spun up model.mat);

md = transientrestart(md);




md.timestepping.final_year = md.timestepping.start_year + number_of_years;

save(run_name,'md','-v7.3');









