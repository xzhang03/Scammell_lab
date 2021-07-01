% Inputs
photometrypath = 'I:\For Stephen From Chris\CC_OX_DTR_GRAB_5HT_MCHERRY_8_AM_Photometry_data.csv';
scoringpath = 'I:\For Stephen From Chris\CC_OX_DTR_GRAB_5HT_MCHERRY_8_Baseline_AM_SS_OPT.xlsx';

% Debug
% [scoringdata, scoringtext, ~] = xlsread(scoringpath);
% [photometrydata, ~, ~] = xlsread(photometrypath);

% Define transitions
transitions = {'W', 'NR', 60; 'W', 'NR', 600; 'W', 'NR', 1800;...
               'NR', 'W', 60; 'NR', 'W', 600; 'NR', 'R', 1800;...
               'R', 'W', 10; 'NR', 'R', 60; 'NR', 'R', 600};
           
% Paramters (read photometryplot.m for what they mean)
varargin = {'transitions', transitions, 'savefig', true, 'heatmaprange', [-4 4],...
    'downsample', true, 'fs', 2, 'zscoreperrow', true, 'closefigaftersave', true};

% Run the function
photometryplot(photometrypath, scoringpath, varargin)
