% Inputs
photometrypath = 'I:\For Stephen from Chris 2\OX_DTR_GRAB_5HT_MCHERRY_8_Baseline_06152021AM2hrs_PhotometryMATLAB.csv';
scoringpath = 'I:\For Stephen from Chris 2\CC_OX_DTR_GRAB_5HT_MCHERRY_8_Baseline_AM_SS_Transitions.xlsx';

% Debug
% [scoringdata, scoringtext, ~] = xlsread(scoringpath);
% [photometrydata, ~, ~] = xlsread(photometrypath);
%%
% Define transitions
transitions = {'W', 'NR', [0,inf,0,60]; 'W', 'NR', [0,inf,60,600]; 'W', 'NR', [0,inf,600,inf];...
               'NR', 'W', [0,inf,0,60]; 'NR', 'W', [0,inf,60,600]; 'NR', 'R', [0,inf,600,inf];...
               'R', 'W', [0,inf,0,10]; 'NR', 'R', [0,inf,0,60]; 'NR', 'R', [0,inf,60,inf]};
           
% Paramters (read photometryplot.m for what they mean)
varargin = {'transitions', transitions, 'savefig', true, 'heatmaprange', [-4 4],...
    'downsample', true, 'fs', 5, 'zscoreperrow', true, 'closefigaftersave', true};

% Run the function
photometryplot(photometrypath, scoringpath, varargin)
