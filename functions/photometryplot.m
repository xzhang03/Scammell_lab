function photometryplot(photometrypath, scoringpath, varargin)
%photometryplot is a teaching function to plot scammel lab photometry data
%   photometryplot(photometrypath, scoringpath, varargin)

% Take care of inputs or the lack of, so that you can just F5 the function
% and run it (better to call this function with a script)
if nargin < 3
    varargin = {};
    if nargin < 2
        scoringpath = '';
        if nargin < 1
            photometrypath = '';
        end
    end
end

% Parse input parameters
% The information below are defaults and will be overwritten by the
% varargin variable that is passed to here.
p = inputParser;

addOptional(p, 'defaultpath', 'D:\Data'); % Give default path for where the data are to be found
addOptional(p, 'photometryext', '*.csv'); % Photometry extension
addOptional(p, 'scoringext', '*.xlsx'); % Scoring extension
addOptional(p, 'transitions', {'W', 'NR', 1; 'W', 'NR', 10}); % A cell containing the transition info.

% Photometry variables
addOptional(p, 'tcolumn', 1); % time column
addOptional(p, 'ch1column', 2); % channel 1 column
addOptional(p, 'ch2column', 3); % channel 2 column
addOptional(p, 'usezscore', true); % zscore photometry
addOptional(p, 'zscoreperrow', false); % zscore per row (per transition)

% Plotting variables
addOptional(p, 'subplotrows', 6); % Proportion of heatmap vs average
addOptional(p, 'downsample', true); % Downsample photometry before plotting
addOptional(p, 'fs', 10); % Photometry fs to downsample to
addOptional(p, 'heatmaprange', [-0.1 0.1]); % Heatmap range
addOptional(p, 'linethickness', 1); % Line thickness

% Debugging variables
addOptional(p, 'scoringdata', []);
addOptional(p, 'scoringtext', {});
addOptional(p, 'photometrydata', []);

% Save figures
addOptional(p, 'savefig', true);
addOptional(p, 'closefigaftersave', false);


% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Input/output
% If no scoring path was given
if isempty(scoringpath)
    [fn, fp] = uigetfile(fullfile(p.defaultpath, p.scoringext), 'Choose scoring file.');
    scoringpath = fullfile(fp, fn);
end

% If not photometry path was given
if isempty(photometrypath)
    [fn, fp] = uigetfile(fullfile(p.defaultpath, p.photometryext), 'Choose photometry file.');
    photometrypath = fullfile(fp, fn);
end

% Show in terminal which file is being processed
[~, fnt, ~] = fileparts(photometrypath);
fprintf('======== Processing %s ========\n', fnt)

% Load both files
tic
fprintf('Loading... ')
if ~isempty(p.scoringdata) && ~isempty(p.scoringtext)
    scoringdata = p.scoringdata;
    scoringtext = p.scoringtext;
else
    [scoringdata, scoringtext, ~] = xlsread(scoringpath);
end

if ~isempty(p.photometrydata)
    photometrydata = p.photometrydata;
else
    [photometrydata, ~, ~] = xlsread(photometrypath);
end
t = toc;
fprintf('Done. Elapsed time = %i seconds.\n', round(t));

% Parse scoring
onsets = scoringdata(:,5);
durations = scoringdata(:,4);
labels = scoringtext(2:end,3);
labelpairs = [labels(1:end-1), labels(2:end)];

% Parse photometry
ch1 = photometrydata(:,p.ch1column);
ch2 = photometrydata(:,p.ch2column);
time = photometrydata(:,p.tcolumn);
l = size(photometrydata,1);

% Sampling rate
fs = 1/mode(diff(time));

% Downsample photometry data
if p.downsample
    tic
    fprintf('Downsampling... ')
    time = tcpBin(time, fs, p.fs, 'mean', 1);
    ch1 = tcpBin(ch1, fs, p.fs, 'mean', 1);
    ch2 = tcpBin(ch2, fs, p.fs, 'mean', 1);
    l = length(time);
    fs = p.fs;
    t = toc;
    fprintf('Done. Elapsed time = %i seconds.\n', round(t));
end

% zscore photometry
if p.usezscore
    ch1 = tcpZscore(ch1);
    ch2 = tcpZscore(ch2);
end

%% Loop through and make plots
% Number of plots
nplots = size(p.transitions, 1);

for iplot = 1 : nplots
    % Find the right type of sessions
    ipre = strcmpi(labelpairs(:,1), p.transitions{iplot, 1});
    ipost = strcmpi(labelpairs(:,2), p.transitions{iplot, 2});
    inds = find(ipre & ipost) + 1;
    
    % Onsets of the current transitions(in seconsd)
    consets = onsets(inds);
    nonsets = length(consets);
    
    % Durations of pre- and post-transition bouts
    dur_pre = durations(inds-1);
    dur_post = durations(inds);
    
    % Apply the window of interest to the durations
    dur_pre = min(dur_pre, p.transitions{iplot, 3});
    dur_post = min(dur_post, p.transitions{iplot, 3});
    
    % Turn seconds into indices
    inds_pre = round(dur_pre * fs);
    inds_post = round(dur_post * fs);
    consets_ind = zeros(nonsets,1);
    for i = 1 : nonsets
        [~, consets_ind(i)] = min(abs(time - consets(i)));
    end
    
    % Initialize matrix to load the data
    l_pre = round(p.transitions{iplot, 3} * fs);
    l_post = round(p.transitions{iplot, 3} * fs);
    l_all = l_pre + l_post + 1;
    ch1_mat = nan(nonsets, l_all);
    ch2_mat = nan(nonsets, l_all);
    
    % Loop through each transition and gather data
    for i = 1 : nonsets
        % Indices for current transition
        currind_start = consets_ind(i) - l_pre;
        currind_stop = consets_ind(i) + l_post;
        
        % Since not all transitions fit the lengths (some might be
        % shorter), we need to add empty numbers to the front and the end
        % if necessary
        if currind_start < 1
            prebuffer = nan(1, -currind_start + 1);
            r1 = horzcat(prebuffer, ch1(1 : currind_stop)');
            r2 = horzcat(prebuffer, ch2(1 : currind_stop)');
        elseif currind_stop > l
            postbuffer = nan(1, currind_stop - l);
            r1 = horzcat(ch1(currind_start : l)', postbuffer);
            r2 = horzcat(ch2(currind_start : l)', postbuffer);
        else
            r1 = ch1(currind_start : currind_stop)';
            r2 = ch2(currind_start : currind_stop)';
        end
        
        % Load each row into the datamatrix
        if p.zscoreperrow
            ch1_mat(i,:) = tcpZscore(r1);
            ch2_mat(i,:) = tcpZscore(r2);
        else
            ch1_mat(i,:) = r1;
            ch2_mat(i,:) = r2;
        end
        
    end
    
    %% Making figures
    % Initialize figure
    figure
        
    % Loop through
    for i = 1 : 2
        % Subplot for imagesc
        subplot(p.subplotrows, 2, (2 + i) : 2 : (2 * p.subplotrows));

        % Imagesc
        if i == 1
            imagesc(ch1_mat);
        else
            imagesc(ch2_mat);
        end
        
        % Color range
        if ~isempty(p.heatmaprange)
            colormap(b2r_arbitrary_input(p.heatmaprange(1), p.heatmaprange(2), [1 0 0], [0 0 1], [1 1 1]));
        else
            crange = caxis();
            colormap(b2r_arbitrary_input(crange(1), crange(2), [1 0 0], [0 0 1], [1 1 1]));
        end
        
        % Get xaxis range and ticks
        xrange = get(gca,'xlim');
        set(gca, 'XTick', [1, l_pre+1, l_pre+l_post]);
        
        % Set the xaxis right
        xlabel('Time (s)')
        set(gca,'XTickLabel', num2cell([-p.transitions{iplot, 3}, 0, p.transitions{iplot, 3}]));
        
        % Draw event lines
        hold on
        for j = 1 : nonsets
            % Draw pre-transition event line
            plot([l_pre-inds_pre(j), l_pre+1], [j j], 'g-', 'LineWidth', p.linethickness);
            
            % Draw post-transition event line
            plot([l_pre+1, l_pre+inds_post(j)], [j j], 'm-', 'LineWidth', p.linethickness);
        end
        
        % Draw a zero line
        plot([l_pre+1, l_pre+1], [0.5 nonsets+0.5], 'k-');
        hold off

        % Subplot 2
        % Subplot for average data
        subplot(p.subplotrows, 2, i)

        % Plot average data
        if i == 1
            datavec = nanmean(ch1_mat,1);
        else
            datavec = nanmean(ch2_mat,1);
        end
        plot(datavec);
        
        % Add an y = 0 line and an x = 0 line
        hold on
        plot(xrange, [0 0], 'Color', [0 0 0]);
        plot([l_pre+1, l_pre+1], [min(datavec), max(datavec)], 'k-');
        hold off

        % Align ranges
        xlim(xrange);
        set(gca, 'XTick', [1, l_pre+1, l_pre+l_post]);
        
        % Use channel 1 to set yrange
        if i == 1
            yrange = [min(datavec), max(datavec)];
        end
        ylim(yrange);

        set(gca,'XTickLabel',[]);
        
        % title
        title(sprintf('%s - %s', p.transitions{iplot, 1}, p.transitions{iplot, 2}))
    end
    
    % Save plots
    if p.savefig
        [fp, ~, ~] = fileparts(scoringpath);
        figfp = fullfile(fp, 'Plots');
        if ~exist(figfp, 'dir')
            mkdir(figfp);
        end
        
        % output figure
        figfn = sprintf('%s-%s_%is.png', p.transitions{iplot, 1}, p.transitions{iplot, 2},...
            p.transitions{iplot, 3});
        
        % save
        saveas(gcf, fullfile(figfp, figfn));
        
        if p.closefigaftersave
            close(gcf)
        end
    end
end

end

