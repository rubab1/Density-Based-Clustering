function process3()
% PROCESS3 Plot ROC curves from QTTESTCLUSTER results
%
% PROCESS3 plots receiver operator characteristic (ROC) curves
% and scatter plots of trigger significance in order to compare
% the performance of clustering algorithms on various waveforms
% as tested by QTESTCLUSTER.
%
% See also QTESTCLUSTER.

% Shourov K. Chatterji
% shourov@ligo.caltech.edu
% 2006-Jul-24

% $Id:$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  parameters                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% confidence interval for error bars
confidenceInterval = 0.9;

% read number of block and block livetime from file
[numberOfBlocks, blockLivetime] = dataread('file', 'livetime.txt', '%f %f');

% total analyzed livetime
livetime = numberOfBlocks * blockLivetime;

% read detection window duration from process2.sh
detectionDuration = load('detection.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              initialize display                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open figure window
clf;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              create html index                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to html index file
indexFile = 'index.html';

% open html index file
indexFileID = fopen(indexFile, 'w');

% initialize html index file
fprintf(indexFileID, '<html>\n');
fprintf(indexFileID, '<head>\n');
fprintf(indexFileID, '<title>QTestCluster: %s</title>\n', ...
        strrep(regexprep(pwd, '^.*/', ''), '_', ' '));
fprintf(indexFileID, '</head>\n');
fprintf(indexFileID, '<body>\n');
fprintf(indexFileID, '<h2>QTestCluster: %s</h2>\n', ...
        strrep(regexprep(pwd, '^.*/', ''), '_', ' '));
fprintf(indexFileID, '<p>\n');
fprintf(indexFileID, '<a href="parameters.txt">parameters</a> | ');
fprintf(indexFileID, '<a href="log.txt">log</a>\n');
fprintf(indexFileID, '</p>\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        results file column defintions                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% results file column index for central time
timeColumn = 1;

% results file column index for duration
durationColumn = 2;

% results file column index for central frequency
frequencyColumn = 3;

% results file column index for bandwidth
bandwidthColumn = 4;

% results file column index for normalized energy
energyColumn = 5;

% results file column index for number of tiles
tilesColumn = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             load noise triggers                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load unclustered noise triggers
noiseSingles = load('noise_singles_sparse.txt');

% load clustered noise triggers
noiseClusters = load('noise_clusters_sparse.txt');

% extract unclustered noise trigger total energy
noiseSinglesEnergy = noiseSingles(:, energyColumn);

% extract clustered noise trigger total energy
noiseClustersEnergy = noiseClusters(:, energyColumn);

% extract unclustered noise trigger number of tiles
noiseSinglesTiles = noiseSingles(:, tilesColumn);

% extract clustered noise trigger number of tiles
noiseClustersTiles = noiseClusters(:, tilesColumn);

% extract unclustered noise trigger duration
noiseSinglesDuration = noiseSingles(:, durationColumn);

% extract clustered noise trigger duration
noiseClustersDuration = noiseClusters(:, durationColumn);

% determine significance of unclustered noise triggers
% ********** CONSIDER ALTERNATIVE DEFINITION OF SIGNIFICANCE **********
noiseSinglesSignificance = noiseSinglesEnergy;

% determine significance of clustered noise triggers
% ********** CONSIDER ALTERNATIVE DEFINITION OF SIGNIFICANCE **********
noiseClustersSignificance = noiseClustersEnergy;

% sort unclustered noise triggers by decreasing significance
[ignore, noiseSinglesSortedIndices] = sort(noiseSinglesSignificance);
noiseSinglesSortedIndices = flipud(noiseSinglesSortedIndices);
noiseSinglesEnergy = ...
    noiseSinglesEnergy(noiseSinglesSortedIndices);
noiseSinglesTiles = ...
    noiseSinglesTiles(noiseSinglesSortedIndices);
noiseSinglesDuration = ...
    noiseSinglesDuration(noiseSinglesSortedIndices);
noiseSinglesSignificance = ...
    noiseSinglesSignificance(noiseSinglesSortedIndices);

% sort clustered noise triggers by significance
[ignore, noiseClustersSortedIndices] = sort(noiseClustersSignificance);
noiseClustersSortedIndices = flipud(noiseClustersSortedIndices);
noiseClustersEnergy = ...
    noiseClustersEnergy(noiseClustersSortedIndices);
noiseClustersTiles = ...
    noiseClustersTiles(noiseClustersSortedIndices);
noiseClustersDuration = ...
    noiseClustersDuration(noiseClustersSortedIndices);
noiseClustersSignificance = ...
    noiseClustersSignificance(noiseClustersSortedIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            get list of waveforms                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read list of injection files
injectionFiles = dir('*_injections.txt');

% number of waveforms
numberOfWaveforms = length(injectionFiles);

% extract waveform names from injection file names
waveforms = cell(1, numberOfWaveforms);
for waveformNumber = 1 : numberOfWaveforms,,
  waveforms{waveformNumber} = strrep(injectionFiles(waveformNumber).name, ...
                                     '_injections.txt', '');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             loop over waveforms                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over waveforms
for waveformNumber = 1 : numberOfWaveforms,

  % extract waveform name
  waveform = waveforms{waveformNumber};

  % display status
  fprintf(1, 'processing %s...\n', waveform);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         load waveform properties                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % load injection properties
  injections = load([waveform '_injections.txt']);
  
  % extract injection durations
  injectionsDuration = injections(:, 2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          load waveform triggers                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % load properties of injections detected without clustering
  measuredSingles = load([waveform '_singles_measured.txt']);

  % load properties of injections detected with clustering
  measuredClusters = load([waveform '_clusters_measured.txt']);

  % extract energy of injections detected without clustering
  measuredSinglesEnergy = measuredSingles(:, energyColumn);

  % extract energy injections detected with clustering
  measuredClustersEnergy = measuredClusters(:, energyColumn);

  % extract number of tiles of injections detected without clustering
  measuredSinglesTiles = measuredSingles(:, tilesColumn);

  % extract number of tiles of injections detected with clustering
  measuredClustersTiles = measuredClusters(:, tilesColumn);

  % extract duration of injections detected without clustering
  measuredSinglesDuration = measuredSingles(:, durationColumn);

  % extract duration of injections detected with clustering
  measuredClustersDuration = measuredClusters(:, durationColumn);

  % determine significance of injections detected without clustering
  % ********** CONSIDER ALTERNATIVE DEFINITION OF SIGNIFICANCE **********
  measuredSinglesSignificance = measuredSinglesEnergy;

  % determine significance of injections detected with clustering
  % ********** CONSIDER ALTERNATIVE DEFINITION OF SIGNIFICANCE **********
  measuredClustersSignificance = measuredClustersEnergy;

  % sort injections detected without clustering by decreasing significance
  [ignore, measuredSinglesSortedIndices] = sort(measuredSinglesSignificance);
  measuredSinglesSortedIndices = flipud(measuredSinglesSortedIndices);
  measuredSinglesEnergy = ...
      measuredSinglesEnergy(measuredSinglesSortedIndices);
  measuredSinglesTiles = ...
      measuredSinglesTiles(measuredSinglesSortedIndices);
  measuredSinglesDuration = ...
      measuredSinglesDuration(measuredSinglesSortedIndices);
  measuredSinglesSignificance = ...
      measuredSinglesSignificance(measuredSinglesSortedIndices);

  % sort injections detected without clustering by decreasing significance
  [ignore, measuredClustersSortedIndices] = sort(measuredClustersSignificance);
  measuredClustersSortedIndices = flipud(measuredClustersSortedIndices);
  measuredClustersEnergy = ...
      measuredClustersEnergy(measuredClustersSortedIndices);
  measuredClustersTiles = ...
      measuredClustersTiles(measuredClustersSortedIndices);
  measuredClustersDuration = ...
      measuredClustersDuration(measuredClustersSortedIndices);
  measuredClustersSignificance = ...
      measuredClustersSignificance(measuredClustersSortedIndices);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                          determine false rates                             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % observed number of unclustered noise events and confidence interval
  [singlesFalseEvents, singlesFalseEventsInterval] = ...
      poissfit((1 : length(noiseSinglesSignificance)), 1 - confidenceInterval);

  % observer number of clustered noise events and confidence interval
  [clustersFalseEvents, clustersFalseEventsInterval] = ...
      poissfit((1 : length(noiseClustersSignificance)), 1 - confidenceInterval);

  % estimated unclustered false rate and confidence interval
  singlesFalseRate = singlesFalseEvents' / livetime;
  singlesFalseRateInterval = singlesFalseEventsInterval' / livetime;

  % estimated clustered false rate and confidence interval
  clustersFalseRate = clustersFalseEvents' ./ livetime;
  clustersFalseRateInterval = clustersFalseEventsInterval' / livetime;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      determine detection efficiency                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % number of injections detected by unclustered triggers
  singlesDetected = ...
      sum((ones(length(noiseSinglesSignificance), 1) * ...
           measuredSinglesSignificance') > ...
          (noiseSinglesSignificance * ...
           ones(1, length(measuredSinglesSignificance))), ...
          2);

  % number of injections detected by clustered triggers
  clustersDetected = ...
      sum((ones(length(noiseClustersSignificance), 1) * ...
           measuredClustersSignificance') > ...
          (noiseClustersSignificance * ...
           ones(1, length(measuredClustersSignificance))), ...
          2);

  % esimated detection efficiency and confidence interval without clustering
  [singlesEfficiency, singlesEfficiencyInterval] = ...
      binofit(singlesDetected, numberOfBlocks, 1 - confidenceInterval);

  % esimated detection efficiency and confidence interval with clustering
  [clustersEfficiency, clustersEfficiencyInterval] = ...
      binofit(clustersDetected, numberOfBlocks, 1 - confidenceInterval);

  % accidental probability of detection without clustering
  singlesAccidentalProbability = 1 - ...
      exp(-2 * singlesFalseRate * detectionDuration) .* ...
      mean(exp(-singlesFalseRate * injectionsDuration'), 2) .* ...
      diag(cumsum(exp(-singlesFalseRate * noiseSinglesDuration'), 2)) ./ ...
      (1 : length(singlesFalseRate))';

  % accidental probability of detection with clustering
  clustersAccidentalProbability = 1 - ...
      exp(-2 * clustersFalseRate * detectionDuration) .* ...
      mean(exp(-clustersFalseRate * injectionsDuration'), 2) .* ...
      diag(cumsum(exp(-clustersFalseRate * noiseClustersDuration'), 2)) ./ ...
      (1 : length(clustersFalseRate))';
  
  % corrected detection efficiency without clustering
  singlesEfficiency = ...
      max(singlesEfficiency - singlesAccidentalProbability, 0) ./ ...
      (1 - singlesAccidentalProbability);
  
  % corrected detection efficiency with clustering
  clustersEfficiency = ...
      max(clustersEfficiency - clustersAccidentalProbability, 0) ./ ...
      (1 - clustersAccidentalProbability);
  
  % corrected detection efficiency confidence interval without clustering
  singlesAccidentalProbability = singlesAccidentalProbability * [1 1];
  singlesEfficiencyInterval = ...
      max(singlesEfficiencyInterval - singlesAccidentalProbability, 0) ./ ...
      (1 - singlesAccidentalProbability);

  % corrected detection efficiency confidence interval with clustering
  clustersAccidentalProbability = clustersAccidentalProbability * [1 1];
  clustersEfficiencyInterval = ...
      max(clustersEfficiencyInterval - clustersAccidentalProbability, 0) ./ ...
      (1 - clustersAccidentalProbability);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                            write roc results                               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % path to unclustered roc results file
  singlesROCFile = [waveform '_singles_roc.txt'];

  % path to clustered roc results file
  clustersROCFile = [waveform '_clusters_roc.txt'];
  
  % open unclustered roc results file
  singlesROCFileID = fopen(singlesROCFile, 'w');

  % open clustered roc results file
  clustersROCFileID = fopen(clustersROCFile, 'w');

  % write unclustered roc results
  fprintf(singlesROCFileID, ...
          '%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n', ...
          [singlesFalseRate singlesFalseRateInterval ...
           singlesEfficiency singlesEfficiencyInterval]');

  % write clustered roc results
  fprintf(clustersROCFileID, ...
          '%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n', ...
          [clustersFalseRate clustersFalseRateInterval ...
           clustersEfficiency clustersEfficiencyInterval]');

  % close unclustered roc results file
  fclose(singlesROCFileID);

  % close unclustered roc results file
  fclose(clustersROCFileID);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                             plot roc curves                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;

  % set font size
  set(gca, 'FontSize', 18);

  % set axis to semilogy
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'linear');

  % set standard axis
  axis([10^floor(log10(1 / livetime)) 1e-1 0 1]);

  % retain axis properties
  hold on;

  % plot detection effiency confidence intervals
  verticalErrorBars = ...
      semilogx([singlesFalseRate singlesFalseRate]', ...
               singlesEfficiencyInterval', ...
               [clustersFalseRate clustersFalseRate]', ...
               clustersEfficiencyInterval');
  set(verticalErrorBars, 'Color', 0.75 * [1 1 1]);
  set(verticalErrorBars, 'LineWidth', 1);

  % plot false rate confidence intervales
  horizontalErrorBars = ...
      semilogx(singlesFalseRateInterval', ...
               [singlesEfficiency singlesEfficiency]', ...
               clustersFalseRateInterval', ...
               [clustersEfficiency clustersEfficiency]');
  set(horizontalErrorBars, 'Color', 0.75 * [1 1 1]);
  set(horizontalErrorBars, 'LineWidth', 1);

  % plot data points
  dataPoints = semilogx(singlesFalseRate, singlesEfficiency, 'r.-', ...
                        clustersFalseRate, clustersEfficiency, 'b.-');
  set(dataPoints, 'MarkerSize', 10);

  % label axes
  xlabel('False rate [Hz]');
  ylabel('Detection efficiency');
  title([strrep(waveform, '_', ' ') ' ROC']);

  % label curves
  labels = legend(dataPoints, 'singles', 'clusters', 2);
  set(labels, 'FontSize', 16);

  % flush figure
  drawnow;

  % export figure as png image
  print('-dpng', [waveform '_roc.png']);

  % produce thumbnail image
  unix(['convert -resize 400x400 ' waveform '_roc.png ' ...
        waveform '_roc_small.png']);

  % append figure to html index file
  fprintf(indexFileID, ...
          '<a href="%s_roc.png"><img src="%s_roc_small.png"></a>\n', ...
          waveform, waveform);
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         plot event significance                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;
  
  % set font size
  set(gca, 'FontSize', 18);

  % produce scatter plot
  dataPoints = ...
      loglog(noiseSinglesTiles, noiseSinglesEnergy, 'm.', ...
             noiseClustersTiles, noiseClustersEnergy, 'r.', ...
             measuredSinglesTiles, measuredSinglesEnergy, 'c.', ...
             measuredClustersTiles, measuredClustersEnergy, 'b.');
  set(dataPoints, 'Marker', 'o');
  set(dataPoints, 'MarkerSize', 6);

  % set standard axis
  axis([1e0 1e3 1e0 1e4]);

  % label axes
  xlabel('Number of tiles in cluster');
  ylabel('Total cluster energy');
  title([strrep(waveform, '_', ' ') ' significance']);

  % label trigger types
  labels = legend(dataPoints, ...
                  'noise singles', 'noise clusters', ...
                  'detected singles', 'detected clusters', ...
                  4);
  set(labels, 'FontSize', 16);

  % flush figure
  drawnow;

  % export figure as png image
  print('-dpng', [waveform '_significance.png']);

  % produce thumbnail image
  unix(['convert -resize 400x400 ' waveform '_significance.png ' ...
        waveform '_significance_small.png']);

  % append figure to html index file
  fprintf(indexFileID, ...
          ['<a href="%s_significance.png">' ...
           '<img src="%s_significance_small.png"></a>\n<br />'], ...
          waveform, waveform);

  hold off;
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                             plot roc breakdown                             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;
  drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  hold on;

  % set font size
  set(gca, 'FontSize', 18);
  % set axis to loglog
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'log');

  % produce plot
  dataPoints = ...
      loglog(noiseSinglesEnergy, singlesFalseRate, 'm.',...
      noiseClustersEnergy, clustersFalseRate, 'c.');
  set(dataPoints, 'Marker', 'o');
  set(dataPoints, 'MarkerSize', 6);

  % label axes
  xlabel('Energy');
  ylabel('False Rate');
  title([strrep(waveform, '_', ' ') ' significance-falserate']);
  
  % label trigger types
  labels = legend(dataPoints, 'singles', 'clusters', 2);
  set(labels, 'FontSize', 16);

  % flush figure
  drawnow;

  hold off;
  
  % export figure as png image
  print('-dpng', [waveform '_significance-falserate.png']);

  % produce thumbnail image
  unix(['convert -resize 400x400 ' waveform '_significance-falserate.png ' ...
        waveform '_significance-falserate_small.png']);
  
  % append figure to html index file
  fprintf(indexFileID, ...
          ['<a href="%s_significance-falserate.png">' ...
           '<img src="%s_significance-falserate_small.png"></a>\n'], ...
          waveform, waveform);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;

  hold on;
  
  % set font size
  set(gca, 'FontSize', 18);
  % set axis to semilogx
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'linear');

  % produce plot
  dataPoints = semilogx(noiseSinglesEnergy, singlesEfficiency, 'r.',...
      noiseClustersEnergy, clustersEfficiency, 'b.');
  set(dataPoints, 'Marker', 'o');
  set(dataPoints, 'MarkerSize', 6);
  % label trigger types
  labels = legend(dataPoints, 'singles', 'clusters', 2);
  set(labels, 'FontSize', 16);
  
  % label axes
  xlabel('Energy');
  ylabel('Efficiency');
  title([strrep(waveform, '_', ' ') ' significance-efficiency']);

  hold off;
  
  % flush figure
  drawnow;

  % export figure as png image
  print('-dpng', [waveform '_significance-efficiency.png']);

  % produce thumbnail image
  unix(['convert -resize 400x400 ' waveform '_significance-efficiency.png ' ...
        waveform '_significance-efficiency_small.png']);

  % append figure to html index file
  fprintf(indexFileID, ...
          ['<a href="%s_significance-efficiency.png">' ...
           '<img src="%s_significance-efficiency_small.png"></a>\n<br />'], ...
          waveform, waveform);
         
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           end loop over waveforms                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over waveforms
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            close html index file                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% finalize html index file
fprintf(indexFileID, '<hr>\n');
fprintf(indexFileID, 'Created by user %s on %s at %s<br />\n', ...
        getenv('USER'), datestr(clock, 29), datestr(clock, 13));
fprintf(indexFileID, '</body>\n');
fprintf(indexFileID, '</html>\n');

% close html index file
fclose(indexFileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              modified POISSFIT                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambdahat, lambdaci] = poissfit(x,alpha)
%POISSFIT Parameter estimates and confidence intervals for Poisson data.
%   POISSFIT(X) Returns the estimate of the parameter of the Poisson
%   distribution give the data X.
%
%   [LAMBDAHAT, LAMBDACI] = POISSFIT(X,ALPHA) gives MLEs and 100(1-ALPHA)
%   percent confidence intervals given the data. By default, the
%   optional parameter ALPHA = 0.05 corresponding to 95% confidence
%   intervals.
%
%   See also POISSCDF, POISSINV, POISSPDF, POISSRND, POISSTAT, MLE.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 2.9 $  $Date: 2002/01/17 21:31:38 $

% ------------------------------------------------------------------------------
% This is a modified version of poissfit that can process a vector of single
% measurements.
%
% Shourov K. Chatterji
% shourov@ligo.caltech.edu
% 2006-Jul-24
% ------------------------------------------------------------------------------

if nargin < 2
    alpha = 0.05;
end

% Initialize params to zero.
[m, n] = size(x);
% if min(m,n) == 1
%    x = x(:);
%    m = max(m,n);
%    n = 1;
% end

lambdahat = mean(x, 1);

if nargout > 1,
   lsum = m*lambdahat;
   k = find(lsum < 100);
   if any(k)
      lb(k) = chi2inv(alpha/2, 2*lsum(k))/2;
      ub(k) = chi2inv(1-alpha/2, 2*(lsum(k)+1))/2;
   end
   k = find(lsum >= 100);
   if any(k)
      lb(k) = norminv(alpha/2,lsum(k),sqrt(lsum(k)));
      ub(k) = norminv(1 - alpha/2,lsum(k),sqrt(lsum(k)));
   end

   lambdaci = [lb;ub]/m;
end
