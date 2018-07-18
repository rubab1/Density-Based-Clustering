function process3_new()
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
confidenceInterval = 0.6826;

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
  set(gca, 'FontSize', 24);

  % set axis to semilogy
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'linear');

  % set standard axis
  axis([10^floor(log10(1 / livetime)) 1e-1 0 1]);

  % set figure size
  set(gca, 'Position', [0.13 0.14 0.77 0.78]);

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
  dataPoints = semilogx(clustersFalseRate, clustersEfficiency, 'ks-', ...
                        singlesFalseRate, singlesEfficiency, 'ko-');
  set(dataPoints(1), 'LineWidth', 1.5);
  set(dataPoints(1), 'MarkerSize', 8);
  set(dataPoints(1), 'Color', [0.75 0.75 0.75]);
  set(dataPoints(1), 'MarkerEdgeColor', get(dataPoints(1), 'Color'));
  set(dataPoints(1), 'MarkerFaceColor', get(dataPoints(1), 'Color'));
  set(dataPoints(2), 'LineWidth', 1.5);
  set(dataPoints(2), 'MarkerSize', 8);
  set(dataPoints(2), 'Color', [0.25 0.25 0.25]);
  set(dataPoints(2), 'MarkerEdgeColor', get(dataPoints(2), 'Color'));
  set(dataPoints(2), 'MarkerFaceColor', get(dataPoints(2), 'Color'));

  % label axes
  xlabel('False rate [Hz]');
  ylabel('Detection efficiency');
  title('Sinusoidal Gaussian signals with SNR of 10');

  % label curves
  labels = legend(dataPoints, 'clusters', 'singles', 2);
  set(labels, 'FontSize', 20);

  % display grid
  grid on;
  set(gca, 'XMinorGrid', 'off');
  set(gca, 'YMinorGrid', 'off');

  % flush figure
  drawnow;

  % export figure as eps image
  print('-depsc2', [waveform '_roc.eps']);

  % release figure hold
  hold off;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         plot event significance                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;

  % set font size
  set(gca, 'FontSize', 24);

  % set axis to loglog
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'log');

  % set figure size
  set(gca, 'Position', [0.13 0.14 0.77 0.78]);

  % retain axis properties
  hold on;

  % produce scatter plot
  dataPoints = loglog(noiseSinglesTiles, noiseSinglesEnergy, 'kv', ...
                      measuredSinglesTiles, measuredSinglesEnergy, 'k^', ...
                      noiseClustersTiles, noiseClustersEnergy, 'ko', ...
                      measuredClustersTiles, measuredClustersEnergy, 'ks');
  set(dataPoints(1), 'LineWidth', 1.5);
  set(dataPoints(1), 'MarkerSize', 8);
  set(dataPoints(1), 'Color', [0.00 0.00 0.00]);
  set(dataPoints(1), 'MarkerEdgeColor', get(dataPoints(1), 'Color'));
  set(dataPoints(1), 'MarkerFaceColor', get(dataPoints(1), 'Color'));
  set(dataPoints(2), 'LineWidth', 1.5);
  set(dataPoints(2), 'MarkerSize', 8);
  set(dataPoints(2), 'Color', [0.50 0.50 0.50]);
  set(dataPoints(2), 'MarkerEdgeColor', get(dataPoints(2), 'Color'));
  set(dataPoints(2), 'MarkerFaceColor', get(dataPoints(2), 'Color'));
  set(dataPoints(3), 'LineWidth', 1.5);
  set(dataPoints(3), 'MarkerSize', 8);
  set(dataPoints(3), 'Color', [0.25 0.25 0.25]);
  set(dataPoints(3), 'MarkerEdgeColor', get(dataPoints(3), 'Color'));
  set(dataPoints(3), 'MarkerFaceColor', get(dataPoints(3), 'Color'));
  set(dataPoints(4), 'LineWidth', 1.5);
  set(dataPoints(4), 'MarkerSize', 8);
  set(dataPoints(4), 'Color', [0.75 0.75 0.75]);
  set(dataPoints(4), 'MarkerEdgeColor', get(dataPoints(4), 'Color'));
  set(dataPoints(4), 'MarkerFaceColor', get(dataPoints(4), 'Color'));

  % set standard axis
  axis([1e0 1e3 1e0 1e4]);

  % label axes
  xlabel('Number of tiles in cluster');
  ylabel('Total cluster energy');
  title('Sinusoidal Gaussian signals with SNR of 10');

  % label trigger types
  labels = legend(dataPoints, ...
                  'noise singles', 'detected singles', ...
                  'noise clusters', 'detected clusters', ...
                  4);
  set(labels, 'FontSize', 20);

  % display grid
  grid on;
  % set(gca, 'XMinorGrid', 'off');
  % set(gca, 'YMinorGrid', 'off');

  % flush figure
  drawnow;

  % export figure as png image
  print('-depsc2', [waveform '_significance.eps']);

  % release figure hold
  hold off;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                plot false rate vs. significance threshold                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;

  % set font size
  set(gca, 'FontSize', 24);

  % set axis to loglog
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'log');

  % set figure size
  set(gca, 'Position', [0.13 0.14 0.77 0.78]);

  % retain axis properties
  hold on;

  % produce plot
  dataPoints = loglog(noiseClustersEnergy, clustersFalseRate, 'ks-', ...
                      noiseSinglesEnergy, singlesFalseRate, 'ko-');
  set(dataPoints(1), 'LineWidth', 1.5);
  set(dataPoints(1), 'MarkerSize', 8);
  set(dataPoints(1), 'Color', [0.75 0.75 0.75]);
  set(dataPoints(1), 'MarkerEdgeColor', get(dataPoints(1), 'Color'));
  set(dataPoints(1), 'MarkerFaceColor', get(dataPoints(1), 'Color'));
  set(dataPoints(2), 'LineWidth', 1.5);
  set(dataPoints(2), 'MarkerSize', 8);
  set(dataPoints(2), 'Color', [0.25 0.25 0.25]);
  set(dataPoints(2), 'MarkerEdgeColor', get(dataPoints(2), 'Color'));
  set(dataPoints(2), 'MarkerFaceColor', get(dataPoints(2), 'Color'));

  % label axes
  xlabel('Energy threshold');
  ylabel('False Rate');
  title('Sinusoidal Gaussian signals with SNR of 10');

  % label trigger types
  labels = legend(dataPoints, 'clusters', 'singles', 1);
  set(labels, 'FontSize', 20);

  % display grid
  grid on;
  set(gca, 'XTick', 10.^(0 : 1 : 4));
  set(gca, 'XMinorGrid', 'off');
  set(gca, 'YMinorGrid', 'off');

  % flush figure
  drawnow;

  % export figure as png image
  print('-depsc2', [waveform '_significance-falserate.eps']);

  % release figure hold
  hold off;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %           plot detection efficiency vs. significance threshold             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % clear figure
  clf;

  % set font size
  set(gca, 'FontSize', 24);

  % set axis to semilogx
  set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'linear');

  % set figure size
  set(gca, 'Position', [0.13 0.14 0.77 0.78]);

  % retain axis properties
  hold on;

  % produce plot
  dataPoints = semilogx(noiseClustersEnergy, clustersEfficiency, 'ks-', ...
                        noiseSinglesEnergy, singlesEfficiency, 'ko-');
  set(dataPoints(1), 'LineWidth', 1.5);
  set(dataPoints(1), 'MarkerSize', 8);
  set(dataPoints(1), 'Color', [0.75 0.75 0.75]);
  set(dataPoints(1), 'MarkerEdgeColor', get(dataPoints(1), 'Color'));
  set(dataPoints(1), 'MarkerFaceColor', get(dataPoints(1), 'Color'));
  set(dataPoints(2), 'LineWidth', 1.5);
  set(dataPoints(2), 'MarkerSize', 8);
  set(dataPoints(2), 'Color', [0.25 0.25 0.25]);
  set(dataPoints(2), 'MarkerEdgeColor', get(dataPoints(2), 'Color'));
  set(dataPoints(2), 'MarkerFaceColor', get(dataPoints(2), 'Color'));

  % label trigger types
  labels = legend(dataPoints, 'clusters', 'singles', 1);
  set(labels, 'FontSize', 20);

  % label axes
  xlabel('Energy threshold');
  ylabel('Detection efficiency');
  title('Sinusoidal Gaussian signals with SNR of 10');

  % display grid
  grid on;
  set(gca, 'XTick', 10.^(0 : 1 : 4));
  set(gca, 'XMinorGrid', 'off');
  set(gca, 'YMinorGrid', 'off');

  % flush figure
  drawnow;

  % export figure as png image
  print('-depsc2', [waveform '_significance-efficiency.eps']);

  % release figure hold
  hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           end loop over waveforms                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over waveforms
end

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
