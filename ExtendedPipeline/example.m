% tile signal space
disp('tiling signal space...');
timeRange = 32;
frequencyRange = [64 1024];
qRange = [3.5 100];
sampleFrequency = 4096;
maximumMismatch = 0.1;
tiling = qtile(timeRange, qRange, frequencyRange, sampleFrequency, ...
               maximumMismatch);

% load detector data
disp('loading data...');
frameCacheFile = 'framecache.txt';
frameCache = loadframecache(frameCacheFile);
frameTypes = 'RDS_R_L3';
channelNames = 'H1:LSC-DARM_ERR';

%% Newest 12
% centerTime = 873910913.500000000; % 5.0 Mpc
% centerTime = 873350371.500000000; % 5.0 Mpc
% centerTime = 871356787.500000000; % 5.0 Mpc **
% centerTime = 868998719.500000000; % 5.0 Mpc
% centerTime = 868776347.500000000; % 5.0 Mpc
% centerTime = 868504980.500000000; % 5.0 Mpc
% centerTime = 867994624.500000000; % 5.0 Mpc *
% centerTime = 866088116.500000000; % 5.0 Mpc
% centerTime = 864872058.500000000; % 5.0 Mpc
% centerTime = 864358601.500000000; % 5.0 Mpc
% centerTime = 864282553.500000000; % 5.0 Mpc
% centerTime = 864224339.500000000; % 5.0 Mpc

%% Oldest 10
% centerTime = 841706022.500000000; % 5.0 Mpc
% centerTime = 841685117.500000000; % 5.0 Mpc
% centerTime = 841167503.500000000; % 5.0 Mpc ***
% centerTime = 840620331.500000000; % 5.0 Mpc
% centerTime = 840270025.500000000; % 5.0 Mpc
% centerTime = 839817610.500000000; % 5.0 Mpc ***
% centerTime = 839518469.500000000; % 5.0 Mpc ***
 centerTime = 838098134.500000000; % 5.0 Mpc ****
% centerTime = 835652713.500000000; % 5.0 Mpc
% centerTime = 833108974.500000000; % 5.0 Mpc

startTime = centerTime - timeRange / 2;
stopTime = centerTime + timeRange / 2;
[data, sampleFrequencies] = qreaddata(frameCache, channelNames, frameTypes, ...
                                      startTime, stopTime);

% resample data
disp('resampling data...');
data = qresample(data, sampleFrequencies, sampleFrequency);

% condition data
disp('conditioning data...');
data = qcondition(data, tiling);

% apply q transform
disp('transforming data...');
outlierFactor = 2.0;
transforms = qtransform(data, tiling, [], outlierFactor);

% threshold on significance
disp('thresholding on significance...');
falseRate = 1e0;
timeRange = 6 * [-1 +1] / 2;
frequencyRange = [64 512];
significants = qthreshold(transforms, tiling, startTime, falseRate, ...
                          centerTime, timeRange, frequencyRange);

% find most significant tile
[value, index] = max(significants{1}.normalizedEnergy);

% report SNR of most significant tile
fprintf(1, 'SNR of most significant tile: %.2f\n', sqrt(2 * (value - 1)));

% display spectrogram
disp('plotting spectrogram...');
qRange = significants{1}.q(index);
normalizedEnergyRange = [];
handle = qspectrogram(transforms, tiling, startTime, centerTime, ...
                      timeRange, frequencyRange, qRange, ...
                      normalizedEnergyRange, channelNames);
set(get(handle(1), 'Title'), 'String', '');
drawnow;
print -depsc2 001_injection.eps

% select events
disp('selecting non-overlapping tiles...');
durationInflation = 0.5;
bandwidthInflation = 0.5;
events = qselect(significants, durationInflation, bandwidthInflation);

% events{1}.normalizedEnergy = events{1}.normalizedEnergy * ...
%     durationInflation * bandwidthInflation;

% display eventgram
disp('plotting eventgram...');
normalizedEnergyRange = [];
handle = qeventgram(events, tiling, startTime, centerTime, ...
                    timeRange, frequencyRange, ...
                    durationInflation, bandwidthInflation, ...
                    normalizedEnergyRange, channelNames);
set(get(handle(1), 'Title'), 'String', '');
drawnow;
print -depsc2 002_tiles.eps

% measure distances
disp('measuring tile distances...');
tiles = [events{1}.time' events{1}.frequency' events{1}.q'];
distances = pdist(tiles, @qmetric);

% cluster tiles
disp('clustering tiles...');
tree = linkage(distances, 'single');
clusterDistance = 2.5;
clusters = cluster(tree, 'cutoff', clusterDistance, 'criterion', 'distance');

% sort clusters
disp('sorting clusters...');
uniqueClusters = sort(unique(clusters));
clusterNormalizedEnergy = zeros(size(uniqueClusters));
for clusterIndex = 1 : length(uniqueClusters),
  clusterNormalizedEnergy(clusterIndex) = ...
      sum(events{1}.normalizedEnergy(find(clusters == clusterIndex)));
end
[maximumClusterEnergy, maximumClusterIndex] = max(clusterNormalizedEnergy);
maximumClusterTileEnergy = max(events{1}.normalizedEnergy(find(clusters == maximumClusterIndex)));
maximumTileEnergy = max(events{1}.normalizedEnergy);

% display clustering
disp('plotting clusters...');
clf;
set(gca, 'FontSize', 16);
axis([timeRange log2(frequencyRange)]);
set(gca, 'YTick', log2(min(frequencyRange)) : log2(max(frequencyRange)));
set(gca, 'YTickLabel', 2.^(log2(min(frequencyRange)) : log2(max(frequencyRange))));
xlabel('Time [seconds]');
ylabel('Frequency [Hz]');
hold on;
uniqueClusters = unique(clusters);
markers = {'ro', 'gs', 'bv', 'c^', 'm<', 'k>', ...
           'rd', 'go', 'bs', 'cv', 'm^', 'k<', ...
           'r>', 'gd', 'bo', 'cs', 'mv', 'k^', ...
           'r<', 'g>', 'bd', 'co', 'ms', 'kv', ...
           'r^', 'g<', 'b>', 'cd', 'mo', 'ks', ...
           'rv', 'g^', 'b<', 'c>', 'md', 'ko', ...
           'rs', 'gv', 'b^', 'c<', 'm>', 'kd'};
% markers = {'ro', 'go', 'bo', 'co', 'mo', 'ko', ...
%            'rs', 'gs', 'bs', 'cs', 'ms', 'ks', ...
%            'rv', 'gv', 'bv', 'cv', 'mv', 'kv', ...
%            'r^', 'g^', 'b^', 'c^', 'm^', 'k^', ...
%            'r<', 'g<', 'b<', 'c<', 'm<', 'k<', ...
%            'r>', 'g>', 'b>', 'c>', 'm>', 'k>', ...
%            'rd', 'gd', 'bd', 'cd', 'md', 'kd'};
markers = [markers markers markers markers];
indices = find(clusters == maximumClusterIndex);
handle = plot(tiles(indices, 1) - centerTime, ...
              log2(tiles(indices, 2)), ...
              'ko');
set(handle, 'Color', 0.9 * [1 1 1]);
set(handle, 'MarkerSize', 30);
set(handle, 'MarkerFaceColor', get(handle, 'Color'));
set(handle, 'MarkerEdgeColor', 'none');
for clusterIndex = 1 : length(uniqueClusters),
  indices = find(clusters == uniqueClusters(clusterIndex));
  handle = plot(tiles(indices, 1) - centerTime, ...
                log2(tiles(indices, 2)), ...
                markers{clusterIndex});
  set(handle, 'MarkerSize', 10);
  set(handle, 'MarkerFaceColor', get(handle, 'Color'));
  set(handle, 'MarkerEdgeColor', 'none');
end
drawnow;
print -depsc2 003_hiclustered.eps

% measure distances
disp('measuring tile distances...');
distances = qdistance(events, durationInflation, bandwidthInflation);

% cluster tiles
disp('clustering tiles...');
clusteringRadius = 8;
numberThreshold = 4;
clusters = qcluster2(events, distances, clusteringRadius, numberThreshold);

clusters = clusters{1}.cluster;

% sort clusters
disp('sorting clusters...');
uniqueClusters = sort(unique(clusters));
clusterNormalizedEnergy = zeros(size(uniqueClusters));
for clusterIndex = 1 : length(uniqueClusters),
  clusterNormalizedEnergy(clusterIndex) = ...
      sum(events{1}.normalizedEnergy(find(clusters == clusterIndex)));
end
[maximumClusterEnergy, maximumClusterIndex] = max(clusterNormalizedEnergy);
maximumClusterTileEnergy = max(events{1}.normalizedEnergy(find(clusters == maximumClusterIndex)));
maximumTileEnergy = max(events{1}.normalizedEnergy);

% display clustering
disp('plotting clusters...');
clf;
set(gca, 'FontSize', 16);
axis([timeRange log2(frequencyRange)]);
set(gca, 'YTick', log2(min(frequencyRange)) : log2(max(frequencyRange)));
set(gca, 'YTickLabel', 2.^(log2(min(frequencyRange)) : log2(max(frequencyRange))));
xlabel('Time [seconds]');
ylabel('Frequency [Hz]');
hold on;
uniqueClusters = unique(clusters);
markers = {'ro', 'gs', 'bv', 'c^', 'm<', 'k>', ...
           'rd', 'go', 'bs', 'cv', 'm^', 'k<', ...
           'r>', 'gd', 'bo', 'cs', 'mv', 'k^', ...
           'r<', 'g>', 'bd', 'co', 'ms', 'kv', ...
           'r^', 'g<', 'b>', 'cd', 'mo', 'ks', ...
           'rv', 'g^', 'b<', 'c>', 'md', 'ko', ...
           'rs', 'gv', 'b^', 'c<', 'm>', 'kd'};
% markers = {'ro', 'go', 'bo', 'co', 'mo', 'ko', ...
%            'rs', 'gs', 'bs', 'cs', 'ms', 'ks', ...
%            'rv', 'gv', 'bv', 'cv', 'mv', 'kv', ...
%            'r^', 'g^', 'b^', 'c^', 'm^', 'k^', ...
%            'r<', 'g<', 'b<', 'c<', 'm<', 'k<', ...
%            'r>', 'g>', 'b>', 'c>', 'm>', 'k>', ...
%            'rd', 'gd', 'bd', 'cd', 'md', 'kd'};
markers = [markers markers markers markers];
indices = find(clusters == maximumClusterIndex);
handle = plot(tiles(indices, 1) - centerTime, ...
              log2(tiles(indices, 2)), ...
              'ko');
set(handle, 'Color', 0.9 * [1 1 1]);
set(handle, 'MarkerSize', 30);
set(handle, 'MarkerFaceColor', get(handle, 'Color'));
set(handle, 'MarkerEdgeColor', 'none');
for clusterIndex = 1 : length(uniqueClusters),
  indices = find(clusters == uniqueClusters(clusterIndex));
  handle = plot(tiles(indices, 1) - centerTime, ...
                log2(tiles(indices, 2)), ...
                markers{clusterIndex});
  set(handle, 'MarkerSize', 10);
  set(handle, 'MarkerFaceColor', get(handle, 'Color'));
  set(handle, 'MarkerEdgeColor', 'none');
end
drawnow;
print -depsc2 004_dbclustered.eps
