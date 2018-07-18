


% sort clusters
disp('sorting clusters...');
uniqueClusters = sort(unique(tempClusters));


% Read cluster info
for clusterIndex = 1:clusters{channelNumber}.numberOfClusters,
    
    % Read time, frequency,and normalized energy info of each tile
    i = find(tempClusters == clusterIndex);
    t = events{channelNumber}.time(i);
    f = events{channelNumber}.frequency(i);
    z = events{channelNumber}.normalizedEnergy(i);

    % Total number of tiles and sum of tile-energy
    clusters{channelNumber}.cluster{clusterIndex}.n = length(i);
    clusters{channelNumber}.cluster{clusterIndex}.z.total = sum(z);

    % Information about "most significant tile" of the cluster
    iMax = find((z == max(z)), 1, 'first');
    clusters{channelNumber}.cluster{clusterIndex}.t.mst = t(iMax);
    clusters{channelNumber}.cluster{clusterIndex}.f.mst = f(iMax);
    clusters{channelNumber}.cluster{clusterIndex}.z.mst = z(iMax);

    % Info about the cluster     
    clusters{channelNumber}.cluster{clusterIndex}.t.start = min(t);
    clusters{channelNumber}.cluster{clusterIndex}.t.end = max(t);
    clusters{channelNumber}.cluster{clusterIndex}.t.middle = (min(t) + max(t)) / 2;
    clusters{channelNumber}.cluster{clusterIndex}.t.duration = max(t) - min(t);
    clusters{channelNumber}.cluster{clusterIndex}.t.mean = mean(t);
    clusters{channelNumber}.cluster{clusterIndex}.t.median = median(t);
    clusters{channelNumber}.cluster{clusterIndex}.t.std = std(t);
    clusters{channelNumber}.cluster{clusterIndex}.f.min = min(f);
    clusters{channelNumber}.cluster{clusterIndex}.f.max = max(f);
    clusters{channelNumber}.cluster{clusterIndex}.f.middle = (min(f) + max(f)) / 2;
    clusters{channelNumber}.cluster{clusterIndex}.f.bandwidth = max(f) - min(f);
    clusters{channelNumber}.cluster{clusterIndex}.f.mean = mean(f);
    clusters{channelNumber}.cluster{clusterIndex}.f.median = median(f);
    clusters{channelNumber}.cluster{clusterIndex}.f.std = std(f);
    clusters{channelNumber}.cluster{clusterIndex}.z.min = min(z);
    clusters{channelNumber}.cluster{clusterIndex}.z.max = max(z);
    clusters{channelNumber}.cluster{clusterIndex}.z.middle = (min(z) + max(z)) / 2;
    clusters{channelNumber}.cluster{clusterIndex}.z.range = max(z) - min(z);
    clusters{channelNumber}.cluster{clusterIndex}.z.mean = mean(z);
    clusters{channelNumber}.cluster{clusterIndex}.z.median = median(z);
    clusters{channelNumber}.cluster{clusterIndex}.z.std = std(z);
    
end

% disp('Press any key for next iteration');
% pause;
% Eavaluate criteria to check if to run this module again
% If yes, call this module within itself from here,
% else, produce report.

% Prepare arrays for reporting
n = []; tmst = []; tstart = []; tend = []; tmiddle = []; tduration = [];
tmedian = []; tmean = []; tstd = []; fmst = []; fmin = []; fmax = [];
fmiddle = []; fbandwidth = []; fmedian = []; fmean = []; fstd = []; 
ztotal = []; zmst = []; zmin = []; zmax = []; zmiddle = []; zrange = [];
zmedian = []; zmean = []; zstd = [];

for clusterIndex = 1:clusters{channelNumber}.numberOfClusters,
    n = [n clusters{channelNumber}.cluster{clusterIndex}.n];
    ztotal = [ztotal clusters{channelNumber}.cluster{clusterIndex}.z.total];
    tmst = [tmst clusters{channelNumber}.cluster{clusterIndex}.t.mst];
    fmst = [fmst clusters{channelNumber}.cluster{clusterIndex}.f.mst];
    zmst = [zmst clusters{channelNumber}.cluster{clusterIndex}.z.mst];
    tstart = [tstart clusters{channelNumber}.cluster{clusterIndex}.t.start];
    tend = [tend clusters{channelNumber}.cluster{clusterIndex}.t.end];
    tmiddle = [tmiddle clusters{channelNumber}.cluster{clusterIndex}.t.middle];
    tduration = [tduration clusters{channelNumber}.cluster{clusterIndex}.t.duration];
    tmean = [tmean clusters{channelNumber}.cluster{clusterIndex}.t.mean];
    tmedian = [tmedian clusters{channelNumber}.cluster{clusterIndex}.t.median];
    tstd = [tstd clusters{channelNumber}.cluster{clusterIndex}.t.std];
    fmin = [fmin clusters{channelNumber}.cluster{clusterIndex}.f.min];
    fmax = [fmax clusters{channelNumber}.cluster{clusterIndex}.f.max];
    fmiddle = [fmiddle clusters{channelNumber}.cluster{clusterIndex}.f.middle];
    fbandwidth = [fbandwidth clusters{channelNumber}.cluster{clusterIndex}.f.bandwidth];
    fmean = [fmean clusters{channelNumber}.cluster{clusterIndex}.f.mean];
    fmedian = [fmedian clusters{channelNumber}.cluster{clusterIndex}.f.median];
    fstd = [fstd clusters{channelNumber}.cluster{clusterIndex}.f.std];
    zmin = [zmin clusters{channelNumber}.cluster{clusterIndex}.z.min];
    zmax = [zmax clusters{channelNumber}.cluster{clusterIndex}.z.max];
    zmiddle = [zmiddle clusters{channelNumber}.cluster{clusterIndex}.z.middle];
    zrange = [zrange clusters{channelNumber}.cluster{clusterIndex}.z.range];
    zmean = [zmean clusters{channelNumber}.cluster{clusterIndex}.z.mean];
    zmedian = [zmedian clusters{channelNumber}.cluster{clusterIndex}.z.median];
    zstd = [zstd clusters{channelNumber}.cluster{clusterIndex}.z.std];
end

% print results
disp('printing results...');
[maximumClusterEnergy, maximumClusterIndex] = max(ztotal);
maximumClusterTileEnergy = zmst(find(ztotal == max(ztotal), 1 , 'first'));
maximumTileEnergy = max(events{channelNumber}.normalizedEnergy);
fprintf(1, '\n');
fprintf(1, '          number of clusters:  %d\n', ...
        clusters{channelNumber}.numberOfClusters);
fprintf(1, '         maximum cluster snr:  %.1f\n', ...
        sqrt(2 * maximumClusterEnergy));
fprintf(1, '    maximum cluster tile snr:  %.1f\n', ...
        sqrt(2 * maximumClusterTileEnergy));
fprintf(1, '            maximum tile snr:  %.1f\n', ...
        sqrt(2 * maximumTileEnergy));
fprintf(1, '          improvement factor:  %.1f\n', ...
        sqrt(maximumClusterEnergy / maximumClusterTileEnergy));
fprintf(1, '\n');

% display clustering
disp('plotting clusters...');
clf;
set(gca, 'FontSize', 16);
axis([-4 +4 6 10]);
set(gca, 'YTick', 6:10);
set(gca, 'YTickLabel', 2.^(6:10));
xlabel('Time [seconds]');
ylabel('Frequency [Hz]');
title(strrep(channelName, '_', '\_'));
hold on;
markers = {'ro', 'go', 'bo', 'co', 'mo', 'ko', ...
     'rs', 'gs', 'bs', 'cs', 'ms', 'ks', ...
     'rd', 'gd', 'bd', 'cd', 'md', 'kd', ...
     'rv', 'gv', 'bv', 'cv', 'mv', 'kv', ...
     'r^', 'g^', 'b^', 'c^', 'm^', 'k^', ...
     'r<', 'g<', 'b<', 'c<', 'm<', 'k<', ...
     'r>', 'g>', 'b>', 'c>', 'm>', 'k>', ...
     'rp', 'gp', 'bp', 'cp', 'mp', 'kp', ...
     'rh', 'gh', 'bh', 'ch', 'mh', 'kh'};
markers = [markers markers markers markers];
indices = find(tempClusters == maximumClusterIndex);
handle = plot(tiles(indices, 1) - centerTime, ...
        log2(tiles(indices, 2)), ...
        'ko');
set(handle, 'Color', 0.9 * [1 1 1]);
set(handle, 'MarkerSize', 30);
set(handle, 'MarkerFaceColor', get(handle, 'Color'));
for clusterIndex = 1 : length(uniqueClusters),
    indices = find(tempClusters == uniqueClusters(clusterIndex));
    handle = plot(tiles(indices, 1) - centerTime, ...
              log2(tiles(indices, 2)), ...
              markers{clusterIndex});
    set(handle, 'MarkerSize', 10);
    set(handle, 'MarkerFaceColor', get(handle, 'Color'));
end

% produce ASCII report
disp('producing ASCII report');     

n = n'; tmst = tmst'; tstart = tstart'; tend = tend'; tmiddle = tmiddle';
tduration = tduration'; tmedian = tmedian'; tmean = tmean'; tstd = tstd';
fmst = fmst'; fmin = fmin'; fmax = fmax'; fmiddle = fmiddle';
fbandwidth = fbandwidth'; fmedian = fmedian'; fmean = fmean'; fstd = fstd';
ztotal = ztotal'; zmst = zmst'; zmin = zmin'; zmax = zmax'; zmiddle = zmiddle';
zrange = zrange'; zmedian = zmedian'; zmean = zmean'; zstd = zstd';

filename = sprintf('%s', ['Clusters_channel_' strcat(channelName) '_Number_' num2str(channelNumber) '.txt']);
filematrix  = [n ztotal tmst fmst zmst tstart tend tmiddle tduration tmedian tmean tstd, fmin ... 
    fmax fmiddle fbandwidth fmedian fmean fstd zmin zmax zmiddle zrange zmedian zmean zstd];
save(filename, 'filematrix', '-ASCII', '-DOUBLE');


return;
