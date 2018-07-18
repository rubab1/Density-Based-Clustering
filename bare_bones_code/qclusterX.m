function clusters = qclusterX(events, durationInflation, bandwidthInflation, channelNumber, channelName, centerTime)

warning off all;

% measure distances
disp('measuring tile distances...');
distances = qdistance(events, durationInflation, bandwidthInflation);

distMat = squareform(distances{channelNumber}.distance);

% measure distances
% disp('measuring tile distances...');
% tiles = [events{channelNumber}.time' events{channelNumber}.frequency' ...
%    events{channelNumber}.q' events{channelNumber}.normalizedEnergy'];
% distances = pdist(tiles, @qmetric3);

clf;
sortMat = sort(distMat, 'ascend');
plot(sort((sortMat(5,:)), 'descend'));
xlabel('Tile Number');
ylabel('Distance to fourth-distant-neighbor');
title('4-distance graph');
disp('getting 4-distance, press any key');
pause;

clear sortMat;

% Initialize the point structure
point = cell(1);
% Neigborhood Radius
NR = 8;
% Neighbor Number
NN = 4;
point{1}.NN = NN;
% Number of Tiles
N = length(distMat(:,1));
point{1}.N = N;
% Matrix with tile properties
tiles = [events{channelNumber}.time' events{channelNumber}.frequency' ...
    events{channelNumber}.q' events{channelNumber}.normalizedEnergy'];

for i = 1:N
    point{i}.neighbor.I = find(distMat(:,i)<=NR&distMat(:,i)~=0); % Get index of neighbors, but not of self
    point{i}.neighbor.n = length(point{i}.neighbor.I); % Get number of neighbors
    if point{i}.neighbor.n >= NN % Can it be a seed
        point{i}.ID = 0;
        [zz, Index] = sort(tiles(point{i}.neighbor.I, 4), 'descend');
        point{i}.neighbor.I = point{i}.neighbor.I(Index);
    else
        point{i}.ID = NaN;
    end
end

clear distMat;

% Index of tiles sorted along Normalized Energy
[zz,Index] = sort(tiles(:,4), 'descend');

% Set recursion limit
point{1}.rLimit = 500;
set(0,'RecursionLimit', point{1}.rLimit);

% Initiate Cluster ID
clusterID = 0;
for i = 1:N
    ii = Index(i); % Pick higher significance seeds first
    if (point{ii}.ID == 0) % If point is possible seed not in other cluster already
        clusterID  = clusterID + 1;
        point{ii}.ID = clusterID;
        [k, point] = expandCluster(ii,point,3);
        K = sprintf('Number of recursion %d for cluster %d', k, clusterID);
        disp(K);
    end
end

if 0,
for i = 1:N
    if ~(point{i}.ID > 0)
        clusterID  = clusterID + 1;
        point{i}.ID = clusterID;
    end
end
end

% Reset recursion limit
set(0,'RecursionLimit', 100);

%%%%%%%%%%%%% Clustering ends and reporting starts here %%%%%%%%%%%%%%%


% sort clusters
disp('sorting clusters...');
tempClusters = 0;
for i = 1:N
    tempClusters(i) = point{i}.ID;
end
uniqueClusters = sort(unique(tempClusters(tempClusters>0)));

% Initialize the clusters structure
clusters{channelNumber} = cell(1);

clusters{channelNumber}.id = 'Clusters of Unique Significant Tiles';
clusters{channelNumber}.numberOfClusters = length(uniqueClusters);

if clusters{channelNumber}.numberOfClusters > 0
    for i = clusters{channelNumber}.numberOfClusters
        tempClusters(tempClusters == uniqueClusters(i)) = i;
    end
    % Read cluster info
    for clusterIndex = 1:clusters{channelNumber}.numberOfClusters

        % Read time, frequency,and normalized energy info of each tile
        i = find(tempClusters == clusterIndex);
        t = events{channelNumber}.time(i);
        f = events{channelNumber}.frequency(i);
        z = events{channelNumber}.normalizedEnergy(i);

        % Tile indices, total number of tiles and sum of tile-energy
        clusters{channelNumber}.cluster{clusterIndex}.indices = i;
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
    % Eavaluate criteria to check if to cluster the clusters
    % if not, produce report.

    % Prepare arrays for reporting
    n = []; tmst = []; tstart = []; tend = []; tmiddle = []; tduration = [];
    tmedian = []; tmean = []; tstd = []; fmst = []; fmin = []; fmax = [];
    fmiddle = []; fbandwidth = []; fmedian = []; fmean = []; fstd = []; 
    ztotal = []; zmst = []; zmin = []; zmax = []; zmiddle = []; zrange = [];
    zmedian = []; zmean = []; zstd = [];

    for clusterIndex = 1:clusters{channelNumber}.numberOfClusters,
        n(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.n;
        ztotal(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.total;
        tmst(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.mst;
        fmst(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.mst;
        zmst(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.mst;
        tstart(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.start;
        tend(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.end;
        tmiddle(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.middle;
        tduration(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.duration;
        tmean(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.mean;
        tmedian(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.median;
        tstd(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.t.std;
        fmin(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.min;
        fmax(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.max;
        fmiddle(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.middle;
        fbandwidth(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.bandwidth;
        fmean(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.mean;
        fmedian(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.median;
        fstd(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.f.std;
        zmin(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.min;
        zmax(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.max;
        zmiddle(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.middle;
        zrange(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.range;
        zmean(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.mean;
        zmedian(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.median;
        zstd(clusterIndex) = clusters{channelNumber}.cluster{clusterIndex}.z.std;
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
    
end

return;
