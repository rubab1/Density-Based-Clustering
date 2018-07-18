function clusters = qcluster2(mosaics, distances, clusteringRadius, ...
                              numberThreshold)
% QCLUSTER2 Density based clustering algorithm
%
% QCLUSTER2 implements a density based clustering algorithm similar the
% algorithm DBSCAN.
%
% usage:  clusters = qcluster2(mosaics, distances, clusteringRadius, ...
%                              numberThreshold);
%
%   mosaics            .
%   distances          .
%   clusteringRadius   .
%   numberThreshold    .
%
%   clusters           .
%
% If not specified, a clusteringRadius of 10 and a numberThreshold of 4 are
% assumed.
%
% See also .

% Rubab Khan
% rmk2109@columbia.edu
%
% Shourov Chatterji
% shourov@ligo.caltech.edu
%
% 2006-Jul-13

% $Id:$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 4, nargin));

if nargin < 3,
  clusteringRadius = 6;
end
if nargin < 4,
  numberThreshold = 4;
end

% if input distances are not in a cell array,
if ~iscell(distances),

  % insert distances into a single cell
  distances = mat2cell(distances, size(distances, 1), size(distances, 2));

% otherwise, continue
end

% force one dimensional cell array
distances = distances(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(distances);

% validate significant event structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(mosaics{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
end

% validate distance structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(distances{channelNumber}.id, ...
             'Discrete Q-transform distance structure'),
    error('input argument is not a discrete Q transform distance structure');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        initialize results structures                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of significant tile distances
clusters = cell(size(distances));

% begin loop over channels
for channelNumber = 1 : numberOfChannels

  % insert structure identification string
  clusters{channelNumber}.id = 'Discrete Q-transform cluster structure';

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                     intiialize clustering algorithm                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % number of tiles
  numberOfTiles = length(mosaics{channelNumber}.time);

  % convert distance structure to matrix
  distanceMatrix = squareform(distances{channelNumber}.distance);

  % initialize the tiles structure
  tiles = cell(1, numberOfTiles);

  tiles{1}.maximumNumberOfRecursions = 500;
  tiles{1}.numberThreshold = numberThreshold;

  % begin loop over tiles
  for tileNumber = 1 : numberOfTiles,

    % find indices of neighboring tiles
    tiles{tileNumber}.neighborTileNumbers = ...
        find((distanceMatrix(:, tileNumber) <= clusteringRadius) & ...
             (distanceMatrix(:, tileNumber) ~= 0));

    % number of neighboring tiles
    tiles{tileNumber}.numberOfNeighbors = ...
        length(tiles{tileNumber}.neighborTileNumbers);

    % if neighboring tiles exceed critical density,
    if tiles{tileNumber}.numberOfNeighbors >= numberThreshold,

      % identify tile as potential seed
      tiles{tileNumber}.clusterNumber = 0;

      % sort neighboring tiles by increasing normalized energy
      [ignore, sortedNeighborIndices] = ...
          sort(mosaics{channelNumber}.normalizedEnergy( ...
              tiles{tileNumber}.neighborTileNumbers));
      clear ignore;
      sortedNeighborIndices = sortedNeighborIndices(end : -1 : 1);
      tiles{tileNumber}.neighborTileNumbers = ...
          tiles{tileNumber}.neighborTileNumbers(sortedNeighborIndices);

    % otherwise,
    else

      % identify as not a potential seed
      tiles{tileNumber}.clusterNumber = NaN;

    % continue
    end

  % end loop over tiles
  end

  % free distance matrix memory
  clear distanceMatrix;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                        cluster significant tiles                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % sort significant tiles by increasing normalized energy
  [ignore, sortedTileIndices] = sort(mosaics{channelNumber}.normalizedEnergy);
  clear ignore;
  sortedTileIndices = sortedTileIndices(end : -1 : 1);

  % initialize number of clusters
  clusterNumber = 0;

  % begin loop over sorted tiles
  for sortedTileNumber = 1 : numberOfTiles,

    % find tile number corersponding to sorted tile number
    tileNumber = sortedTileIndices(sortedTileNumber);

    % if current tile has not been processed,
    if (tiles{tileNumber}.clusterNumber == 0),

      % create a new cluster
      clusterNumber  = clusterNumber + 1;

      % assign current tile to new cluster
      tiles{tileNumber}.clusterNumber = clusterNumber;

      % set recursion limit
      set(0, 'RecursionLimit', tiles{1}.maximumNumberOfRecursions);

      % find other tiles in the cluster
      [numberOfRecursions, tiles] = qcluster2b(tileNumber, tiles, 3);

      % reset recursion limit
      set(0, 'RecursionLimit', 100);

      % report status
      % fprintf(1, 'QCLUSTER2: cluster %d required %d recursions\n', ...
      %             clusterNumber, numberOfRecursions);

    % otherwise, skip to the next tile
    end

  % end loop over sorted tiles
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                             collect results                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % collect cluster numbers
  for tileNumber = 1 : numberOfTiles
    clusters{channelNumber}.cluster(tileNumber) = ...
        tiles{tileNumber}.clusterNumber;
  end

  % identify unique cluster numbers
  clusteredTileIndices = find(~isnan(clusters{channelNumber}.cluster));
  uniqueClusters = ...
      unique(clusters{channelNumber}.cluster(clusteredTileIndices));

  % determine number of clusters
  clusters{channelNumber}.numberOfClusters = length(uniqueClusters);

  % compress cluster numbers
  for clusterNumber = 1 : clusters{channelNumber}.numberOfClusters,
    tilesInCluster = find(clusters{channelNumber}.cluster == ...
                          uniqueClusters(clusterNumber));
    clusters{channelNumber}.cluster(tilesInCluster) = clusterNumber;
  end

  % report results
  % fprintf(1, '----------------------------------\n');
  % fprintf(1, 'QCLUSTER2: %#03d clusters\n', ...
  %         clusters{channelNumber}.numberOfClusters);
  % for clusterNumber = 1 : clusters{channelNumber}.numberOfClusters,
  %   tilesInCluster = ...
  %       find(clusters{channelNumber}.cluster == clusterNumber);
  %   clusterEnergy = ...
  %       sum(mosaics{channelNumber}.normalizedEnergy(tilesInCluster));
  %   fprintf(1, 'QCLUSTER2: cluster %#03d: %#03d tiles, %#09.3e energy\n', ...
  %           clusterNumber, length(tilesInCluster), clusterEnergy);
  % end
  % fprintf(1, '----------------------------------\n');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               return clusters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
