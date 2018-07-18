function [numberOfRecursions, tiles] = qcluster2b(tileNumber, tiles, numberOfRecursions)
% QCLUSTER2B Helper function for density based clustering algorithm
%
% QCLUSTER2B is an auxiliary function for the QCLUSTER3 density based
% clustering algorithm.  It implements recursive clustering of tiles
% and should only be called by QCLUSTER3.
%
% usage: [numberOfRecursions, tiles] = qcluster2b(tileNumber, tiles, numberOfRecursions);
%
%   tileNumber          .
%   tiles               .
%   numberOfRecursions  .
%
%   tiles               .
%   numberOfRecursions  .
%
%
%
% See also QCLUSTER3.

% Rubab Khan
% rmk2109@columbia.edu
%
% Shourov K. Chatterji
% shourov@ligo.caltech.edu
%
% 2006-Jul-13

% $Id:$

% apply default arguments
if nargin < 3
  numberOfRecursions = 3;
end

% if recursion limit is exceeded
if (numberOfRecursions >= tiles{1}.maximumNumberOfRecursions),

  % reset current tile
  tiles{tileNumber}.clusterNumber = 0;
  
  % return to calling function
  return;

% otherwise continue
end

% increment recursion number counter
numberOfRecursions = numberOfRecursions + 1;

% determine number of significant tiles
numberOfTiles = length(tiles);

% begin loop over neighboring tiles
for neighborNumber = 1 : tiles{tileNumber}.numberOfNeighbors,

  % determine tile number for neighbor tile
  neighborTileNumber = tiles{tileNumber}.neighborTileNumbers(neighborNumber);

  % determine current cluster number
  currentClusterNumber = tiles{tileNumber}.clusterNumber;

  % determine cluster number for neighbor tile
  oldClusterNumber = tiles{neighborTileNumber}.clusterNumber;

  % if neighbor tile has not been processed,
  if oldClusterNumber == 0

    % assign neighbor tile to current cluster
    tiles{neighborTileNumber}.clusterNumber = currentClusterNumber;

    % continue to build cluster
    [numberOfRecursions, tiles] = ...
        qcluster2b(neighborTileNumber, tiles, numberOfRecursions);

  % if neighbor tile is already in another cluster
  elseif ((oldClusterNumber > 0) && ...
          (oldClusterNumber ~= currentClusterNumber)),

      % merge current cluster into old cluster
      for testTileNumber = 1 : numberOfTiles,
        if (tiles{testTileNumber}.clusterNumber == currentClusterNumber),
          tiles{testTileNumber}.clusterNumber = oldClusterNumber;
        end
      end

      % report status
      if (tiles{neighborTileNumber}.numberOfNeighbors < ...
          tiles{1}.numberThreshold),
        mergeType = 'border tile';
      else
        mergeType = 'maximum recursion tile';
      end
      % fprintf(1, 'QCLUSTER2B: cluster %d merged with cluster %d via %s\n', ...
      %         currentClusterNumber, oldClusterNumber, mergeType);

  % otherwise,
  else

    % assign neighbor tile as border tile of current cluster
    tiles{neighborTileNumber}.clusterNumber = currentClusterNumber;

  % continue
  end

% end loop over neighboring tiles
end

% return to calling function
return;
