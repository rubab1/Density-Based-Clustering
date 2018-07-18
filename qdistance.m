function distances = qdistance(significants, durationInflation, ...
                               bandwidthInflation);
% QDISTANCE Compute distances between significant Q transform tiles
%
% QDISTANCE computes the distance between significant Q transform
% tiles produced by QTHREHSOLD or QSELECT.
%
% distances = qdistance(significants, durationInflation, bandwidthInflation);
%
%   distances            cell array of significant tiles distances
%
%   significants         cell array of significant tiles properties
%   durationInflation    multiplicative scale factor for duration
%   bandwidthInflation   multiplicative scale factor for bandwidth
%
% The distance is defined to be the dimensionless euclidean distance
% between tiles in the time-frequency plane after normalizing by the
% mean duration and bandwidth of each tile pair.
%
% QDISTANCE expects a cell array of Q transform event structures with
% one cell per channel.  The event structures must contain at least
% the following fields, which describe the properties of statistically
% significant tiles used to compute distance.
% 
%   time                 center time of tile [gps seconds]
%   frequency            center frequency of tile [Hz]
%   q                    quality factor of tile []
%   normalizedEnergy     normalized energy of tile []
%
% QDISTANCE returns a cell array of Q transform distance structures
% with one cell per cahnnel.  In addition to a structure identifier,
% the distance structures contains the following single field.
%
%   distance             pairwise distance between tiles
%
% Distances are returned in the same format as the PDIST function.
% In particular, distances are reported as row vectors of length
% N * (N - 1) / 2, where N is the number of significant tiles for
% a given channel.  This row vector is arranged in the order of
% (1,2), (1,3), ..., (1,N), (2,3), ..., (2,N), ..., (N-1, N).  Use
% the SQUAREFORM function to convert distances into a matrix format.
%
% The optional durationInflation and bandwidthInflation arguments are
% multiplicative scale factors that are applied to the duration and
% bandwidth of significant tiles prior to determining their distance.
% If not specified, these parameters both default to unity such that
% the resulting tiles have unity time-frequency area.
%
% See also QTHRESHOLD, QSELECT, QCLUSTER, SQUAREFORM, and PDIST.

% Rubab Khan
% rmk2109@columbia.edu
%
% Shourov Chatterji
% shourov@ligo.caltech.edu
%
% 2006-Jul-07

% $Id:$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 3, nargin));

% default tile inflation factors
if nargin < 2,
  durationInflation = 1.0;
end
if nargin < 3,
  bandwidthInflation = 1.0;
end

% if input events are not in a cell array,
if ~iscell(significants),
  
  % insert significant events into a single cell
  significants = mat2cell(significants, size(significants, 1), ...
                          size(significants, 2));

% otherwise, continue
end

% force one dimensional cell array
significants = significants(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(significants);

% validate significant event structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(significants{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       initialize distances structures                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of significant tile distances
distances = cell(size(significants));

% begin loop over channels
for channelNumber = 1 : numberOfChannels

  % insert structure identification string
  distances{channelNumber}.id = 'Discrete Q-transform distance structure';
  
% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      create pairwise list of tiles                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % number of significant tiles
  numberOfSignificants = length(significants{channelNumber}.time);

  % number of unique significant tile pairs
  numberOfPairs = numberOfSignificants * (numberOfSignificants - 1) / 2;

  % build list of pairwise indices
  [pairIndices1, pairIndices2] = meshgrid(1 : numberOfSignificants);
  pairIndices1 = tril(pairIndices1, -1);
  pairIndices2 = tril(pairIndices2, -1);
  pairIndices1 = pairIndices1(:);
  pairIndices2 = pairIndices2(:);
  pairIndices1 = pairIndices1(find(pairIndices1));
  pairIndices2 = pairIndices2(find(pairIndices2));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                        determine tile properties                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % extract significant tile properties
  time1 = significants{channelNumber}.time(pairIndices1);
  time2 = significants{channelNumber}.time(pairIndices2);
  frequency1 = significants{channelNumber}.frequency(pairIndices1);
  frequency2 = significants{channelNumber}.frequency(pairIndices2);
  q1 = significants{channelNumber}.q(pairIndices1);
  q2 = significants{channelNumber}.q(pairIndices2);
  normalizedEnergy1 = significants{channelNumber}.normalizedEnergy(pairIndices1);
  normalizedEnergy2 = significants{channelNumber}.normalizedEnergy(pairIndices2);

  % determine significant tile dimensions
  bandwidth1 = 2 * sqrt(pi) * frequency1 ./ q1;
  bandwidth2 = 2 * sqrt(pi) * frequency2 ./ q2;
  duration1 = 1 ./ bandwidth1;
  duration2 = 1 ./ bandwidth2;

  % apply tile inflation factors
  duration1 = duration1 * durationInflation;
  duration2 = duration2 * durationInflation;
  bandwidth1 = bandwidth1 * bandwidthInflation;
  bandwidth2 = bandwidth2 * bandwidthInflation;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                         determine tile distances                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % determine time and frequency distance scales
  timeScale = (duration1 .* normalizedEnergy1 + ...
               duration2 .* normalizedEnergy2) / ...
              (normalizedEnergy1 + normalizedEnergy2);
  frequencyScale = (bandwidth1 .* normalizedEnergy1 + ...
                    bandwidth2 .* normalizedEnergy2) / ...
                   (normalizedEnergy1 + normalizedEnergy2);

  % compute normalized time and frequency distance between tiles
  timeDistance = abs(time2 - time1) / timeScale;
  frequencyDistance = abs(frequency2 - frequency1) / frequencyScale;

  % determine normalized euclidean distance between tiles
  distances{channelNumber}.distance = ...
      sqrt(timeDistance.^2 + 30 * frequencyDistance.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over channels                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   return statistically significant events                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
