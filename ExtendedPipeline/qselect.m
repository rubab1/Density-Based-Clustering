function events = qselect(significants, durationInflation, ...
                          bandwidthInflation, maximumEvents, ...
                          channelNames, debugLevel);
% QSELECT Identify statistically significant events in Discrete Q transforms
%
% QSELECT selects statistically significant events from the set of statitically
% significant Q transform tiles.  Events are defined by the properties of their
% most significant tile and are identified by exluding the less significant of
% any overlapping tiles.  The input significant tiles are first sorted by
% decreasing normalized energy.  Starting with the most significant tile, tiles
% are discarded if they overlap with a more significant tile.  The remaining set
% of tiles comprises a minimal set of tiles that describes an event.
%
% QSELECT returns a cell array of significant event properties
%
% usage: events = qselect(significants, durationInflation, ...
%                         bandwidthInflation, maximumEvents, ...
%                         channelNames, debugLevel);
%
%   significants         cell array of significant tiles properties
%   durationInflation    multiplicative scale factor for duration
%   bandwidthInflation   multiplicative scale factor for bandwidth
%   maximumEvents        maximum allowable number of events
%   channelNames         cell array of channel names
%   debugLevel           verboseness of debug output
%
%   events               cell array of significant event properties
%
% The optional durationInflation and bandwidthInflation arguments are
% multiplicative scale factors that are applied to the duration and bandwidth of
% significant tiles prior to testing for overlap.  If not specified, these
% parameters both default to unity such that the resulting tiles have unity
% time-frequency area.  The normalized energy of the resulting tiles are scaled
% by the product of the duration and bandwidth inflation factors to avoid over
% counting the total energy of clusters of tiles.  Likewise, the amplitude of
% the resulting tiles is scaled by the square root of the product of the
% duration and bandwidth inflation factors.
%
% The optional maximumEvents argument provides a safety mechanism to limit the
% total number of events returned by QSELECT.  If this maximum number of events
% is exceeded, an overflow flag is set, only the maximumEvents most significant
% events are returned, and a warning is issued if debugLevel is set to unity or
% higher.  By default, maximumEvents is set to infinity and debugLevel is set to
% unity.
%
% QSELECT both expects and returns a cell array of Q transform event structures
% with one cell per channel.  The event structures contain the following
% required fields, which describe the properties of statistically significant
% tiles.  Additional fields such as amplitude, phase, or coherent transform
% properties are optional and are retained along with the required fields.
%
%   time                 center time of tile [gps seconds]
%   frequency            center frequency of tile [Hz]
%   duration             duration of tile [seconds]
%   bandwidth            bandwidth of tile [Hz]
%   normalizedEnergy     normalized energy of tile []
%
% The event structures also contain the following flag which indicated if the
% maximum number of significant tiles of significant events was exceeded.
%
%   overflowFlag         boolean overflow flag
%
% See also QTILE, QCONDITION, QTRANSFORM, QTHRESHOLD, QEXAMPLE, and QPIPELINE.

% Shourov K. Chatterji
% shourov@ligo.caltech.edu

% Leo C. Stein
% lstein@ligo.mit.edu

% $Id: qselect.m,v 1.9 2007/05/07 16:41:54 shourov Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(1, 6, nargin));

% apply default arguments
if (nargin < 2) || isempty(durationInflation),
  durationInflation = 1.0;
end
if (nargin < 3) || isempty(bandwidthInflation),
  bandwidthInflation = 1.0;
end
if (nargin < 4) || isempty(maximumEvents),
  maximumEvents = Inf;
end
if (nargin < 5) || isempty(channelNames),
  channelNames = [];
end
if (nargin < 6) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
if ~iscell(significants),
  significants = mat2cell(significants, size(significants, 1), ...
                          size(significants, 2));
end
if ~isempty(channelNames) & ~iscell(channelNames),
  channelNames = mat2cell(channelNames, size(channelNames, 1), ...
                          size(channelNames, 2));
end

% force one dimensional cell arrays
significants = significants(:);
channelNames = channelNames(:);

% determine number of channels
numberOfChannels = length(significants);

% provide default channel names
if isempty(channelNames),
  channelNames = cell(numberOfChannels, 1);
  for channelNumber = 1 : numberOfChannels,
    channelNames{channelNumber} = ['X' int2str(channelNumber) ':STRAIN'];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate significant event structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(significants{channelNumber}.id, ...
             'Discrete Q-transform event structure'),
    error('input argument is not a discrete Q transform event structure');
  end
end

% validate number of channels
if length(channelNames) ~= numberOfChannels,
  error('number of channelNames is inconsistent with number of transforms');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            initialize statistically significant events structures            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty array of significant event indices
eventIndices = cell(size(significants));

% create empty cell array of significant event structures
events = cell(size(significants));

% begin loop over channels
for channelNumber = 1 : numberOfChannels

  % insert structure identification string
  events{channelNumber}.id = 'Discrete Q-transform event structure';
  
  % propogate overflow flag
  events{channelNumber}.overflowFlag = ...
      significants{channelNumber}.overflowFlag;

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over channels                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                   sort by decreasing normalized energy                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % sort tile indices by normalized energy
  [ignore, sortedIndices] = ...
      sort(significants{channelNumber}.normalizedEnergy);

  % sort by decreasing normalized energy
  sortedIndices = fliplr(sortedIndices);

  % get list of field names available for this channel
  fieldNames = fieldnames(significants{channelNumber});

  % remove id and overflowFlag fields from list of field names
  fieldNames = fieldNames(~strcmp(fieldNames, 'id') & ...
                          ~strcmp(fieldNames, 'overflowFlag'));

  % number of field names available for this channel
  numberOfFields = length(fieldNames);

  % reorder significant tile properties by decreasing normalized energy
  for fieldNumber = 1 : numberOfFields,
    significants{channelNumber}.(fieldNames{fieldNumber}) = ...
      significants{channelNumber}.(fieldNames{fieldNumber})(sortedIndices);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                           find tile boundaries                             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % vector of significant tile bandwidths
  % significants{channelNumber}.bandwidth = 2 * sqrt(pi) * ...
  %     significants{channelNumber}.frequency ./ significants{channelNumber}.q;

  % vector of significant tile durations
  % significants{channelNumber}.duration = ...
  %     1 ./ significants{channelNumber}.bandwidth;

  % vector of significant tile start times
  minimumTimes = significants{channelNumber}.time - ...
      durationInflation * ...
      significants{channelNumber}.duration / 2;

  % vector of significant tile stop times
  maximumTimes = significants{channelNumber}.time + ...
      durationInflation * ...
      significants{channelNumber}.duration / 2;

  % vector of significant tile lower frequencies
  minimumFrequencies = significants{channelNumber}.frequency - ...
      bandwidthInflation * ...
      significants{channelNumber}.bandwidth / 2;

  % vector of significant tile upper frequencies
  maximumFrequencies = significants{channelNumber}.frequency + ...
      bandwidthInflation * ...
      significants{channelNumber}.bandwidth / 2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      compress significant tile list                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % number of significant tiles in list
  numberOfTiles = length(significants{channelNumber}.time);

  % if input significant tile list is empty,
  if numberOfTiles == 0,

    % set empty event list
    for fieldNumber = 1 : numberOfFields,
      events{channelNumber}.(fieldNames{fieldNumber}) = [];
    end

    % skip to next channel
    continue;

  % otherwise, continue
  end

  % initialize event list
  eventIndices{channelNumber} = 1;

  % begin loop over significant tiles
  for tileIndex = 2 : numberOfTiles,

    % determine if current tile overlaps any events
    overlap = any((minimumTimes(tileIndex) < ...
                   maximumTimes(eventIndices{channelNumber})) & ...
                  (maximumTimes(tileIndex) > ...
                   minimumTimes(eventIndices{channelNumber})) & ...
                  (minimumFrequencies(tileIndex) < ...
                   maximumFrequencies(eventIndices{channelNumber})) & ...
                  (maximumFrequencies(tileIndex) > ...
                   minimumFrequencies(eventIndices{channelNumber})));

    % if tile does not overlap with any event,
    if ~overlap,

      % append it to the list of events
      eventIndices{channelNumber} = [eventIndices{channelNumber} tileIndex];

    % otherwise, continue
    end

  % end loop over significant tiles
  end

  % extract events from significant tiles
  for fieldNumber = 1 : numberOfFields,
    events{channelNumber}.(fieldNames{fieldNumber}) = ...
      significants{channelNumber}.(fieldNames{fieldNumber})(eventIndices{channelNumber});
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                    check for excessive number of events                    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % determine number of significant tiles in channel
  numberOfEvents = length(events{channelNumber}.time);

  % if maximum allowable number of significant tiles is exceeded
  if numberOfEvents > maximumEvents,

    % issue warning
    if debugLevel > 0,
      fprintf(1, 'WARNING: %s: maximum number of events exceeded.\n', ...
              channelNames{channelNumber});
    end

    % set overflow flag
    events{channelNumber}.overflowFlag = 1;

    % indices of most significant tiles
    maximumIndices = 1 : maximumEvents;

    % truncate lists of significant event properties
    for fieldNumber = 1 : numberOfFields,
      events{channelNumber}.(fieldNames{fieldNumber}) = ...
        events{channelNumber}.(fieldNames{fieldNumber})(maximumIndices);
    end

  % otherwise continue
  end

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
