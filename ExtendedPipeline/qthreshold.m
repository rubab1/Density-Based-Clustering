function [significants, channelNames] = ...
    qthreshold(transforms, tiling, startTime, falseEventRate, ...
               referenceTime, timeRange, frequencyRange, qRange, ...
               maximumSignificants, analysisMode, falseVetoRate, ...
               uncertaintyFactor, correlationFactor, applyVeto, ...
               channelNames, debugLevel)
% QTHRESHOLD Identify statistically significant tiles in Discrete Q transforms
%
% QTHRESHOLD identifies discrete Q transform coefficients whose magnitudes
% exceed the threshold that approximately yields the specified single channel
% false rate assuming ideal white noise.
%
% usage: significants = qthreshold(transforms, tiling, startTime, ...
%                                  falseEventRate, referenceTime, ...
%                                  timeRange, frequencyRange, qRange, ...
%                                  maximumSignificants, analysisMode, ...
%                                  falseVetoRate, uncertaintyFactor, ...
%                                  correlationFactor, applyVeto, ...
%                                  channelNames, debugLevel);
%
%   transforms           cell array of input Q transform structures
%   tiling               discrete Q transform tiling structure from QTILE
%   startTime            GPS start time of Q transformed data
%   falseEventRate       desired white noise false event rate [Hz]
%   referenceTime        reference time for time range to threshold on
%   timeRange            vector range of relative times to threshold on
%   frequencyRange       vector range of frequencies to threshold on
%   qRange               scalar Q or vector range of Qs to threshold on
%   maximumSignificants  maximum allowable number of significant tiles
%   analysisMode         string name of analysis mode to implement
%   falseVetoRate        desired white noise veto rate [Hz]
%   uncertaintyFactor    squared calibration uncertainty factor
%   correlationFactor    ***** TO BE DEFINED *****
%   applyVeto            boolean flag to apply same time veto
%   channelNames         cell array of thresholded channel names
%   debugLevel           verboseness of debug output
%
%   significants         cell array of Q transform event structures
%   channelNames         cell array of transform channel names
%
% QTHRESHOLD returns a cell array of Q transform event structures that
% contain the properties of the identified statistically significant tiles
% for each channel.  The event structure contains the following fields.
%
%   time                 center time of tile [gps seconds]
%   frequency            center frequency of tile [Hz]
%   q                    quality factor of tile []
%   duration             duration of tile [seconds]
%   bandwidth            bandwidth of tile [Hz]
%   normalizedEnergy     normalized energy of tile []
%   amplitude            signal amplitude of tile [Hz^-1/2]
%   overflowFlag         boolean overflow flag
%
% For collocated transform data, the following field is also returned.
%
%   incoherentEnergy     incoherent energy of tile []
%
% The user can focus on a subset of the times and frequencies available in
% the transform data by specifying a desired range of central times,
% central frequencies, and Qs to threshold on.  Ranges should be specified
% as a two component vector, consisting of a minimum and maximum value.
% Alternatively, if only a single Q is specified, QTHRESHOLD is only
% applied to the time-frequency plane which has the nearest value of Q in a
% logarithmic sense to the requested value.
%
% To determine the range of central times to threshold on, QTHRESHOLD
% requires the start time of the transformed data in addition to a
% reference time and a relative time range.  Both the start time and
% reference time should be specified as absolute quantities, while the
% range of times to plot should be specified relative to the requested
% reference time.
%
% By default, QTHRESHOLD is applied to all available frequencies and Qs,
% and the reference time and relative time range arguments are set to
% exclude data potentially corrupted by filter transients as identified by
% the transient duration field of the tiling structure.  The default value
% can be obtained for any argument by passing the empty matrix [].
%
% The threshold is set to yield the specified false event rate when applied
% to all available frequencies and Qs, and is not modified to account for
% restricted ranges.  It is also only a rough estimate, and the result
% false event rate may vary significantly depending on the quality of the
% data.
%
% If provided, the optional analysisMode string is used by QTHRESHOLD to
% determine which channels are signal channels, which channels are veto
% channels, and which channels to report results for.  Analysis modes which
% implement null stream based vetoes must also specify a target white noise
% false veto rate, and a squared calibration uncertainty factor.  If the
% optional applyVeto flag is set to true, then QTHRESHOLD also applies the
% veto to the signal channel in addition to simply thresholding the veto
% channel.
%
% The optional maximumSignificants argument provides a safety mechanism to
% limit the total number of events returned by QTHRESHOLD.  If this maximum
% number of significants is exceeded, the overflow flag is set, only the
% maximumSignificants most significant tiles are returned, and a warning is
% issued if debugLevel is set to 1 or higher.  By default, maximumSignificants
% is set to infinity and debugLevel is set to unity.
%
% See also QTILE, QCONDITION, QTRANSFORM, QSELECT, QEXAMPLE, and QPIPELINE.

% Shourov K. Chatterji
% shourov@ligo.caltech.edu

% Leo C. Stein
% lstein@ligo.mit.edu

% $Id: qthreshold.m,v 1.18 2007/10/21 12:05:14 shourov Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(4, 16, nargin));

% apply default arguments
if (nargin < 5) || isempty(referenceTime),
  referenceTime = startTime + tiling.duration / 2;
end
if (nargin < 6) || isempty(timeRange),
  timeRange = 0.5 * (tiling.duration - 2 * tiling.transientDuration) * [-1 +1];
end
if (nargin < 7) || isempty(frequencyRange),
  frequencyRange = [-Inf +Inf];
end
if (nargin < 8) || isempty(qRange),
  qRange = [-Inf +Inf];
end
if (nargin < 9) || isempty(maximumSignificants),
  maximumSignificants = Inf;
end
if (nargin < 10) || isempty(analysisMode),
  analysisMode = 'single';
end
if (nargin < 11) || isempty(falseVetoRate),
  falseVetoRate = 0;
end
if (nargin < 12) || isempty(uncertaintyFactor),
  uncertaintyFactor = 0;
end
if (nargin < 13) || isempty(correlationFactor),
  correlationFactor = 0;
end
if (nargin < 14) || isempty(applyVeto),
  applyVeto = false;
end
if (nargin < 15) || isempty(channelNames),
  channelNames = [];
end
if (nargin < 16) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
if ~iscell(transforms),
  transforms = mat2cell(transforms, size(transforms, 1), ...
                        size(transforms, 2));
end
if ~isempty(channelNames) & ~iscell(channelNames),
  channelNames = mat2cell(channelNames, size(channelNames, 1), ...
                          size(channelNames, 2));
end

% force one dimensional cell arrays
transforms = transforms(:);
channelNames = channelNames(:);

% determine number of channels
numberOfChannels = length(transforms);

% provide default channel names
if isempty(channelNames),
  channelNames = cell(numberOfChannels, 1);
  for channelNumber = 1 : numberOfChannels,
    channelNames{channelNumber} = ['X' int2str(channelNumber) ':STRAIN'];
  end
end

% force ranges to be monotonically increasing column vectors
timeRange = unique(timeRange(:));
frequencyRange = unique(frequencyRange(:));
qRange = unique(qRange(:));

% if only a single Q is requested, find nearest Q plane
if length(qRange) == 1,
  [ignore, qPlane] = min(abs(log(tiling.qs / qRange)));
  qRange = tiling.qs(qPlane) * [1 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate tiling structure
if ~strcmp(tiling.id, 'Discrete Q-transform tile structure'),
  error('input argument is not a discrete Q transform tiling structure');
end

% validate transform structures
for channelNumber = 1 : numberOfChannels,
  if ~strcmp(transforms{channelNumber}.id, ...
             'Discrete Q-transform transform structure'),
    error('input argument is not a discrete Q transform structure');
  end
end

% validate number of channels
if length(channelNames) ~= numberOfChannels,
  error('number of channelNames is inconsistent with number of transforms');
end

% Check for two component range vectors
if length(timeRange) ~= 2,
  error('Time range must be two component vector [tmin tmax].');
end
if length(frequencyRange) ~= 2,
  error('Frequency range must be two component vector [fmin fmax].');
end
if length(qRange) > 2,
  error('Q range must be scalar or two component vector [Qmin Qmax].');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         normalized energy threshold                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% approximate number of statistically independent tiles per second
independentsRate = tiling.numberOfIndependents / tiling.duration;

% apply emperically determined correction factor
independentsRate = independentsRate * 1.5;

% probability associated with desired false event rate
falseEventProbability = falseEventRate / independentsRate;

% probability associated with desired false veto rate
falseVetoProbability = falseVetoRate / independentsRate;

% normalized energy threshold for desired false event rate
eventThreshold = -log(falseEventProbability);

% normalized energy threshold for desired false veto rate
if falseVetoProbability == 0,
  vetoThreshold = Inf;
else
  vetoThreshold = -log(falseVetoProbability);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             apply analysis mode                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% switch on analysis mode
switch lower(analysisMode),

  % handle case of single detector and coincident analysis
  case {'single', 'coincident'},

    % threshold on all signal channels individually
    for channelNumber = 1 : numberOfChannels,
      outputChannels{channelNumber}.channelName = channelNames{channelNumber};
      outputChannels{channelNumber}.signalChannels = channelNumber;
      outputChannels{channelNumber}.signalField = 'normalizedEnergies';
      outputChannels{channelNumber}.nullChannels = [];
      outputChannels{channelNumber}.nullField = [];
      outputChannels{channelNumber}.referenceChannels = [];
      outputChannels{channelNumber}.referenceField = [];
    end

  % handle case of h1h2 analysis
  case 'h1h2',

    % threshold on signal channel and apply null stream veto
    outputChannels{1}.channelName = channelNames{1};
    outputChannels{1}.signalChannels = 1;
    outputChannels{1}.signalField = 'normalizedEnergies';
    outputChannels{1}.nullChannels = 3;
    outputChannels{1}.nullField = 'normalizedEnergies';
    outputChannels{1}.referenceChannels = 4;
    outputChannels{1}.referenceField = 'normalizedEnergies';

    % threshold on null stream veto channel
    outputChannels{2}.channelName = channelNames{3};
    outputChannels{2}.signalChannels = [];
    outputChannels{2}.signalField = 'normalizedEnergies';
    outputChannels{2}.nullChannels = 3;
    outputChannels{2}.nullField = 'normalizedEnergies';
    outputChannels{2}.referenceChannels = 4;
    outputChannels{2}.referenceField = 'normalizedEnergies';

  % handle unknown analysis mode
  otherwise,
    error(['unknown analysis mode "' analysisMode '"']);

% end switch on analysis mode
end

% number of output channels
numberOfOutputChannels = length(outputChannels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             initialize statistically significant event structure             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of significant event structures
significants = cell(numberOfOutputChannels, 1);

% create empty cell array of output channel names
channelNames = cell(numberOfOutputChannels, 1);

% begin loop over channels
for outputChannelNumber = 1 : numberOfOutputChannels

  % insert structure identification string
  significants{outputChannelNumber}.id = 'Discrete Q-transform event structure';

  % initialize result vectors
  significants{outputChannelNumber}.time = [];
  significants{outputChannelNumber}.duration = [];
  significants{outputChannelNumber}.frequency = [];
  significants{outputChannelNumber}.bandwidth = [];
  significants{outputChannelNumber}.q = [];
  significants{outputChannelNumber}.normalizedEnergy = [];
  significants{outputChannelNumber}.amplitude = [];

  % initialize overflow flag
  significants{outputChannelNumber}.overflowFlag = 0;

  % include incoherent energy for veto channels
  if isempty(outputChannels{outputChannelNumber}.signalChannels),
    significants{outputChannelNumber}.incoherentEnergy = [];
  end

  % fill cell array of output channel names
  channelNames{outputChannelNumber} = ...
      outputChannels{outputChannelNumber}.channelName;

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin loop over Q planes                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over Q planes
for plane = 1 : tiling.numberOfPlanes,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                              threshold on Q                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % skip Q planes outside of requested Q range
  if ((tiling.planes{plane}.q < ...
       min(qRange)) | ...
      (tiling.planes{plane}.q > ...
       max(qRange))),
    continue;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      begin loop over frequency rows                        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % begin loop over frequency rows
  for row = 1 : tiling.planes{plane}.numberOfRows,

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    threshold on central frequency                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % skip frequency rows outside of requested frequency range
    if ((tiling.planes{plane}.rows{row}.frequency < ...
         min(frequencyRange)) | ...
        (tiling.planes{plane}.rows{row}.frequency > ...
         max(frequencyRange))),
      continue;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       begin loop over channels                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % begin loop over channels
    for outputChannelNumber = 1 : numberOfOutputChannels,

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                  extract output channel details                        %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % extract output channel structure
      outputChannel = outputChannels{outputChannelNumber};

      % extract output channel details
      channelName = outputChannel.channelName;
      signalChannels = outputChannel.signalChannels;
      signalField = outputChannel.signalField;
      nullChannels = outputChannel.nullChannels;
      nullField = outputChannel.nullField;
      referenceChannels = outputChannel.referenceChannels;
      referenceField = outputChannel.referenceField;

      % number of null channels
      numberOfNullChannels = length(nullChannels);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                    threshold on significance                           %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % if signal channel,
      if ~isempty(signalChannels),

        % threshold on signal significance
        significantTileIndices = find( ...
            transforms{signalChannels}.planes{plane}.rows{row} ...
              .(signalField) >=  ...
            eventThreshold);

        % apply same tile veto
        if applyVeto,
          for nullChannelNumber = 1 : numberOfNullChannels,
            nullChannel = nullChannels(nullChannelNumber);
            referenceChannel = referenceChannels(nullChannelNumber);
            keepIndices = find( ...
                transforms{nullChannel}.planes{plane}.rows{row} ...
                  .(nullField)(significantTileIndices) < ...
                vetoThreshold + uncertaintyFactor * ...
                transforms{referenceChannel}.planes{plane}.rows{row} ...
                  .(referenceField)(significantTileIndices));
            significantTileIndices = significantTileIndices(keepIndices);
          end
        end

      % if veto channel,
      else

        % threshold on veto significance
        significantTileIndices = find( ...
            transforms{nullChannels}.planes{plane}.rows{row} ...
              .(nullField) >= ...
            vetoThreshold + uncertaintyFactor * ...
            transforms{referenceChannels}.planes{plane}.rows{row} ...
              .(referenceField));

      % end test for veto or signal channel
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                    threshold on central time                           %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % skip tiles outside requested time range
      keepIndices = ...
          (tiling.planes{plane}.rows{row}.times(significantTileIndices) >= ...
           (referenceTime - startTime + min(timeRange))) & ...
          (tiling.planes{plane}.rows{row}.times(significantTileIndices) <= ...
           (referenceTime - startTime + max(timeRange)));
      significantTileIndices = significantTileIndices(keepIndices);

      % number of statistically significant tiles in frequency row
      numberOfSignificants = length(significantTileIndices);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %      append significant tile properties to event structure             %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % append center times of significant tiles in row
      significants{outputChannelNumber}.time = ...
          [significants{outputChannelNumber}.time ...
           tiling.planes{plane}.rows{row}.times(significantTileIndices) + ...
           startTime];

      % append center frequencies of significant tiles in row
      significants{outputChannelNumber}.frequency = ...
          [significants{outputChannelNumber}.frequency ...
           tiling.planes{plane}.rows{row}.frequency * ...
           ones(1, numberOfSignificants)];

      % append qs of significant tiles in row
      significants{outputChannelNumber}.q = ...
          [significants{outputChannelNumber}.q ...
           tiling.planes{plane}.q * ...
           ones(1, numberOfSignificants)];

      % append durations of significant tiles in row
      significants{outputChannelNumber}.duration = ...
          [significants{outputChannelNumber}.duration ...
           tiling.planes{plane}.rows{row}.duration * ...
           ones(1, numberOfSignificants)];

      % append bandwidths of significant tiles in row
      significants{outputChannelNumber}.bandwidth = ...
          [significants{outputChannelNumber}.bandwidth ...
           tiling.planes{plane}.rows{row}.bandwidth * ...
           ones(1, numberOfSignificants)];

      % if signal channel,
      if ~isempty(signalChannels),
        
        % append normalized energies of significant tiles in row
        significants{outputChannelNumber}.normalizedEnergy = ...
            [significants{outputChannelNumber}.normalizedEnergy ...
             (transforms{signalChannels}.planes{plane}.rows{row} ...
                .(signalField)(significantTileIndices))];

        % append amplitudes of significant tiles in row
        significants{outputChannelNumber}.amplitude = ...
            [significants{outputChannelNumber}.amplitude ...
             sqrt(transforms{signalChannels}.planes{plane}.rows{row} ...
                  .energies(significantTileIndices) - ...
                  transforms{signalChannels}.planes{plane}.rows{row}.meanEnergy)];

      % if veto channel,
      else

        % append normalized energies of significant tiles in row
        significants{outputChannelNumber}.normalizedEnergy = ...
            [significants{outputChannelNumber}.normalizedEnergy ...
             (transforms{nullChannels}.planes{plane}.rows{row} ...
                .(nullField)(significantTileIndices))];

        % append amplitudes of significant tiles in row
        significants{outputChannelNumber}.amplitude = ...
            [significants{outputChannelNumber}.amplitude ...
             sqrt(transforms{nullChannels}.planes{plane}.rows{row} ...
                  .energies(significantTileIndices) - ...
                  transforms{nullChannels}.planes{plane}.rows{row}.meanEnergy)];

        % append incoherent energies of significant tiles in row,
        significants{outputChannelNumber}.incoherentEnergy = ...
            [significants{outputChannelNumber}.incoherentEnergy ...
             (transforms{referenceChannels}.planes{plane}.rows{row} ...
                .(referenceField)(significantTileIndices))];

      % end test for veto or signal channel
      end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        end loop over channels                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % end loop over channels
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                       end loop over frequency rows                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % end loop over frequency rows
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            end loop over Q planes                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over Q planes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               check for excessive number of significant tiles                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for outputChannelNumber = 1 : numberOfOutputChannels

  % determine number of significant tiles in channel
  numberOfSignificants = length(significants{outputChannelNumber}.time);

  % if maximum allowable number of significant tiles is exceeded
  if numberOfSignificants > maximumSignificants,

    % issue warning
    if debugLevel > 0,
      fprintf(1, 'WARNING: %s: maximum number of significants exceeded.\n', ...
              outputChannels{outputChannelNumber}.channelName);
    end

    % set overflow flag
    significants{outputChannelNumber}.overflowFlag = 1;

    % sort significant tiles by normalized energy
    [ignore, maximumIndices] = ...
        sort(significants{outputChannelNumber}.normalizedEnergy);

    % find indices of most significant tiles
    maximumIndices = maximumIndices(end - maximumSignificants + 1 : end);

    % get list of field names available for this channel
    fieldNames = fieldnames(significants{outputChannelNumber});

    % remove id and overflowFlag fields from list of field names
    fieldNames = fieldNames(~strcmp(fieldNames, 'id') & ...
                            ~strcmp(fieldNames, 'overflowFlag'));

    % number of field names available for this channel
    numberOfFields = length(fieldNames);

    % extract most significant tile properties
    for fieldNumber = 1 : numberOfFields,
      significants{outputChannelNumber}.(fieldNames{fieldNumber}) = ...
          significants{outputChannelNumber} ...
            .(fieldNames{fieldNumber})(maximumIndices);
    end

  % otherwise continue
  end

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    return statistically significant tiles                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return
