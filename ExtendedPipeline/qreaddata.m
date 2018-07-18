function [data, sampleFrequencies] = ...
    qreaddata(frameCache, channelNames, frameTypes, startTime, stopTime, ...
              timeShifts, debugLevel)
% QREADDATA Read multiple channels of data from frame files
%
% QREADDATA finds and retrieves multiple channels of time series data from a set
% of frame files.  The data are specified by the frame file type, channel name,
% start time, and duration or stop time.  The necessary frame files are located
% using a file name caching scheme.
%
% QREADDATA is build on top of the READFRAMEDATA set of high level functions
% for reading data from frame files.
%
% usage:
%
%   [data, sampleFrequencies] = ...
%     qreaddata(frameCache, channelNames, frameTypes, startTime, stopTime, ...
%               timeShifts, debugLevel);
%
%   frameCache           frame file cache structure from LOADFRAMECACHE
%   channelNames         cell array of channel names
%   frameTypes           cell array of frame types
%   startTime            gps start time of data to extract
%   stopTime             gps stop time (or duration) of data to extract
%   timeShifts           vector of time shifts [seconds]
%   debugLevel           verboseness of debug output
%
%   data                 cell array of extracted data
%   sampleFrequencies    vector of sample frequencies
%
% The data is returned as a cell array of row vectors in the same order as in
% the cell array of channel names.
%
% QREADDATA retrieves data from the requested start time up to, but not
% including, the requested stop time, such that stop minus start seconds are
% retrieved.  Alternatively, the desired duration in seconds may be specified
% instead of the GPS stop time parameter.
%
% The optional time shift argument should be a vector with one element per
% channel that specifies an additional time shift to apply to the start time for
% each detector.  By convention, a positive time shift corresponds to a delay of
% the corresponding time series.  If no time offsets are specified, a default
% time shift of zero is assumed for all detectors.
%
% To avoid the effects of round-off error, the requested start time, stop
% time, and time shifts should all be aligned with samples of the input data
% stream.
%
% If it is unable to load the requested data, QREADDATA returns empty result
% vectors and zero sample frequency for the failed channel as well as issuing a
% warning if debugLevel is set to 1 or higher.  By default, a debugLevel of
% unity is assumed.
%
% See also READFRAMEDATA, LOADFRAMECACHE, CREATEFRAMECACHE, FRGETVECT,
% QPIPELINE, and QEXAMPLE.

% Shourov K. Chatterji
% shourov@ligo.caltech.edu

% $Id: qreaddata.m,v 1.6 2007/04/03 13:18:21 shourov Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            hard coded parameters                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% permit redundant frame data
allowRedundantFlag = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for sufficient command line arguments
error(nargchk(5, 7, nargin));

% apply default arguments
if (nargin < 6) || isempty(timeShifts),
  timeShifts = [];
end
if (nargin < 7) || isempty(debugLevel),
  debugLevel = 1;
end

% force cell arrays
if ~iscell(channelNames),
  channelNames = mat2cell(channelNames, size(channelNames, 1), ...
                          size(channelNames, 2));
end
if ~iscell(frameTypes),
  frameTypes = mat2cell(frameTypes, size(frameTypes, 1), ...
                        size(frameTypes, 2));
end

% provide default time shifts
if isempty(timeShifts),
  timeShifts = zeros(size(channelNames));
end

% force one dimensional cell arrays and vectors
channelNames = channelNames(:);
frameTypes = frameTypes(:);
timeShifts = timeShifts(:);

% determine number of channels
numberOfChannels = length(channelNames);

% convert duration to absolute stop time
if stopTime < startTime,
  stopTime = startTime + stopTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate frame types
if length(frameTypes) ~= numberOfChannels,
  error('number of frame types is inconsistent with number of channels');
end

% validate time shifts
if length(timeShifts) ~= numberOfChannels,
  error('number of time shifts is inconsistent with number of channels');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         initialize result structures                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize result structures
data = cell(numberOfChannels, 1);
sampleFrequencies = zeros(numberOfChannels, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  read data                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  % switch on frametype
  switch upper(frameTypes{channelNumber}),

    % if no data is requested,
    case 'NONE',

     % return empty result
     data{channelNumber} = [];
     sampleFrequencies(channelNumber) = 0;

    % if whitenoise is requested,
    case {'WHITE', 'WHITENOISE'},

      % generate white noise data
      sampleFrequencies(channelNumber) = 16384;
      data{channelNumber} = sqrt(sampleFrequencies(channelNumber) / 2) * ...
                            normrnd(0, 1, 1, ...
                                    sampleFrequencies(channelNumber) * ...
                                    (stopTime - startTime));

    % otherwise,
    otherwise,

      % read requested data
      [data{channelNumber}, sampleFrequencies(channelNumber)] = ...
          readframedata(frameCache, channelNames{channelNumber}, ...
                        frameTypes{channelNumber}, ...
                        startTime - timeShifts(channelNumber), ...
                        stopTime - timeShifts(channelNumber), ...
                        allowRedundantFlag, debugLevel);

  % end switch on frametype
  end

  % force row vector data
  data{channelNumber} = data{channelNumber}(:)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          return to calling function                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
