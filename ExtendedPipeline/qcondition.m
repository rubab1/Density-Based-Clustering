function [data, calibration, coefficients] = ...
    qcondition(data, tiling, responseFunctions);
% QCONDITION Calibrate, high pass filter and whiten time series data
%
% QCONDITION high pass filters and whitens time series data prior to analysis
% by the Q transform.  The data is first zero-phase high pass filtered at the
% minimum analysis frequency specified in the tiling structure.  The resulting
% high pass filtered data is then whitened by zero-phase linear prediction at
% a frequency resolution equal to the minimum analysis bandwidth requested in
% the tiling structure.  Note that the resulting whitened data is returned as
% a frequency series, not a time series.  In addition, the returned frequency
% series extend from zero frequency to the Nyquist frequency.  As a result,
% they are of length N / 2 + 1, where N is the length of the individual input
% time series.
%
% To recover the characteristic amplitude of candidate signals, QCONDITION also
% returns calibration factors that take into account the effect of the high pass
% and whitening filters on the magnitude of Q transform coefficients, as well
% as any calibration information that is provided.  QCONDITION also returns
% the effective frequency domain coefficients of the combined high pass and
% whitening filters for each channel.
%
% usage:
%
% [data, calibration, coefficients] = ...
%     qcondition(data, tiling, responseFunctions);
%
%   data                  cell array of input time series
%   tiling                discrete Q transform tiling structure from QTILE
%   responseFunctions     cell array of calibration response functions
%
%   data                  cell array of conditioned output frequency series
%   calibration           cell array of Q transform calibration structures
%   coefficients          cell array of frequency domain filter coefficients
%
% The resulting calibration structures are parallel to the structure returned
% by QTILE and contain the following supplemental fields for each frequency
% row.
%
%   calibrationFactor     row specific amplitude correction factors
%
% If provided, responseFunctions should consist of a cell array of complex
% valued frequency response vectors, with one cell per channel.  These vectors
% should extend from zero frequency to the Nyquist frequency at a frequency
% resolution equal to the reciprocal of the data duration.  As a result they
% should be of length N / 2 + 1, where N is the length of the individual input
% time series.  By default, QCONDITION assumes that the input data is already
% calibrated and that responseFunctions is unity.
%
% See also QTILE, QTRANSFORM, QTHRESHOLD, QSELECT, and QPIPELINE.

% Shourov K. Chatterji
% shourov@ligo.mit.edu

% $Id: qcondition.m,v 1.10 2007/10/21 12:05:14 shourov Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        process command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify correct number of input arguments
error(nargchk(2, 3, nargin));

% apply default arguments
if (nargin < 3) || isempty(responseFunctions),
  responseFunctions = [];
end

% provide default response functions
if isempty(responseFunctions),
  responseFunctions = cell(size(data));
end

% force cell arrays
if ~iscell(data),
  data = mat2cell(data, size(data, 1), size(data, 2));
end
if ~iscell(responseFunctions),
  responseFunctions = mat2cell(responseFunctions, ...
                               size(responseFunctions, 1), ...
                               size(responseFunctions, 2));
end

% force one dimensional cell arrays
data = data(:);
responseFunctions = responseFunctions(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       validate command line arguments                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of channels
numberOfChannels = length(data);

% validate tiling structure
if ~strcmp(tiling.id, 'Discrete Q-transform tile structure'),
  error('input argument is not a discrete Q transform tiling structure');
end

% determine required data lengths
dataLength = tiling.sampleFrequency * tiling.duration;
halfDataLength = dataLength / 2 + 1;

% validate data length and force row vectors
for channelNumber = 1 : numberOfChannels,
  data{channelNumber} = data{channelNumber}(:).';
  if length(data{channelNumber}) ~= dataLength,
    error('data length not consistent with tiling');
  end
end

% validate calibration response channels
if length(responseFunctions) ~= numberOfChannels,
  error('number of channels not consistent with calibration response');
end

% validate calibration response length and force row vectors
for channelNumber = 1 : numberOfChannels,
  if ~isempty(responseFunctions{channelNumber}),
    responseFunctions{channelNumber} = ...
        responseFunctions{channelNumber}(:).';
    if length(responseFunctions{channelNumber}) ~= halfDataLength,
      error('calibration response length not consistent with tiling');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                initialize Q transform calibration structures                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create empty cell array of Q transform calibration structures
calibration = cell(size(data));

% create empty cell array of conditioning filter coefficients
coefficients = cell(size(data));

% begin loop over channels
for channelNumber = 1 : numberOfChannels,

  % insert structure identification string
  calibration{channelNumber}.id = 'Discrete Q-transform calibration structure';

  % create empty cell array of Q plane structures
  calibration{channelNumber}.planes = cell(size(tiling.planes));

  % begin loop over Q planes
  for plane = 1 : tiling.numberOfPlanes,

    % create empty cell array of frequency row structures
    calibration{channelNumber}.planes{plane}.rows = ...
        cell(size(tiling.planes{plane}.numberOfRows));

  % end loop over Q planes
  end

% end loop over channels
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             hardcoded parameters                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% high pass filter order
hpfOrder = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              derived parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nyquist frequency
nyquistFrequency = tiling.sampleFrequency / 2;

% linear predictor error filter order
lpefOrder = ceil(tiling.sampleFrequency * tiling.whiteningDuration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 window data                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% design window
% tukeyWindow = ...
%     tukeywin(dataLength, 2 * ...
%              max(tiling.transientDuration - tiling.whiteningDuration, 0) / ...
%              tiling.duration).';

% apply window
% for channelNumber = 1 : numberOfChannels,
%   data{channelNumber} = data{channelNumber} .* tukeyWindow;
% end

% release window memory
% clear tukeyWindow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               high pass filter                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if high pass filtering is requested,
if tiling.highPassCutoff > 0,

  % design high pass filter
  [hpfZeros, hpfPoles, hpfGain] = butter(hpfOrder, tiling.highPassCutoff / ...
                                         nyquistFrequency, 'high');
  hpfSOS = zp2sos(hpfZeros, hpfPoles, hpfGain);

  % apply high pass filter
  for channelNumber = 1 : numberOfChannels,
    data{channelNumber} = sosfiltfilt(hpfSOS, data{channelNumber});
  end

  % magnitude response of high pass filter
  minimumFrequencyStep = 1 / tiling.duration;
  frequencies = 0 : minimumFrequencyStep : nyquistFrequency;
  hpfArgument = (frequencies / tiling.highPassCutoff).^(2 * hpfOrder);
  hpfResponse = hpfArgument ./ (1 + hpfArgument);

% end test for high pass filtering
end

% supress high pass filter transients
for channelNumber = 1 : numberOfChannels,
  data{channelNumber}(1 : lpefOrder) = ...
      zeros(1, lpefOrder);
  data{channelNumber}(dataLength - lpefOrder + 1 : dataLength) = ...
      zeros(1, lpefOrder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            fast fourier transform                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fourier transform
for channelNumber = 1 : numberOfChannels,
  data{channelNumber} = fft(data{channelNumber});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        linear predictor error filter                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if linear predictive whitening is requested,
if tiling.whiteningDuration > 0,

  % auto-correlation estimates
  for channelNumber = 1 : numberOfChannels,
    autoCorrelation{channelNumber} = real(data{channelNumber}).^2 + ...
                                     imag(data{channelNumber}).^2;
  end
  for channelNumber = 1 : numberOfChannels,
    autoCorrelation{channelNumber} = real(ifft(autoCorrelation{channelNumber}));
  end

  % solve for lpef coefficients
  filterLength = lpefOrder + 1;
  for channelNumber = 1 : numberOfChannels,
    lpefCoefficients{channelNumber} = ...
        levinson(autoCorrelation{channelNumber}(1 : filterLength), lpefOrder);
  end

  % release auto-correlation vector memory
  clear autoCorrelation;

  % create zero-phase lpef
  for channelNumber = 1 : numberOfChannels,
    lpefCoefficients{channelNumber} = [lpefCoefficients{channelNumber} ...
                                       zeros(1, dataLength - filterLength)];
  end
  for channelNumber = 1 : numberOfChannels,
    lpefCoefficients{channelNumber} = fft(lpefCoefficients{channelNumber});
  end
  for channelNumber = 1 : numberOfChannels,
    lpefCoefficients{channelNumber} = abs(lpefCoefficients{channelNumber});
  end

  % extract one sided frequency domain data and lpef coefficients
  for channelNumber = 1 : numberOfChannels,
    data{channelNumber} = data{channelNumber}(1 : halfDataLength);
    lpefCoefficients{channelNumber} = ...
        lpefCoefficients{channelNumber}(1 : halfDataLength);
  end

  % apply zero-phase lpef
  for channelNumber = 1 : numberOfChannels,
    data{channelNumber} = lpefCoefficients{channelNumber} .* data{channelNumber};
  end

  % renormalize to unity spectral density
  for channelNumber = 1 : numberOfChannels,
    normalization = sqrt(dataLength * nyquistFrequency * halfDataLength / ...
                         sum(abs(data{channelNumber}).^2));
    data{channelNumber} = ...
        normalization * data{channelNumber};
    lpefCoefficients{channelNumber} = ...
        normalization * lpefCoefficients{channelNumber};
  end

% if linear predictive whitening is not requested,
else

  % extract one sided frequency domain data
  for channelNumber = 1 : numberOfChannels,
    data{channelNumber} = data{channelNumber}(1 : halfDataLength);
  end

% end test for linear predictive whitening
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             apply calibration                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute calibration correction factors
for channelNumber = 1 : numberOfChannels,
  coefficients{channelNumber} = ones(size(data{channelNumber}));
  if tiling.highPassCutoff > 0,
    coefficients{channelNumber} = ...
        coefficients{channelNumber} .* hpfResponse;
  end
  if tiling.whiteningDuration > 0,
    coefficients{channelNumber} = ...
        coefficients{channelNumber} .* lpefCoefficients{channelNumber};
  end
  if ~isempty(responseFunctions{channelNumber}),
    coefficients{channelNumber} = coefficients{channelNumber} ./ ...
        abs(responseFunctions{channelNumber});
  end
  % smoothingBandwidth = 100 / tiling.duration;
  % smoothingLength = ceil(smoothingBandwidth / minimumFrequencyStep);
  % coefficients{channelNumber} = ...
  %     smooth(coefficients{channelNumber}, smoothingLength, 'moving').';
  for plane = 1 : tiling.numberOfPlanes,
    for row = 1 : tiling.planes{plane}.numberOfRows,
      calibrationIntegrand = (tiling.planes{plane}.rows{row}.window ./ ...
          coefficients{channelNumber}(tiling.planes{plane}.rows{row}.dataIndices)).^2;
      calibration{channelNumber}.planes{plane}.rows{row}.calibrationFactor = ...
          sqrt(sum(calibrationIntegrand * minimumFrequencyStep) / 2) * ...
          dataLength / tiling.planes{plane}.rows{row}.numberOfTiles;
    end
  end
end

% apply calibration phase
for channelNumber = 1 : numberOfChannels,
  if ~isempty(responseFunctions{channelNumber}),
    data{channelNumber} = data{channelNumber} .* ...
        exp(sqrt(-1) * angle(responseFunctions{channelNumber}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    return                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to calling function
return;
