function [pathGainCell, aux] = buildFramePathGain(numFrame, numSat, numUser, pathOpt)
%BUILDFRAMEPATHGAIN Build frame-wise complex path-gain cells.
% Generates a 1xNf cell array of Ns x Nu complex path-gain matrices for
% multi-frame snapshot simulation. The function is intended to replace the
% ad-hoc localBuildFramePathGain helpers in the performance scripts while
% keeping the path-gain container format accepted by genMultiFrameSnapshots.
%
%Syntax:
%  pathGainCell = buildFramePathGain(numFrame, numSat, numUser)
%
%  pathGainCell = buildFramePathGain(numFrame, numSat, numUser, ampSpec)
%
%  [pathGainCell, aux] = buildFramePathGain(numFrame, numSat, numUser, pathOpt)
%
%Inputs:
%  numFrame        - number of frames
%
%  numSat          - number of satellites
%
%  numUser         - number of users
%
%  ampSpec         - optional amplitude specification. This shorthand is
%                    equivalent to pathOpt.ampSpec. Supported layouts:
%                    - scalar                : shared by all links
%                    - Nu x 1 / 1 x Nu       : one amplitude per user
%                    - Ns x 1 / 1 x Ns       : one amplitude per satellite
%                    - Ns x Nu               : one amplitude per link
%                    - 1 x Nf cell           : one supported layout per frame
%
%  pathOpt         - optional options structure
%    .ampSpec         - amplitude specification, default 1
%    .phaseMode       - frame-phase mode
%                       - 'shared'      : one phase per link shared by all
%                                         frames (default)
%                       - 'independent' : one phase per frame and per link
%                       - 'zero'        : zero phase for all blocks
%                       - 'custom'      : use phaseSpec directly
%
%    .phaseSpec       - optional phase specification in radians.
%                       Supported layouts match ampSpec. When phaseMode is
%                       'shared' or 'custom', a non-cell input is reused by
%                       all frames. When phaseMode is 'independent' or
%                       'custom', a 1xNf cell may be used to specify one
%                       phase matrix per frame.
%
%    .ampJitterStdDb - optional per-frame amplitude jitter standard
%                       deviation in dB. The jitter is applied
%                       multiplicatively and independently to each frame.
%                       default 0
%
%    .frameMask      - optional valid-block mask
%                       - scalar logical : applied to all sat-frame blocks
%                       - 1 x Nf logical : one mask per frame
%                       - Ns x Nf logical: one mask per satellite-frame
%                       Blocks masked out are set to zero.
%
%Outputs:
%  pathGainCell     - 1 x Nf cell, each cell is an Ns x Nu complex matrix
%
%  aux              - auxiliary structure with fields
%    .numFrame
%    .numSat
%    .numUser
%    .phaseMode
%    .frameMask
%    .ampCell
%    .phaseCell
%    .pathGainMat   : Ns x Nu x Nf complex array
%
%Notes:
%  - The default output uses unit amplitudes and a frame-shared random
%    phase per satellite-user link, which is usually the appropriate
%    default for continuous-phase experiments.
%  - To reproduce the older helper that generated one random phase per
%    frame, use pathOpt.phaseMode = 'independent'.
%
%See also:
%  genMultiFrameSnapshots, genMultiSatSnapshots

if nargin < 4 || isempty(pathOpt)
  pathOpt = struct();
elseif ~isstruct(pathOpt)
  pathOpt = struct('ampSpec', pathOpt);
end

validateattributes(numFrame, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
  mfilename, 'numFrame');
validateattributes(numSat, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
  mfilename, 'numSat');
validateattributes(numUser, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
  mfilename, 'numUser');

pathOpt = localParsePathOpt(pathOpt, numFrame, numSat, numUser);

ampCell = localBuildAmpCell(pathOpt.ampSpec, pathOpt.ampJitterStdDb, numFrame, numSat, numUser);
phaseCell = localBuildPhaseCell(pathOpt.phaseMode, pathOpt.phaseSpec, numFrame, numSat, numUser);

pathGainCell = cell(1, numFrame);
pathGainMat = zeros(numSat, numUser, numFrame);

for iFrame = 1:numFrame
  currentGain = ampCell{iFrame} .* exp(1j * phaseCell{iFrame});

  satMask = pathOpt.frameMask(:, iFrame);
  if any(~satMask)
    currentGain(~satMask, :) = 0;
  end

  pathGainCell{iFrame} = currentGain;
  pathGainMat(:, :, iFrame) = currentGain;
end

aux = struct();
aux.numFrame = numFrame;
aux.numSat = numSat;
aux.numUser = numUser;
aux.phaseMode = pathOpt.phaseMode;
aux.frameMask = pathOpt.frameMask;
aux.ampCell = ampCell;
aux.phaseCell = phaseCell;
aux.pathGainMat = pathGainMat;
end


function pathOpt = localParsePathOpt(pathOpt, numFrame, numSat, numUser)
%LOCALPARSEPATHOPT Parse path-gain generation options.

if ~isfield(pathOpt, 'ampSpec') || isempty(pathOpt.ampSpec)
  pathOpt.ampSpec = 1;
end

if ~isfield(pathOpt, 'phaseMode') || isempty(pathOpt.phaseMode)
  pathOpt.phaseMode = 'shared';
end
pathOpt.phaseMode = localNormalizePhaseMode(pathOpt.phaseMode);

if ~isfield(pathOpt, 'phaseSpec')
  pathOpt.phaseSpec = [];
end

if ~isfield(pathOpt, 'ampJitterStdDb') || isempty(pathOpt.ampJitterStdDb)
  pathOpt.ampJitterStdDb = 0;
end
validateattributes(pathOpt.ampJitterStdDb, {'numeric'}, ...
  {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'pathOpt.ampJitterStdDb');

if ~isfield(pathOpt, 'frameMask') || isempty(pathOpt.frameMask)
  pathOpt.frameMask = true(numSat, numFrame);
else
  pathOpt.frameMask = localNormalizeFrameMask(pathOpt.frameMask, numFrame, numSat);
end

localValidateAmpSpec(pathOpt.ampSpec, numFrame, numSat, numUser);
localValidatePhaseSpec(pathOpt.phaseSpec, numFrame, numSat, numUser, pathOpt.phaseMode);
end


function phaseMode = localNormalizePhaseMode(phaseMode)
%LOCALNORMALIZEPHASEMODE Normalize path phase mode string.

phaseMode = lower(string(phaseMode));
phaseMode = char(phaseMode);

validMode = {'shared', 'independent', 'zero', 'custom'};
if ~ismember(phaseMode, validMode)
  error('buildFramePathGain:InvalidPhaseMode', ...
    'pathOpt.phaseMode must be ''shared'', ''independent'', ''zero'', or ''custom''.');
end
end


function frameMask = localNormalizeFrameMask(frameMask, numFrame, numSat)
%LOCALNORMALIZEFRAMEMASK Normalize optional satellite-frame valid mask.

if isscalar(frameMask)
  frameMask = repmat(logical(frameMask), numSat, numFrame);
  return;
end

if ~islogical(frameMask)
  frameMask = logical(frameMask);
end

if isequal(size(frameMask), [1, numFrame])
  frameMask = repmat(frameMask, numSat, 1);
  return;
end

if isequal(size(frameMask), [numSat, numFrame])
  return;
end

error('buildFramePathGain:InvalidFrameMaskSize', ...
  ['pathOpt.frameMask must be a scalar, a 1xNf logical vector, or an ', ...
   'NsxNf logical matrix.']);
end


function localValidateAmpSpec(ampSpec, numFrame, numSat, numUser)
%LOCALVALIDATEAMPSPEC Validate amplitude specification.

if iscell(ampSpec)
  if numel(ampSpec) ~= numFrame
    error('buildFramePathGain:InvalidAmpCellLength', ...
      'Cell ampSpec must contain exactly numFrame entries.');
  end

  for iFrame = 1:numFrame
    localExpandAmpSpec(ampSpec{iFrame}, numSat, numUser, sprintf('ampSpec{%d}', iFrame));
  end
  return;
end

localExpandAmpSpec(ampSpec, numSat, numUser, 'ampSpec');
end


function localValidatePhaseSpec(phaseSpec, numFrame, numSat, numUser, phaseMode)
%LOCALVALIDATEPHASESPEC Validate phase specification when provided.

if isempty(phaseSpec)
  if strcmp(phaseMode, 'custom')
    error('buildFramePathGain:MissingPhaseSpec', ...
      'pathOpt.phaseSpec must be provided when phaseMode is ''custom''.');
  end
  return;
end

if iscell(phaseSpec)
  if numel(phaseSpec) ~= numFrame
    error('buildFramePathGain:InvalidPhaseCellLength', ...
      'Cell phaseSpec must contain exactly numFrame entries.');
  end

  for iFrame = 1:numFrame
    localExpandPhaseSpec(phaseSpec{iFrame}, numSat, numUser, sprintf('phaseSpec{%d}', iFrame));
  end
  return;
end

localExpandPhaseSpec(phaseSpec, numSat, numUser, 'phaseSpec');
end


function ampCell = localBuildAmpCell(ampSpec, ampJitterStdDb, numFrame, numSat, numUser)
%LOCALBUILDAMPCELL Build one amplitude matrix per frame.

ampCell = cell(1, numFrame);

if iscell(ampSpec)
  for iFrame = 1:numFrame
    ampCell{iFrame} = localExpandAmpSpec(ampSpec{iFrame}, numSat, numUser, ...
      sprintf('ampSpec{%d}', iFrame));
  end
else
  ampBase = localExpandAmpSpec(ampSpec, numSat, numUser, 'ampSpec');
  for iFrame = 1:numFrame
    ampCell{iFrame} = ampBase;
  end
end

if ampJitterStdDb <= 0
  return;
end

jitterScale = ampJitterStdDb / 20;
for iFrame = 1:numFrame
  ampCell{iFrame} = ampCell{iFrame} .* 10 .^ (jitterScale * randn(numSat, numUser));
end
end


function phaseCell = localBuildPhaseCell(phaseMode, phaseSpec, numFrame, numSat, numUser)
%LOCALBUILDPHASECELL Build one phase matrix per frame.

phaseCell = cell(1, numFrame);

switch phaseMode
  case 'zero'
    for iFrame = 1:numFrame
      phaseCell{iFrame} = zeros(numSat, numUser);
    end

  case 'shared'
    if isempty(phaseSpec)
      phaseBase = 2 * pi * rand(numSat, numUser);
    else
      phaseBase = localExpandPhaseSpec(phaseSpec, numSat, numUser, 'phaseSpec');
    end

    for iFrame = 1:numFrame
      phaseCell{iFrame} = phaseBase;
    end

  case 'independent'
    if isempty(phaseSpec)
      for iFrame = 1:numFrame
        phaseCell{iFrame} = 2 * pi * rand(numSat, numUser);
      end
    elseif iscell(phaseSpec)
      for iFrame = 1:numFrame
        phaseCell{iFrame} = localExpandPhaseSpec(phaseSpec{iFrame}, numSat, numUser, ...
          sprintf('phaseSpec{%d}', iFrame));
      end
    else
      phaseBase = localExpandPhaseSpec(phaseSpec, numSat, numUser, 'phaseSpec');
      for iFrame = 1:numFrame
        phaseCell{iFrame} = phaseBase;
      end
    end

  case 'custom'
    if iscell(phaseSpec)
      for iFrame = 1:numFrame
        phaseCell{iFrame} = localExpandPhaseSpec(phaseSpec{iFrame}, numSat, numUser, ...
          sprintf('phaseSpec{%d}', iFrame));
      end
    else
      phaseBase = localExpandPhaseSpec(phaseSpec, numSat, numUser, 'phaseSpec');
      for iFrame = 1:numFrame
        phaseCell{iFrame} = phaseBase;
      end
    end

  otherwise
    error('buildFramePathGain:UnsupportedPhaseMode', ...
      'Unsupported phase mode: %s.', phaseMode);
end
end


function ampMat = localExpandAmpSpec(ampSpec, numSat, numUser, varName)
%LOCALEXPANDAMPSPEC Expand one amplitude specification to Ns x Nu.

if isempty(ampSpec)
  ampMat = ones(numSat, numUser);
  return;
end

validateattributes(ampSpec, {'numeric'}, {'real', 'finite'}, mfilename, varName);

if isscalar(ampSpec)
  if ampSpec < 0
    error('buildFramePathGain:NegativeAmplitude', ...
      '%s must be nonnegative.', varName);
  end
  ampMat = repmat(ampSpec, numSat, numUser);
  return;
end

if isvector(ampSpec)
  ampVec = ampSpec(:);

  if numel(ampVec) == numUser
    if any(ampVec < 0)
      error('buildFramePathGain:NegativeAmplitude', ...
        '%s must be nonnegative.', varName);
    end
    ampMat = repmat(ampVec.', numSat, 1);
    return;
  end

  if numel(ampVec) == numSat
    if any(ampVec < 0)
      error('buildFramePathGain:NegativeAmplitude', ...
        '%s must be nonnegative.', varName);
    end
    ampMat = repmat(ampVec, 1, numUser);
    return;
  end

  error('buildFramePathGain:InvalidAmplitudeVectorLength', ...
    '%s vector length must equal numSat or numUser.', varName);
end

if isequal(size(ampSpec), [numSat, numUser])
  if any(ampSpec(:) < 0)
    error('buildFramePathGain:NegativeAmplitude', ...
      '%s must be nonnegative.', varName);
  end
  ampMat = ampSpec;
  return;
end

error('buildFramePathGain:InvalidAmplitudeSize', ...
  '%s must be scalar, vector, or an NsxNu matrix.', varName);
end


function phaseMat = localExpandPhaseSpec(phaseSpec, numSat, numUser, varName)
%LOCALEXPANDPHASESPEC Expand one phase specification to Ns x Nu.

if isempty(phaseSpec)
  phaseMat = zeros(numSat, numUser);
  return;
end

validateattributes(phaseSpec, {'numeric'}, {'real', 'finite'}, mfilename, varName);

if isscalar(phaseSpec)
  phaseMat = repmat(phaseSpec, numSat, numUser);
  return;
end

if isvector(phaseSpec)
  phaseVec = phaseSpec(:);

  if numel(phaseVec) == numUser
    phaseMat = repmat(phaseVec.', numSat, 1);
    return;
  end

  if numel(phaseVec) == numSat
    phaseMat = repmat(phaseVec, 1, numUser);
    return;
  end

  error('buildFramePathGain:InvalidPhaseVectorLength', ...
    '%s vector length must equal numSat or numUser.', varName);
end

if isequal(size(phaseSpec), [numSat, numUser])
  phaseMat = phaseSpec;
  return;
end

error('buildFramePathGain:InvalidPhaseSize', ...
  '%s must be scalar, vector, or an NsxNu matrix.', varName);
end
