% Regression check for CP-U warm-start candidate construction.
% This script focuses on one narrow contract only:
%   1) the unknown-rate branch must expose the local warm-start tags through
%      optimInfo.candidateVariant;
%   2) the parsed fdRate anchors must match the expected local candidate set
%      around the current unknown-rate seed;
%   3) the anchors must stay inside fdRateRange and preserve stable order.
clear(); close all; clc;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running regressionMfUnknownWarmStartSet ...\n');

mfOpt = struct();
mfOpt.useLogObjective = false;
mfOpt.phaseMode = 'continuous';
mfOpt.fdRateMode = 'unknown';
mfOpt.steeringMode = 'framewise';
mfOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
mfOpt.initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.05; -0.04];
mfOpt.initDoaHalfWidth = [0.10; 0.10];
mfOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-10);

[modelUnknown, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, mfOpt);
[initParamUnknown, ~] = buildDoaDopplerMfInit(modelUnknown, []);
[~, optimInfoUnknown] = solveDoaDopplerMfBranches(modelUnknown, initParamUnknown, false);

if ~isfield(optimInfoUnknown, 'candidateVariant') || isempty(optimInfoUnknown.candidateVariant)
  error('regressionMfUnknownWarmStartSet:MissingCandidateVariant', ...
    'CP-U did not expose any continuation candidateVariant tags.');
end
if ~isfield(optimInfoUnknown, 'candidateObjective') || isempty(optimInfoUnknown.candidateObjective)
  error('regressionMfUnknownWarmStartSet:MissingCandidateObjective', ...
    'CP-U did not expose candidateObjective for the warm-start set.');
end

candidateVariant = reshape(string(optimInfoUnknown.candidateVariant), 1, []);
candidateObjective = reshape(optimInfoUnknown.candidateObjective, 1, []);
if numel(candidateVariant) ~= numel(candidateObjective)
  error('regressionMfUnknownWarmStartSet:CandidateSizeMismatch', ...
    'candidateVariant and candidateObjective must have the same length.');
end
if ~all(contains(candidateVariant, "unknownWarmLocal_fdRate"))
  error('regressionMfUnknownWarmStartSet:UnexpectedCandidateTag', ...
    'All CP-U continuation candidates must expose unknownWarmLocal_fdRate tags.');
end

fdRateSeed = initParamUnknown(end);
[fdRateCandActual, tagResolution] = localParseWarmRateFromTag(candidateVariant);
fdRateCandExpect = localBuildExpectedWarmSet(modelUnknown.fdRateRange, fdRateSeed);
fdRateCandExpectTag = localRoundToResolution(fdRateCandExpect, tagResolution);

if numel(fdRateCandActual) ~= numel(fdRateCandExpectTag)
  error('regressionMfUnknownWarmStartSet:CandidateCountMismatch', ...
    ['The parsed CP-U warm-start candidate count (%d) does not match the ', ...
     'expected local set (%d).'], numel(fdRateCandActual), numel(fdRateCandExpectTag));
end

rateDiffMax = max(abs(fdRateCandActual - fdRateCandExpectTag));
rateTol = max(5e-4, 0.51 * tagResolution);
if rateDiffMax > rateTol
  error('regressionMfUnknownWarmStartSet:CandidateValueMismatch', ...
    ['The CP-U warm-start fdRate anchors do not match the expected local ', ...
     'candidate set after applying tag precision (max diff %.3e, tol %.3e).'], ...
     rateDiffMax, rateTol);
end

if any(fdRateCandActual < modelUnknown.fdRateRange(1) - 1e-12) || ...
    any(fdRateCandActual > modelUnknown.fdRateRange(2) + 1e-12)
  error('regressionMfUnknownWarmStartSet:CandidateOutOfRange', ...
    'At least one CP-U warm-start fdRate anchor left model.fdRateRange.');
end
if any(~isfinite(candidateObjective))
  error('regressionMfUnknownWarmStartSet:NonFiniteCandidateObjective', ...
    'At least one CP-U warm-start candidateObjective is non-finite.');
end

fprintf('  fdRate seed (Hz/s)      : %.6f\n', fdRateSeed);
fprintf('  candidate count         : %d\n', numel(fdRateCandActual));
fprintf('  candidate fdRate (Hz/s) : %s\n', localFormatNumericRow(fdRateCandActual));
fprintf('  max fdRate diff         : %.3e\n', rateDiffMax);
fprintf('PASS: regressionMfUnknownWarmStartSet\n');


function fdRateCand = localBuildExpectedWarmSet(fdRateRange, fdRateSeed)
%LOCALBUILDEXPECTEDWARMSET Rebuild the CP-U local warm-start anchors.

fdMin = fdRateRange(1);
fdMax = fdRateRange(2);
rangeSpan = fdMax - fdMin;
if ~(isfinite(rangeSpan) && rangeSpan >= 0)
  error('regressionMfUnknownWarmStartSet:InvalidFdRateRange', ...
    'fdRateRange must be finite and ordered.');
end

if rangeSpan == 0
  fdRateCand = fdMin;
  return;
end
if ~isfinite(fdRateSeed)
  fdRateSeed = 0.5 * (fdMin + fdMax);
end

localHalfWidth = min(0.08 * rangeSpan, max(25, 0.03 * max(abs(fdRateSeed), rangeSpan)));
fdRateCand = [fdRateSeed - localHalfWidth, fdRateSeed, fdRateSeed + localHalfWidth];
fdRateCand = min(max(fdRateCand, fdMin), fdMax);
fdRateCand = localUniqueStableTol(fdRateCand, 1e-9);
end


function [fdRateCand, resolution] = localParseWarmRateFromTag(candidateVariant)
%LOCALPARSEWARMRATEFROMTAG Parse fdRate anchors from solveVariant tags.

prefix = "unknownWarmLocal_fdRate";
fdRateCand = nan(size(candidateVariant));
resolution = inf;
for iCand = 1:numel(candidateVariant)
  tagUse = string(candidateVariant(iCand));
  if ~startsWith(tagUse, prefix)
    error('regressionMfUnknownWarmStartSet:InvalidVariantTag', ...
      'Unexpected CP-U continuation tag: %s.', tagUse);
  end
  valueText = char(extractAfter(tagUse, prefix));
  fdRateCand(iCand) = str2double(valueText);
  if ~isfinite(fdRateCand(iCand))
    error('regressionMfUnknownWarmStartSet:InvalidVariantValue', ...
      'Unable to parse a finite fdRate from tag %s.', tagUse);
  end
  dotIdx = find(valueText == '.', 1, 'last');
  if isempty(dotIdx)
    resolution = min(resolution, 1);
  else
    numDecimal = numel(valueText) - dotIdx;
    resolution = min(resolution, 10^(-numDecimal));
  end
end
if ~isfinite(resolution)
  resolution = 1e-3;
end
end


function valueVec = localRoundToResolution(valueVec, resolution)
%LOCALROUNDTORESOLUTION Round values to the printed tag precision.

if ~(isfinite(resolution) && resolution > 0)
  return;
end
valueVec = resolution * round(valueVec ./ resolution);
end


function uniqueVec = localUniqueStableTol(valVec, tol)
%LOCALUNIQUESTABLETOL Stable unique with tolerance.

valVec = reshape(valVec, 1, []);
uniqueVec = zeros(1, 0);
for iVal = 1:numel(valVec)
  if isempty(uniqueVec) || all(abs(valVec(iVal) - uniqueVec) > tol)
    uniqueVec(end + 1) = valVec(iVal); %#ok<AGROW>
  end
end
end


function textOut = localFormatNumericRow(valueVec)
%LOCALFORMATNUMERICROW Format one numeric row for compact logging.

valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6f', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one noiseless dynamic multi-sat case.

rng(253);

numUsr = 1;
numFrame = 7;
refFrameIdx = ceil(numFrame / 2);
frameIntvlSec = 1 / 750;

sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
numSym = 32;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
E = referenceEllipsoid('sphere');

utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);
usrLla = [37.78; 36.59; 0];

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1);
truth.latlonTrueDeg = usrLla(1:2, 1);

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);

snapOpt = struct();
snapOpt.spatial.model = 'dynamic';
snapOpt.spatial.refFrameIdx = sceneSeq.refFrameIdx;
snapOpt.phase.timeModel = 'global';
snapOpt.phase.frameModel = 'shared';
snapOpt.phase.sharedPhase = 2 * pi * rand(sceneSeq.numSat, sceneSeq.numUser);
snapOpt.wave.delayModel = 'phaseOnly';
snapOpt.wave.timeRef = 'zero';
snapOpt.wave.carrierPhaseModel = 'none';
snapOpt.precomp.linkParamCell = linkParamCell;

[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, 0, repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame), snapOpt);

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, [41, 41], searchRange, E, ...
  struct('sceneSeq', sceneSeq));

fdRange = localExpandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = localExpandRangeToTruth([-1e4, 0], ...
  [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.truth = truth;
fixture.viewMs = viewMs;
fixture.rxSigCell = rxSigCell;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
end


function value = localExpandRangeToTruth(defaultRange, truthVals, fracPad, absPad)
%LOCALEXPANDRANGETOTRUTH Expand a nominal range so that truth is included.

truthVals = truthVals(isfinite(truthVals));
if isempty(truthVals)
  value = defaultRange;
  return;
end

minTruth = min(truthVals);
maxTruth = max(truthVals);
spanTruth = max(maxTruth - minTruth, eps);
pad = max(absPad, fracPad * spanTruth);
value = [min(defaultRange(1), minTruth - pad), ...
         max(defaultRange(2), maxTruth + pad)];
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
