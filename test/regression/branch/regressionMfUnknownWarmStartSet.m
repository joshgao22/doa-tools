function regressionMfUnknownWarmStartSet(varargin)
% Regression check for CP-U outer warm-start candidate construction.
% This regression focuses on one narrow contract only:
%   1) buildUnknownInitCandidateSet must expose the documented outer start
%      sources through startTag;
%   2) the outer candidate set must retain both fromCpK and fromStatic;
%   3) each outer candidate must carry finite initParam and a local DoA box;
%   4) outer-source identity must be preserved even when the two sources
%      collapse to the same numeric fdRate initializer.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();
fixture.dynMsKnownOpt.verbose = verbose;
fixture.dynMsUnknownOpt.verbose = verbose;

fprintf('Running regressionMfUnknownWarmStartSet ...\n');

initParamStaticMs = buildDynamicInitParamFromCase(fixture.bestStaticMsCase, true, fixture.truth.fdRateFit);
if isempty(initParamStaticMs)
  error('regressionMfUnknownWarmStartSet:MissingStaticInit', ...
    'Failed to build the static known-rate initializer.');
end

caseDynMsKnown = runDynamicDoaDopplerCase('MS-MF-CP-K', 'multi', ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.fdRange, fixture.fdRateRange, false, ...
  fixture.dynMsKnownOpt, true, [], initParamStaticMs);
initParamMsUnknownCpK = buildDynamicInitParamFromCase( ...
  caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamMsUnknownStatic = buildDynamicInitParamFromCase( ...
  fixture.bestStaticMsCase, false, fixture.truth.fdRateFit);

initCandUnknown = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, fixture.bestStaticMsCase, ...
  initParamMsUnknownCpK, initParamMsUnknownStatic, ...
  fixture.dynMsUnknownOpt.initDoaHalfWidth, fixture.staticMsOpt.initDoaHalfWidth);

if isempty(initCandUnknown) || ~iscell(initCandUnknown)
  error('regressionMfUnknownWarmStartSet:MissingInitCandidateSet', ...
    'buildUnknownInitCandidateSet must return a non-empty cell candidate set.');
end

numCand = numel(initCandUnknown);
startTagList = strings(numCand, 1);
fdRateInit = nan(numCand, 1);
initLen = nan(numCand, 1);
halfWidthNorm = nan(numCand, 1);

for iCand = 1:numCand
  cand = initCandUnknown{iCand};
  if ~isstruct(cand)
    error('regressionMfUnknownWarmStartSet:InvalidCandidateType', ...
      'Each CP-U outer warm-start candidate must be a struct.');
  end
  if ~isfield(cand, 'startTag') || isempty(cand.startTag)
    error('regressionMfUnknownWarmStartSet:MissingStartTag', ...
      'Each CP-U outer warm-start candidate must expose startTag.');
  end
  if ~isfield(cand, 'initParam') || isempty(cand.initParam)
    error('regressionMfUnknownWarmStartSet:MissingInitParam', ...
      'Each CP-U outer warm-start candidate must carry initParam.');
  end
  if ~isfield(cand, 'initDoaParam') || isempty(cand.initDoaParam)
    error('regressionMfUnknownWarmStartSet:MissingInitDoaParam', ...
      'Each CP-U outer warm-start candidate must carry initDoaParam.');
  end
  if ~isfield(cand, 'initDoaHalfWidth') || isempty(cand.initDoaHalfWidth)
    error('regressionMfUnknownWarmStartSet:MissingInitDoaHalfWidth', ...
      'Each CP-U outer warm-start candidate must carry initDoaHalfWidth.');
  end

  startTagList(iCand) = string(cand.startTag);
  initParamUse = reshape(cand.initParam, [], 1);
  fdRateInit(iCand) = initParamUse(end);
  initLen(iCand) = numel(initParamUse);
  halfWidthNorm(iCand) = norm(reshape(cand.initDoaHalfWidth, [], 1));

  if ~all(isfinite(initParamUse))
    error('regressionMfUnknownWarmStartSet:NonFiniteInitParam', ...
      'Each CP-U outer warm-start initParam must be finite.');
  end
  if ~all(isfinite(cand.initDoaParam(:)))
    error('regressionMfUnknownWarmStartSet:NonFiniteInitDoaParam', ...
      'Each CP-U outer warm-start initDoaParam must be finite.');
  end
  if ~all(isfinite(cand.initDoaHalfWidth(:))) || any(cand.initDoaHalfWidth(:) <= 0)
    error('regressionMfUnknownWarmStartSet:InvalidInitDoaHalfWidth', ...
      'Each CP-U outer warm-start initDoaHalfWidth must be finite and positive.');
  end
end

requiredStartTag = ["fromCpK"; "fromStatic"];
if numCand ~= 2 || ~all(ismember(requiredStartTag, startTagList))
  error('regressionMfUnknownWarmStartSet:UnexpectedStartTagSet', ...
    'CP-U outer warm-start candidates must contain exactly fromCpK and fromStatic.');
end

if any(initLen < 4)
  error('regressionMfUnknownWarmStartSet:ShortInitParam', ...
    'CP-U outer warm-start initParam must include DoA, fdRef, and fdRate.');
end
isCollapsedFdRateInit = range(fdRateInit) <= 0;

fprintf('  outer start tags        : %s\n', localFormatStringRow(startTagList));
fprintf('  outer fdRate init (Hz/s): %s\n', localFormatNumericRow(fdRateInit));
fprintf('  initParam length        : %s\n', localFormatNumericRow(initLen));
fprintf('  DoA half-width norm     : %s\n', localFormatNumericRow(halfWidthNorm));
fprintf('  collapsed fdRate init   : %d\n', isCollapsedFdRateInit);
fprintf('PASS: regressionMfUnknownWarmStartSet\n');

end

function textOut = localFormatStringRow(valueVec)
%LOCALFORMATSTRINGROW Format one string row for compact logging.

valueVec = reshape(string(valueVec), 1, []);
textOut = ['[', strjoin(cellstr(valueVec), ', '), ']'];
end


function textOut = localFormatNumericRow(valueVec)
%LOCALFORMATNUMERICROW Format one numeric row for compact logging.

valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6f', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact multi-sat dynamic fixture.

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
gridSize = [41, 41];
searchRange = [32.78, 42.78; 31.59, 41.59];
weightSweepAlpha = [0; 0.25; 0.5; 1];

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
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci);
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = sceneRef.satIdx(otherSatIdxLocal);

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1);
truth.latlonTrueDeg = usrLla(1:2, 1);
truth.fdRefTrueHz = truth.fdRefFit;
truth.fdRateTrueHzPerSec = truth.fdRateFit;
truth.refSatIdxGlobal = sceneRef.satIdx(refSatIdxLocal);

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

pathGainCell = repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame);
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots( ...
  sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, 0, pathGainCell, snapOpt);
rxSigRef = rxSigCell{sceneSeq.refFrameIdx};

sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(rxSigCell, refSatIdxLocal, 'multiFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeqRefOnly));

sceneOtherOnly = selectSatScene(sceneRef, otherSatIdxLocal);
rxSigOtherOnly = selectRxSigBySat(rxSigRef, otherSatIdxLocal, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
  gridSize, searchRange, E);

viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeq));

fdRange = localExpandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = localExpandRangeToTruth([-1e4, 0], ...
  [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);

doaOnlyOpt = struct();
doaOnlyOpt.useLogObjective = true;

staticBaseOpt = struct();
staticBaseOpt.useLogObjective = true;

dynBaseOpt = struct();
dynBaseOpt.useLogObjective = true;
dynBaseOpt.initFdCount = 81;
dynBaseOpt.useAccessMask = false;
dynBaseOpt.phaseMode = 'continuous';
dynBaseOpt.steeringMode = 'framewise';
dynBaseOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx;
dynBaseOpt.debugEnable = true;
dynBaseOpt.debugStoreEvalTrace = false;
dynBaseOpt.debugMaxEvalTrace = 120;

caseBundle = buildDoaDopplerStaticTransitionBundle( ...
  viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, ...
  carrierFreq, waveInfo.sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, false, doaOnlyOpt, staticBaseOpt, ...
  weightSweepAlpha, [0.01; 0.01]);

bestStaticMsCase = caseBundle.bestStaticMsCase;
staticMsOpt = caseBundle.staticMsOpt;

dynMsKnownOpt = dynBaseOpt;
dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsKnownOpt.initDoaHalfWidth = [0.003; 0.003];
dynMsKnownOpt.enableFdAliasUnwrap = true;

dynMsUnknownOpt = dynBaseOpt;
dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:);
dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002];
dynMsUnknownOpt.enableFdAliasUnwrap = true;

fixture = struct();
fixture.truth = truth;
fixture.truthSummary = struct('fdRefTrueHz', truth.fdRefFit, 'fdRateTrueHzPerSec', truth.fdRateFit);
fixture.viewMs = viewMs;
fixture.bestStaticMsCase = bestStaticMsCase;
fixture.staticMsOpt = staticMsOpt;
fixture.dynMsKnownOpt = dynMsKnownOpt;
fixture.dynMsUnknownOpt = dynMsUnknownOpt;
fixture.pilotWave = pilotWave;
fixture.carrierFreq = carrierFreq;
fixture.sampleRate = waveInfo.sampleRate;
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
