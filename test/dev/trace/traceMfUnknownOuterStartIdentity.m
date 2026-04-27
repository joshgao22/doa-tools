% Diagnostic trace for CP-U outer warm-start identity.
% This script checks three concrete trace questions only:
%   1) the outer fromCpK and fromStatic warm starts must remain distinct
%      before entering the CP-U multi-start wrapper;
%   2) if they collapse, the script should report whether the collapse has
%      already happened in the upstream source cases, in the packed
%      buildDynamicInitParamFromCase output, or in the outer candidate-set
%      object returned by buildUnknownInitCandidateSet;
%   3) the printed trace should be detailed enough to guide the next code
%      modification without rerunning the whole pipeline blindly.
clear(); close all; clc;
fixture = localBuildRegressionFixture();

fprintf('Running traceMfUnknownOuterStartIdentity ...\n');

initParamStaticMs = buildDynamicInitParamFromCase(fixture.bestStaticMsCase, true, fixture.truth.fdRateFit);
if isempty(initParamStaticMs)
  error('traceMfUnknownOuterStartIdentity:MissingStaticInit', ...
    'Failed to build the static known-rate initializer.');
end

caseDynMsKnown = runDynamicDoaDopplerCase('MS-MF-CP-K', 'multi', ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.fdRange, fixture.fdRateRange, false, ...
  fixture.dynMsKnownOpt, true, [], initParamStaticMs);
initParamUnknownCpK = buildDynamicInitParamFromCase( ...
  caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initParamUnknownStatic = buildDynamicInitParamFromCase( ...
  fixture.bestStaticMsCase, false, fixture.truth.fdRateFit);

if isempty(initParamUnknownCpK) || isempty(initParamUnknownStatic)
  error('traceMfUnknownOuterStartIdentity:MissingUnknownInit', ...
    'Failed to build one or both CP-U outer warm-start vectors.');
end
if numel(initParamUnknownCpK) ~= numel(initParamUnknownStatic)
  error('traceMfUnknownOuterStartIdentity:InitSizeMismatch', ...
    'The fromCpK and fromStatic initParam vectors must have the same length.');
end
if numel(initParamUnknownCpK) < 4
  error('traceMfUnknownOuterStartIdentity:UnexpectedInitLength', ...
    'Expected CP-U initParam to contain DoA, fdRef, and fdRate.');
end

srcCpK = localBuildSourceState(caseDynMsKnown.estResult, caseDynMsKnown.estResult.fdRateEst);
srcStatic = localBuildSourceState(fixture.bestStaticMsCase.estResult, fixture.truth.fdRateFit);
startCpK = localUnpackUnknownInit(initParamUnknownCpK);
startStatic = localUnpackUnknownInit(initParamUnknownStatic);

srcGap = localComputeStateGap(srcCpK, srcStatic, []);
packedGap = localComputeStateGap(startCpK, startStatic, [initParamUnknownCpK(:), initParamUnknownStatic(:)]);

srcCpKAngleErrDeg = calcLatlonAngleError(caseDynMsKnown.estResult.doaParamEst(:), startCpK.latlon);
srcStaticAngleErrDeg = calcLatlonAngleError(fixture.bestStaticMsCase.estResult.doaParamEst(:), startStatic.latlon);
srcCpKFdRefGapHz = abs(caseDynMsKnown.estResult.fdRefEst - startCpK.fdRef);
srcStaticFdRefGapHz = abs(fixture.bestStaticMsCase.estResult.fdRefEst - startStatic.fdRef);
srcCpKFdRateGapHzPerSec = abs(caseDynMsKnown.estResult.fdRateEst - startCpK.fdRate);
srcStaticFdRateGapHzPerSec = abs(fixture.truth.fdRateFit - startStatic.fdRate);

if srcCpKAngleErrDeg > 1e-9 || srcStaticAngleErrDeg > 1e-9 || ...
    srcCpKFdRefGapHz > 1e-9 || srcStaticFdRefGapHz > 1e-9 || ...
    srcCpKFdRateGapHzPerSec > 1e-9 || srcStaticFdRateGapHzPerSec > 1e-9
  error('traceMfUnknownOuterStartIdentity:SourceMismatch', ...
    ['At least one outer warm-start vector is inconsistent with its ', ...
     'documented source case.']);
end

initCandUnknown = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, fixture.bestStaticMsCase, ...
  initParamUnknownCpK, initParamUnknownStatic, ...
  fixture.dynMsUnknownOpt.initDoaHalfWidth, fixture.staticMsOpt.initDoaHalfWidth);
outerTrace = localDescribeOuterCandidateSet(initCandUnknown, numel(initParamUnknownCpK));
traceStage = localInferCollapseStage(srcGap, packedGap, outerTrace);

fprintf('  source fromCpK latlon (deg)      : %s\n', localFormatNumericRow(srcCpK.latlon));
fprintf('  source fromStatic latlon (deg)   : %s\n', localFormatNumericRow(srcStatic.latlon));
fprintf('  source fromCpK [fdRef,fdRate]    : %s\n', localFormatNumericRow([srcCpK.fdRef, srcCpK.fdRate]));
fprintf('  source fromStatic [fdRef,fdRate] : %s\n', localFormatNumericRow([srcStatic.fdRef, srcStatic.fdRate]));
fprintf('  source DoA gap (deg)             : %.6f\n', srcGap.doaGapDeg);
fprintf('  source fdRef gap (Hz)            : %.6f\n', srcGap.fdRefGapHz);
fprintf('  source fdRate gap (Hz/s)         : %.6f\n', srcGap.fdRateGapHzPerSec);

fprintf('  packed fromCpK latlon (deg)      : %s\n', localFormatNumericRow(startCpK.latlon));
fprintf('  packed fromStatic latlon (deg)   : %s\n', localFormatNumericRow(startStatic.latlon));
fprintf('  packed fromCpK [fdRef,fdRate]    : %s\n', localFormatNumericRow([startCpK.fdRef, startCpK.fdRate]));
fprintf('  packed fromStatic [fdRef,fdRate] : %s\n', localFormatNumericRow([startStatic.fdRef, startStatic.fdRate]));
fprintf('  packed DoA gap (deg)             : %.6f\n', packedGap.doaGapDeg);
fprintf('  packed fdRef gap (Hz)            : %.6f\n', packedGap.fdRefGapHz);
fprintf('  packed fdRate gap (Hz/s)         : %.6f\n', packedGap.fdRateGapHzPerSec);
fprintf('  packed initParam gap norm        : %.6f\n', packedGap.paramGapNorm);

fprintf('  outer candidate object type      : %s\n', class(initCandUnknown));
if outerTrace.hasCandidates
  fprintf('  outer candidate tags             : %s\n', localFormatStringRow(outerTrace.tagList));
  fprintf('  outer candidate count            : %d\n', outerTrace.count);
  fprintf('  outer candidate pair gap norm    : %.6f\n', outerTrace.bestPairGapNorm);
  fprintf('  outer candidate pair DoA gap     : %.6f\n', outerTrace.bestPairDoaGapDeg);
  fprintf('  outer candidate pair fdRef gap   : %.6f\n', outerTrace.bestPairFdRefGapHz);
  fprintf('  outer candidate pair fdRate gap  : %.6f\n', outerTrace.bestPairFdRateGapHzPerSec);
else
  fprintf('  outer candidate count            : 0 (no parseable candidate initParam found)\n');
end
fprintf('  inferred collapse stage          : %s\n', traceStage);

collapseTol = 1e-8 * max(1, norm([initParamUnknownCpK(:); initParamUnknownStatic(:)]));
if packedGap.paramGapNorm <= collapseTol && packedGap.doaGapDeg <= 1e-9 && ...
    packedGap.fdRefGapHz <= 1e-9 && packedGap.fdRateGapHzPerSec <= 1e-9
  error('traceMfUnknownOuterStartIdentity:OuterStartsCollapsed', ...
    ['The fromCpK and fromStatic outer warm starts collapsed before CP-U ', ...
     'solving. Inferred stage: %s.'], traceStage);
end

fprintf('PASS: traceMfUnknownOuterStartIdentity\n');


function state = localBuildSourceState(estResult, fdRateValue)
%LOCALBUILDSOURCESTATE Convert one upstream source case to a compact state.

state = struct();
state.latlon = reshape(estResult.doaParamEst, 2, 1);
state.fdRef = estResult.fdRefEst;
state.fdRate = fdRateValue;
end


function stateGap = localComputeStateGap(stateA, stateB, initPair)
%LOCALCOMPUTESTATEGAP Compute one compact gap summary between two states.

stateGap = struct();
stateGap.doaGapDeg = calcLatlonAngleError(stateA.latlon, stateB.latlon);
stateGap.fdRefGapHz = abs(stateA.fdRef - stateB.fdRef);
stateGap.fdRateGapHzPerSec = abs(stateA.fdRate - stateB.fdRate);
stateGap.paramGapNorm = NaN;
if ~isempty(initPair)
  stateGap.paramGapNorm = norm(initPair(:, 1) - initPair(:, 2));
end
end


function initState = localUnpackUnknownInit(initParam)
%LOCALUNPACKUNKNOWNINIT Unpack one CP-U initParam vector.

initParam = initParam(:);
initState = struct();
initState.latlon = initParam(1:2);
initState.fdRef = initParam(end - 1);
initState.fdRate = initParam(end);
end


function outerTrace = localDescribeOuterCandidateSet(initCandUnknown, initLen)
%LOCALDESCRIBEOUTERCANDIDATESET Best-effort trace of the outer candidate object.

outerTrace = struct();
outerTrace.hasCandidates = false;
outerTrace.count = 0;
outerTrace.tagList = strings(1, 0);
outerTrace.initParamMat = zeros(initLen, 0);
outerTrace.bestPairGapNorm = NaN;
outerTrace.bestPairDoaGapDeg = NaN;
outerTrace.bestPairFdRefGapHz = NaN;
outerTrace.bestPairFdRateGapHzPerSec = NaN;

[tagList, initParamMat] = localExtractCandidateSet(initCandUnknown, initLen);
if isempty(initParamMat)
  return;
end

outerTrace.hasCandidates = true;
outerTrace.count = size(initParamMat, 2);
outerTrace.tagList = tagList;
outerTrace.initParamMat = initParamMat;

if size(initParamMat, 2) >= 2
  gapNorm = inf;
  gapDoa = inf;
  gapFdRef = inf;
  gapFdRate = inf;
  for iCand = 1:(size(initParamMat, 2) - 1)
    for jCand = (iCand + 1):size(initParamMat, 2)
      gapNormUse = norm(initParamMat(:, iCand) - initParamMat(:, jCand));
      stateI = localUnpackUnknownInit(initParamMat(:, iCand));
      stateJ = localUnpackUnknownInit(initParamMat(:, jCand));
      gapDoaUse = calcLatlonAngleError(stateI.latlon, stateJ.latlon);
      gapFdRefUse = abs(stateI.fdRef - stateJ.fdRef);
      gapFdRateUse = abs(stateI.fdRate - stateJ.fdRate);
      if gapNormUse < gapNorm
        gapNorm = gapNormUse;
        gapDoa = gapDoaUse;
        gapFdRef = gapFdRefUse;
        gapFdRate = gapFdRateUse;
      end
    end
  end
  outerTrace.bestPairGapNorm = gapNorm;
  outerTrace.bestPairDoaGapDeg = gapDoa;
  outerTrace.bestPairFdRefGapHz = gapFdRef;
  outerTrace.bestPairFdRateGapHzPerSec = gapFdRate;
end
end


function [tagList, initParamMat] = localExtractCandidateSet(candObj, initLen)
%LOCALEXTRACTCANDIDATESET Best-effort extraction of candidate initParam vectors.

tagList = strings(1, 0);
initParamMat = zeros(initLen, 0);

if isempty(candObj)
  return;
end
if isnumeric(candObj)
  [tagList, initParamMat] = localExtractNumericCandidateSet(candObj, initLen, "cand");
  return;
end
if iscell(candObj)
  for iCell = 1:numel(candObj)
    [tagUse, initUse] = localExtractCandidateEntry(candObj{iCell}, initLen, sprintf('cell%d', iCell));
    if ~isempty(initUse)
      tagList(end + 1) = string(tagUse); %#ok<AGROW>
      initParamMat(:, end + 1) = initUse(:); %#ok<AGROW>
    end
  end
  return;
end
if isstruct(candObj)
  if numel(candObj) > 1
    for iEntry = 1:numel(candObj)
      [tagUse, initUse] = localExtractCandidateEntry(candObj(iEntry), initLen, sprintf('struct%d', iEntry));
      if ~isempty(initUse)
        tagList(end + 1) = string(tagUse); %#ok<AGROW>
        initParamMat(:, end + 1) = initUse(:); %#ok<AGROW>
      end
    end
    return;
  end
  if isfield(candObj, 'startTagList')
    srcTag = reshape(string(candObj.startTagList), 1, []);
  elseif isfield(candObj, 'startTag')
    srcTag = reshape(string(candObj.startTag), 1, []);
  else
    srcTag = strings(1, 0);
  end
  fieldOrder = {'initParamList', 'initParamMat', 'initParam', 'candidateInitParam'};
  for iField = 1:numel(fieldOrder)
    fieldUse = fieldOrder{iField};
    if isfield(candObj, fieldUse)
      valueUse = candObj.(fieldUse);
      if isnumeric(valueUse)
        [tagList, initParamMat] = localExtractNumericCandidateSet(valueUse, initLen, fieldUse);
      elseif iscell(valueUse)
        for iCell = 1:numel(valueUse)
          tagFallback = sprintf('%s%d', fieldUse, iCell);
          if iCell <= numel(srcTag)
            tagFallback = srcTag(iCell);
          end
          [tagUse, initUse] = localExtractCandidateEntry(valueUse{iCell}, initLen, tagFallback);
          if ~isempty(initUse)
            tagList(end + 1) = string(tagUse); %#ok<AGROW>
            initParamMat(:, end + 1) = initUse(:); %#ok<AGROW>
          end
        end
      end
      break;
    end
  end
  if isempty(initParamMat)
    [tagUse, initUse] = localExtractCandidateEntry(candObj, initLen, 'struct1');
    if ~isempty(initUse)
      tagList = string(tagUse);
      initParamMat = initUse(:);
    end
  elseif ~isempty(srcTag) && numel(srcTag) == size(initParamMat, 2)
    tagList = srcTag;
  end
end
end


function [tagUse, initUse] = localExtractCandidateEntry(entryObj, initLen, tagFallback)
%LOCALEXTRACTCANDIDATEENTRY Extract one candidate entry from a mixed object.

tagUse = string(tagFallback);
initUse = [];
if isempty(entryObj)
  return;
end
if isnumeric(entryObj)
  vecUse = entryObj(:);
  if numel(vecUse) == initLen
    initUse = vecUse;
  end
  return;
end
if isstruct(entryObj)
  if isfield(entryObj, 'tag')
    tagUse = string(entryObj.startTag);
  elseif isfield(entryObj, 'name')
    tagUse = string(entryObj.name);
  end
  if isfield(entryObj, 'initParam')
    vecUse = entryObj.initParam(:);
  elseif isfield(entryObj, 'param')
    vecUse = entryObj.param(:);
  else
    vecUse = [];
  end
  if numel(vecUse) == initLen
    initUse = vecUse;
  end
end
end


function [tagList, initParamMat] = localExtractNumericCandidateSet(valueUse, initLen, tagPrefix)
%LOCALEXTRACTNUMERICCANDIDATESET Extract columns or rows from one numeric object.

tagList = strings(1, 0);
initParamMat = zeros(initLen, 0);
if isempty(valueUse)
  return;
end
if isvector(valueUse)
  vecUse = valueUse(:);
  if numel(vecUse) == initLen
    tagList = string(tagPrefix);
    initParamMat = vecUse;
  end
  return;
end
if size(valueUse, 1) == initLen
  initParamMat = valueUse;
elseif size(valueUse, 2) == initLen
  initParamMat = valueUse.';
else
  return;
end
for iCand = 1:size(initParamMat, 2)
  tagList(iCand) = sprintf('%s%d', tagPrefix, iCand); %#ok<AGROW>
end
end


function stageText = localInferCollapseStage(srcGap, packedGap, outerTrace)
%LOCALINFERCOLLAPSESTAGE Infer the earliest stage where the two starts collapse.

stageText = "distinct";
stateTol = 1e-9;
paramTol = 1e-8 * max(1, packedGap.paramGapNorm);
sourceCollapsed = srcGap.doaGapDeg <= stateTol && ...
  srcGap.fdRefGapHz <= stateTol && srcGap.fdRateGapHzPerSec <= stateTol;
packedCollapsed = packedGap.doaGapDeg <= stateTol && ...
  packedGap.fdRefGapHz <= stateTol && packedGap.fdRateGapHzPerSec <= stateTol && ...
  packedGap.paramGapNorm <= max(paramTol, 1e-12);
outerCollapsed = false;
if outerTrace.hasCandidates && outerTrace.count >= 2
  outerCollapsed = outerTrace.bestPairDoaGapDeg <= stateTol && ...
    outerTrace.bestPairFdRefGapHz <= stateTol && ...
    outerTrace.bestPairFdRateGapHzPerSec <= stateTol && ...
    outerTrace.bestPairGapNorm <= max(paramTol, 1e-12);
end

if sourceCollapsed
  stageText = "upstreamSourceCase";
elseif packedCollapsed
  stageText = "buildDynamicInitParamFromCase";
elseif outerCollapsed
  stageText = "buildUnknownInitCandidateSet";
elseif outerTrace.hasCandidates
  stageText = "survivesIntoOuterCandidateSet";
else
  stageText = "packedInitOnly_noCandidateTrace";
end
end


function textOut = localFormatNumericRow(valueVec)
%LOCALFORMATNUMERICROW Format one numeric row for compact logging.

valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6f', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end


function textOut = localFormatStringRow(valueVec)
%LOCALFORMATSTRINGROW Format one string row for compact logging.

valueVec = reshape(string(valueVec), 1, []);
textOut = ['[', strjoin(cellstr(valueVec), ', '), ']'];
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
