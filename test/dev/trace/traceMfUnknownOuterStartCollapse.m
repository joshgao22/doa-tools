% Trace script for CP-U outer warm-start collapse with source provenance.
% This script does not enforce a hard regression contract. It prints the
% source-case states, the packed CP-U initParam vectors, the best-effort
% outer candidate object trace, and a provenance check that tells whether
% fromCpK really comes from the CP-K final case and fromStatic really comes
% from the best static case.
clear(); close all; clc;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running traceMfUnknownOuterStartCollapse ...\n');

initParamStaticMs = buildDynamicInitParamFromCase(fixture.bestStaticMsCase, true, fixture.truth.fdRateFit);
caseDynMsKnown = runDynamicDoaDopplerCase('MS-MF-CP-K', 'multi', ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.fdRange, fixture.fdRateRange, false, ...
  fixture.dynMsKnownOpt, true, [], initParamStaticMs);

srcCpK = localBuildSourceState(caseDynMsKnown.estResult, caseDynMsKnown.estResult.fdRateEst);
srcStatic = localBuildSourceState(fixture.bestStaticMsCase.estResult, fixture.truth.fdRateFit);
initCpK = buildDynamicInitParamFromCase(caseDynMsKnown, false, caseDynMsKnown.estResult.fdRateEst);
initStatic = buildDynamicInitParamFromCase(fixture.bestStaticMsCase, false, fixture.truth.fdRateFit);
packedCpK = localUnpackUnknownInit(initCpK);
packedStatic = localUnpackUnknownInit(initStatic);

initCandUnknown = buildUnknownInitCandidateSet( ...
  caseDynMsKnown, fixture.bestStaticMsCase, ...
  initCpK, initStatic, ...
  fixture.dynMsUnknownOpt.initDoaHalfWidth, fixture.staticMsOpt.initDoaHalfWidth);
outerTrace = localDescribeOuterCandidateSet(initCandUnknown, numel(initCpK));
traceStage = localInferCollapseStage(srcCpK, srcStatic, packedCpK, packedStatic, outerTrace);
provenance = localBuildProvenanceSummary(caseDynMsKnown, fixture.bestStaticMsCase, srcCpK, srcStatic, packedCpK, packedStatic);

fprintf('  source fromCpK latlon (deg)      : %s\n', localFormatNumericRow(srcCpK.latlon));
fprintf('  source fromStatic latlon (deg)   : %s\n', localFormatNumericRow(srcStatic.latlon));
fprintf('  source fromCpK [fdRef,fdRate]    : %s\n', localFormatNumericRow([srcCpK.fdRef, srcCpK.fdRate]));
fprintf('  source fromStatic [fdRef,fdRate] : %s\n', localFormatNumericRow([srcStatic.fdRef, srcStatic.fdRate]));
fprintf('  source DoA gap (deg)             : %.6f\n', calcLatlonAngleError(srcCpK.latlon, srcStatic.latlon));
fprintf('  source fdRef gap (Hz)            : %.6f\n', abs(srcCpK.fdRef - srcStatic.fdRef));
fprintf('  source fdRate gap (Hz/s)         : %.6f\n', abs(srcCpK.fdRate - srcStatic.fdRate));

fprintf('  packed fromCpK latlon (deg)      : %s\n', localFormatNumericRow(packedCpK.latlon));
fprintf('  packed fromStatic latlon (deg)   : %s\n', localFormatNumericRow(packedStatic.latlon));
fprintf('  packed fromCpK [fdRef,fdRate]    : %s\n', localFormatNumericRow([packedCpK.fdRef, packedCpK.fdRate]));
fprintf('  packed fromStatic [fdRef,fdRate] : %s\n', localFormatNumericRow([packedStatic.fdRef, packedStatic.fdRate]));
fprintf('  packed DoA gap (deg)             : %.6f\n', calcLatlonAngleError(packedCpK.latlon, packedStatic.latlon));
fprintf('  packed fdRef gap (Hz)            : %.6f\n', abs(packedCpK.fdRef - packedStatic.fdRef));
fprintf('  packed fdRate gap (Hz/s)         : %.6f\n', abs(packedCpK.fdRate - packedStatic.fdRate));
fprintf('  packed initParam gap norm        : %.6f\n', norm(initCpK(:) - initStatic(:)));

fprintf('  provenance fromCpK <- CP-K final : %d\n', provenance.fromCpKMatchesKnownFinal);
fprintf('  provenance fromStatic <- bestStatic : %d\n', provenance.fromStaticMatchesBestStatic);
fprintf('  CP-K final vs bestStatic DoA gap : %.6f\n', provenance.knownVsStaticDoaGapDeg);
fprintf('  CP-K final vs bestStatic fdRef gap (Hz) : %.6f\n', provenance.knownVsStaticFdRefGapHz);
fprintf('  CP-K final vs bestStatic fdRate gap (Hz/s) : %.6f\n', provenance.knownVsStaticFdRateGapHzPerSec);

fprintf('  outer candidate object type      : %s\n', class(initCandUnknown));
if outerTrace.hasCandidates
  fprintf('  outer candidate tags             : %s\n', localFormatStringRow(outerTrace.tagList));
  fprintf('  outer candidate count            : %d\n', outerTrace.count);
  fprintf('  outer candidate pair gap norm    : %.6f\n', outerTrace.bestPairGapNorm);
  fprintf('  outer candidate pair DoA gap     : %.6f\n', outerTrace.bestPairDoaGapDeg);
  fprintf('  outer candidate pair fdRef gap   : %.6f\n', outerTrace.bestPairFdRefGapHz);
  fprintf('  outer candidate pair fdRate gap  : %.6f\n', outerTrace.bestPairFdRateGapHzPerSec);
else
  fprintf('  outer candidate count            : 0 (candidate object not parseable)\n');
end
fprintf('  inferred collapse stage          : %s\n', traceStage);
fprintf('DONE: traceMfUnknownOuterStartCollapse\n');


function state = localBuildSourceState(estResult, fdRateValue)
state = struct();
state.latlon = reshape(estResult.doaParamEst, 2, 1);
state.fdRef = estResult.fdRefEst;
state.fdRate = fdRateValue;
end

function initState = localUnpackUnknownInit(initParam)
initParam = initParam(:);
initState = struct();
initState.latlon = initParam(1:2);
initState.fdRef = initParam(end - 1);
initState.fdRate = initParam(end);
end

function provenance = localBuildProvenanceSummary(caseDynMsKnown, bestStaticMsCase, srcCpK, srcStatic, packedCpK, packedStatic)
provenance = struct();
stateTol = 1e-9;
knownFinal = localBuildSourceState(caseDynMsKnown.estResult, caseDynMsKnown.estResult.fdRateEst);
bestStaticFdRate = localResolveBestStaticFdRate(bestStaticMsCase, srcStatic);
bestStatic = localBuildSourceState(bestStaticMsCase.estResult, bestStaticFdRate);
provenance.fromCpKMatchesKnownFinal = localStateGapNorm(srcCpK, knownFinal) <= stateTol && ...
  localStateGapNorm(packedCpK, knownFinal) <= stateTol;
provenance.fromStaticMatchesBestStatic = localStateGapNorm(srcStatic, bestStatic) <= stateTol && ...
  localStateGapNorm(packedStatic, bestStatic) <= stateTol;
provenance.knownVsStaticDoaGapDeg = calcLatlonAngleError(knownFinal.latlon, bestStatic.latlon);
provenance.knownVsStaticFdRefGapHz = abs(knownFinal.fdRef - bestStatic.fdRef);
provenance.knownVsStaticFdRateGapHzPerSec = abs(knownFinal.fdRate - bestStatic.fdRate);
provenance.bestStaticFdRateSource = localDescribeBestStaticFdRateSource(bestStaticMsCase);
end

function fdRateUse = localResolveBestStaticFdRate(bestStaticMsCase, srcStatic)
fdRateUse = srcStatic.fdRate;
if isstruct(bestStaticMsCase) && isfield(bestStaticMsCase, 'truth') && isstruct(bestStaticMsCase.truth)
  if isfield(bestStaticMsCase.truth, 'fdRateTrueHzPerSec') && isfinite(bestStaticMsCase.truth.fdRateTrueHzPerSec)
    fdRateUse = bestStaticMsCase.truth.fdRateTrueHzPerSec;
    return;
  end
  if isfield(bestStaticMsCase.truth, 'fdRateFit') && isfinite(bestStaticMsCase.truth.fdRateFit)
    fdRateUse = bestStaticMsCase.truth.fdRateFit;
    return;
  end
end
if isstruct(bestStaticMsCase) && isfield(bestStaticMsCase, 'truthSummary') && isstruct(bestStaticMsCase.truthSummary)
  if isfield(bestStaticMsCase.truthSummary, 'fdRateTrueHzPerSec') && isfinite(bestStaticMsCase.truthSummary.fdRateTrueHzPerSec)
    fdRateUse = bestStaticMsCase.truthSummary.fdRateTrueHzPerSec;
    return;
  end
end
end

function textOut = localDescribeBestStaticFdRateSource(bestStaticMsCase)
textOut = "srcStaticFallback";
if isstruct(bestStaticMsCase) && isfield(bestStaticMsCase, 'truth') && isstruct(bestStaticMsCase.truth)
  if isfield(bestStaticMsCase.truth, 'fdRateTrueHzPerSec') && isfinite(bestStaticMsCase.truth.fdRateTrueHzPerSec)
    textOut = "bestStaticMsCase.truth.fdRateTrueHzPerSec";
    return;
  end
  if isfield(bestStaticMsCase.truth, 'fdRateFit') && isfinite(bestStaticMsCase.truth.fdRateFit)
    textOut = "bestStaticMsCase.truth.fdRateFit";
    return;
  end
end
if isstruct(bestStaticMsCase) && isfield(bestStaticMsCase, 'truthSummary') && isstruct(bestStaticMsCase.truthSummary)
  if isfield(bestStaticMsCase.truthSummary, 'fdRateTrueHzPerSec') && isfinite(bestStaticMsCase.truthSummary.fdRateTrueHzPerSec)
    textOut = "bestStaticMsCase.truthSummary.fdRateTrueHzPerSec";
    return;
  end
end
end

function gapNorm = localStateGapNorm(stateA, stateB)
vecA = [stateA.latlon(:); stateA.fdRef; stateA.fdRate];
vecB = [stateB.latlon(:); stateB.fdRef; stateB.fdRate];
gapNorm = norm(vecA - vecB);
end

function stageText = localInferCollapseStage(srcCpK, srcStatic, packedCpK, packedStatic, outerTrace)
stateTol = 1e-9;
sourceCollapsed = calcLatlonAngleError(srcCpK.latlon, srcStatic.latlon) <= stateTol && ...
  abs(srcCpK.fdRef - srcStatic.fdRef) <= stateTol && ...
  abs(srcCpK.fdRate - srcStatic.fdRate) <= stateTol;
packedCollapsed = calcLatlonAngleError(packedCpK.latlon, packedStatic.latlon) <= stateTol && ...
  abs(packedCpK.fdRef - packedStatic.fdRef) <= stateTol && ...
  abs(packedCpK.fdRate - packedStatic.fdRate) <= stateTol;
outerCollapsed = outerTrace.hasCandidates && outerTrace.count >= 2 && ...
  outerTrace.bestPairDoaGapDeg <= stateTol && ...
  outerTrace.bestPairFdRefGapHz <= stateTol && ...
  outerTrace.bestPairFdRateGapHzPerSec <= stateTol;
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

function outerTrace = localDescribeOuterCandidateSet(initCandUnknown, initLen)
outerTrace = struct();
outerTrace.hasCandidates = false;
outerTrace.count = 0;
outerTrace.tagList = strings(1, 0);
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
if size(initParamMat, 2) >= 2
  bestNorm = inf;
  for iCand = 1:(size(initParamMat, 2) - 1)
    for jCand = (iCand + 1):size(initParamMat, 2)
      vecI = initParamMat(:, iCand);
      vecJ = initParamMat(:, jCand);
      normUse = norm(vecI - vecJ);
      if normUse < bestNorm
        stateI = localUnpackUnknownInit(vecI);
        stateJ = localUnpackUnknownInit(vecJ);
        outerTrace.bestPairGapNorm = normUse;
        outerTrace.bestPairDoaGapDeg = calcLatlonAngleError(stateI.latlon, stateJ.latlon);
        outerTrace.bestPairFdRefGapHz = abs(stateI.fdRef - stateJ.fdRef);
        outerTrace.bestPairFdRateGapHzPerSec = abs(stateI.fdRate - stateJ.fdRate);
        bestNorm = normUse;
      end
    end
  end
end
end

function [tagList, initParamMat] = localExtractCandidateSet(candObj, initLen)
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
  srcTag = strings(1, 0);
  if isfield(candObj, 'startTagList')
    srcTag = reshape(string(candObj.startTagList), 1, []);
  elseif isfield(candObj, 'startTag')
    srcTag = reshape(string(candObj.startTag), 1, []);
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
count = size(initParamMat, 2);
tagList = strings(1, count);
for iEntry = 1:count
  tagList(iEntry) = sprintf('%s%d', tagPrefix, iEntry);
end
end

function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%.6f', x), valueVec, 'UniformOutput', false), ', '), ']'];
end

function textOut = localFormatStringRow(valueVec)
valueVec = reshape(string(valueVec), 1, []);
textOut = ['[', strjoin(cellstr(valueVec), ', '), ']'];
end

function fixture = localBuildRegressionFixture()
rng(253);
numUsr = 1; numFrame = 7; refFrameIdx = ceil(numFrame / 2); frameIntvlSec = 1 / 750;
sampleRate = 122.88e6; symbolRate = 3.84e6; osf = sampleRate / symbolRate; numSym = 32; carrierFreq = 11.7e9; lightSpeed = 299792458; wavelen = lightSpeed / carrierFreq; E = referenceEllipsoid('sphere'); gridSize = [41, 41]; searchRange = [32.78, 42.78; 31.59, 41.59]; weightSweepAlpha = [0; 0.25; 0.5; 1];
utcRef = datetime([2026, 03, 18, 17, 08, 0], 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec); usrLla = [37.78; 36.59; 0];
tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle")); arrUpa = createUpa([4, 4], wavelen / 2);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, [1, 2], [], arrUpa, 15, 55, "satellite", 1, refFrameIdx); sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRef, sceneRef.satPosEci, sceneRef.satVelEci); otherSatIdxLocal = 3 - refSatIdxLocal; otherSatIdxGlobal = sceneRef.satIdx(otherSatIdxLocal);
linkParamCell = cell(1, sceneSeq.numFrame); for iFrame = 1:sceneSeq.numFrame, linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen); end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1); truth.latlonTrueDeg = usrLla(1:2, 1); truth.fdRefTrueHz = truth.fdRefFit; truth.fdRateTrueHzPerSec = truth.fdRateFit; truth.refSatIdxGlobal = sceneRef.satIdx(refSatIdxLocal);
[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', 1); pulseOpt = struct(); pulseOpt.rolloff = 0.25; pulseOpt.span = 8; [pilotWave, waveInfo] = genPilotWaveform(pilotSym, symbolRate, osf, 'rrc', pulseOpt);
snapOpt = struct(); snapOpt.spatial.model = 'dynamic'; snapOpt.spatial.refFrameIdx = sceneSeq.refFrameIdx; snapOpt.phase.timeModel = 'global'; snapOpt.phase.frameModel = 'shared'; snapOpt.phase.sharedPhase = 2 * pi * rand(sceneSeq.numSat, sceneSeq.numUser); snapOpt.wave.delayModel = 'phaseOnly'; snapOpt.wave.timeRef = 'zero'; snapOpt.wave.carrierPhaseModel = 'none'; snapOpt.precomp.linkParamCell = linkParamCell;
pathGainCell = repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame);
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, 0, pathGainCell, snapOpt); rxSigRef = rxSigCell{sceneSeq.refFrameIdx};
sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, refSatIdxLocal); sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx}; rxSigMfRefOnly = selectRxSigBySat(rxSigCell, refSatIdxLocal, 'multiFrame'); viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, gridSize, searchRange, E, struct('sceneSeq', sceneSeqRefOnly));
sceneOtherOnly = selectSatScene(sceneRef, otherSatIdxLocal); rxSigOtherOnly = selectRxSigBySat(rxSigRef, otherSatIdxLocal, 'singleFrame'); viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, gridSize, searchRange, E);
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, gridSize, searchRange, E, struct('sceneSeq', sceneSeq));
fdRange = localExpandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4); fdRateRange = localExpandRangeToTruth([-1e4, 0], [truth.fdRateFit; truth.fdRateSatTrueHzPerSec(:)], 0.1, 5e2);
doaOnlyOpt = struct(); doaOnlyOpt.useLogObjective = true; staticBaseOpt = struct(); staticBaseOpt.useLogObjective = true; dynBaseOpt = struct(); dynBaseOpt.useLogObjective = true; dynBaseOpt.initFdCount = 81; dynBaseOpt.useAccessMask = false; dynBaseOpt.phaseMode = 'continuous'; dynBaseOpt.steeringMode = 'framewise'; dynBaseOpt.steeringRefFrameIdx = sceneSeq.refFrameIdx; dynBaseOpt.debugEnable = true; dynBaseOpt.debugStoreEvalTrace = false; dynBaseOpt.debugMaxEvalTrace = 120;
caseBundle = buildDoaDopplerStaticTransitionBundle(viewRefOnly, viewOtherOnly, viewMs, wavelen, pilotWave, carrierFreq, waveInfo.sampleRate, fdRange, truth, otherSatIdxGlobal, false, doaOnlyOpt, staticBaseOpt, weightSweepAlpha, [0.01; 0.01]);
bestStaticMsCase = caseBundle.bestStaticMsCase; staticMsOpt = caseBundle.staticMsOpt; dynMsKnownOpt = dynBaseOpt; dynMsKnownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:); dynMsKnownOpt.initDoaHalfWidth = [0.01; 0.01]; dynMsKnownOpt.fdRateMode = 'known'; dynMsKnownOpt.fdRateKnown = truth.fdRateFit; dynMsKnownOpt.optimOpt = struct('MaxIterations', 60, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12); dynMsUnknownOpt = dynBaseOpt; dynMsUnknownOpt.initDoaParam = bestStaticMsCase.estResult.doaParamEst(:); dynMsUnknownOpt.initDoaHalfWidth = [0.002; 0.002]; dynMsUnknownOpt.fdRateMode = 'unknown'; dynMsUnknownOpt.optimOpt = struct('MaxIterations', 80, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12);
fixture = struct(); fixture.sceneSeq = sceneSeq; fixture.truth = truth; fixture.viewMs = viewMs; fixture.pilotWave = pilotWave; fixture.sampleRate = waveInfo.sampleRate; fixture.carrierFreq = carrierFreq; fixture.fdRange = fdRange; fixture.fdRateRange = fdRateRange; fixture.bestStaticMsCase = bestStaticMsCase; fixture.staticMsOpt = staticMsOpt; fixture.dynMsKnownOpt = dynMsKnownOpt; fixture.dynMsUnknownOpt = dynMsUnknownOpt;
end

function value = localExpandRangeToTruth(defaultRange, truthVals, fracPad, absPad)
truthVals = truthVals(isfinite(truthVals)); if isempty(truthVals), value = defaultRange; return; end
minTruth = min(truthVals); maxTruth = max(truthVals); spanTruth = max(maxTruth - minTruth, eps); pad = max(absPad, fracPad * spanTruth); value = [min(defaultRange(1), minTruth - pad), max(defaultRange(2), maxTruth + pad)];
end

function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath')); projectRoot = fileparts(fileparts(scriptDir)); addpath(genpath(projectRoot));
end
