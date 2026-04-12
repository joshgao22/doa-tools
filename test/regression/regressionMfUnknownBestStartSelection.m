% Regression check for CP-U multi-start best-start selection.
% This script focuses on one narrow contract only:
%   1) the outer CP-U multi-start wrapper must evaluate both fromCpK and
%      fromStatic warm starts;
%   2) the selected candidate must match the documented selection rule that
%      prefers near-best resolved fits with real release motion;
%   3) the selected tag reported by runDynamicDoaDopplerCase must agree
%      with summarizeDynamicMultiStart.
clear(); close all; clc;

localAddProjectPath();
fixture = localBuildRegressionFixture();

fprintf('Running regressionMfUnknownBestStartSelection ...\n');

initParamStaticMs = buildDynamicInitParamFromCase(fixture.bestStaticMsCase, true, fixture.truth.fdRateFit);
if isempty(initParamStaticMs)
  error('regressionMfUnknownBestStartSelection:MissingStaticInit', ...
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
caseDynMsUnknown = runDynamicDoaDopplerCase('MS-MF-CP-U', 'multi', ...
  fixture.viewMs, fixture.truth, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.fdRange, fixture.fdRateRange, false, ...
  fixture.dynMsUnknownOpt, false, [], initCandUnknown);

multiStartDiag = caseDynMsUnknown.estResult.aux.multiStart;
summaryTable = summarizeDynamicMultiStart(caseDynMsUnknown, fixture.truthSummary);
if isempty(summaryTable) || height(summaryTable) < 2
  error('regressionMfUnknownBestStartSelection:MissingMultiStartRows', ...
    'CP-U multi-start summary must contain both outer warm-start rows.');
end

actualTagList = reshape(string(multiStartDiag.tagList), [], 1);
requiredTag = ["fromCpK"; "fromStatic"];
if ~all(ismember(requiredTag, actualTagList))
  error('regressionMfUnknownBestStartSelection:MissingRequiredTags', ...
    'CP-U multi-start must evaluate both fromCpK and fromStatic candidates.');
end

expectedIdx = localSelectBestDynamicCandidateContract( ...
  multiStartDiag.fvalList, multiStartDiag.isResolvedList, multiStartDiag.tagList, ...
  multiStartDiag.iterationList, multiStartDiag.objImproveList, multiStartDiag.moveNormList);
expectedTag = string(multiStartDiag.tagList(expectedIdx));
actualIdx = multiStartDiag.selectedIdx;
actualTag = string(multiStartDiag.selectedTag);
summarySelectedMask = logical(summaryTable.isSelected);
if nnz(summarySelectedMask) ~= 1
  error('regressionMfUnknownBestStartSelection:InvalidSummarySelection', ...
    'summarizeDynamicMultiStart must mark exactly one selected outer candidate.');
end
summarySelectedTag = string(summaryTable.tag(summarySelectedMask));

if actualIdx ~= expectedIdx || actualTag ~= expectedTag
  error('regressionMfUnknownBestStartSelection:SelectedCandidateMismatch', ...
    ['runDynamicDoaDopplerCase selected %s (idx %d), but the documented ', ...
     'selection rule prefers %s (idx %d).'], actualTag, actualIdx, expectedTag, expectedIdx);
end
if summarySelectedTag ~= actualTag
  error('regressionMfUnknownBestStartSelection:SummaryTagMismatch', ...
    'summarizeDynamicMultiStart selected %s but multiStartDiag selected %s.', ...
    summarySelectedTag, actualTag);
end

fprintf('  candidate tags           : %s\n', localFormatStringRow(actualTagList));
fprintf('  candidate fval           : %s\n', localFormatNumericRow(multiStartDiag.fvalList));
fprintf('  candidate obj improve    : %s\n', localFormatNumericRow(multiStartDiag.objImproveList));
fprintf('  candidate move norm      : %s\n', localFormatNumericRow(multiStartDiag.moveNormList));
fprintf('  selected tag             : %s\n', actualTag);
fprintf('  summary selected tag     : %s\n', summarySelectedTag);
fprintf('PASS: regressionMfUnknownBestStartSelection\n');


function selectedIdx = localSelectBestDynamicCandidateContract(fvalList, isResolvedList, tagList, iterationList, objImproveList, moveNormList)
%LOCALSELECTBESTDYNAMICCANDIDATECONTRACT Rebuild the documented selection rule.

selectedIdx = 1;
fvalList = reshape(fvalList, [], 1);
isResolvedList = reshape(logical(isResolvedList), [], 1);
tagList = reshape(string(tagList), [], 1);
iterationList = reshape(iterationList, [], 1);
objImproveList = reshape(objImproveList, [], 1);
moveNormList = reshape(moveNormList, [], 1);

validResolved = find(isResolvedList & isfinite(fvalList));
if isempty(validResolved)
  validAny = find(isfinite(fvalList));
  if ~isempty(validAny)
    [~, idxRel] = min(fvalList(validAny));
    selectedIdx = validAny(idxRel);
  end
  return;
end

bestFval = min(fvalList(validResolved));
objTol = max([1e-3, 1e-7 * max(abs(bestFval), 1), 250]);
nearBest = validResolved(abs(fvalList(validResolved) - bestFval) <= objTol);
if numel(nearBest) == 1
  selectedIdx = nearBest(1);
  return;
end

scoreImprove = objImproveList(nearBest);
scoreImprove(~isfinite(scoreImprove)) = -inf;
scoreMove = moveNormList(nearBest);
scoreMove(~isfinite(scoreMove)) = -inf;
scoreIter = iterationList(nearBest);
scoreIter(~isfinite(scoreIter)) = -inf;
preferStatic = contains(lower(tagList(nearBest)), 'static');

scoreMat = [scoreImprove(:), scoreMove(:), scoreIter(:), double(preferStatic(:))];
[~, idxRel] = sortrows(scoreMat, [-1, -2, -3, -4]);
selectedIdx = nearBest(idxRel(1));
end


function textOut = localFormatNumericRow(valueVec)
%LOCALFORMATNUMERICROW Format one numeric row for compact logging.

valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end


function textOut = localFormatStringRow(valueVec)
%LOCALFORMATSTRINGROW Format one string row for compact logging.

valueVec = reshape(string(valueVec), 1, []);
cellText = cellstr(valueVec);
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


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
