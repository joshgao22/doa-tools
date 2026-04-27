function [seedCandidateList, periodicCaseCell, periodicSummaryCell] = runSimplePeriodicRefineCandidates( ...
  displayName, satMode, periodicFixture, staticSeedCase, selectedSubsetCase, ...
  toothStepHz, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, truthUse)
%RUNSIMPLEPERIODICREFINECANDIDATES Run narrow frozen same-tooth replay candidates.
% This helper keeps the simplified flow entry focused on orchestration while
% centralizing the periodic replay candidate loop.

arguments
  displayName (1,1) string
  satMode (1,1) string {mustBeMember(satMode,["single","multi"])}
  periodicFixture (1,1) struct
  staticSeedCase (1,1) struct
  selectedSubsetCase (1,1) struct
  toothStepHz (1,1) double
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  optVerbose (1,1) logical
  flowOpt (1,1) struct
  debugTruth (1,1) struct
  truthUse (1,1) struct
end

[viewUse, ~] = localResolveViewAndTruth(satMode, periodicFixture);
seedCandidateList = resolveSimplePeriodicRefineDoaSeedCandidates(flowOpt, satMode, staticSeedCase, selectedSubsetCase);
numFrozenCandidate = numel(seedCandidateList);
periodicCaseCell = cell(numFrozenCandidate, 1);
periodicSummaryCell = cell(numFrozenCandidate, 1);
for iCand = 1:numFrozenCandidate
  [periodicCaseCell{iCand}, periodicSummaryCell{iCand}] = localRunPeriodicReplayCandidate( ...
    displayName, viewUse, truthUse, periodicFixture, selectedSubsetCase, seedCandidateList(iCand), ...
    pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz);
end
end

function [caseUse, summaryUse] = localRunPeriodicReplayCandidate(displayName, viewUse, truthUse, periodicFixture, selectedSubsetCase, ...
  seedCandidate, pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, debugTruth, toothStepHz)
fdRangeUse = [selectedSubsetCase.estResult.fdRefEst - flowOpt.periodicRefineFdHalfWidthHz, ...
  selectedSubsetCase.estResult.fdRefEst + flowOpt.periodicRefineFdHalfWidthHz];
fdRateRangeUse = [selectedSubsetCase.estResult.fdRateEst - flowOpt.periodicRefineFdRateHalfWidthHzPerSec, ...
  selectedSubsetCase.estResult.fdRateEst + flowOpt.periodicRefineFdRateHalfWidthHzPerSec];

dynOpt = flowOpt.dynBaseOpt;
dynOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynOpt.initDoaParam = seedCandidate.doaInitParam(:);
dynOpt.initDoaHalfWidth = seedCandidate.doaHalfWidthDeg(:);
dynOpt.enableFdAliasUnwrap = true;
dynOpt.freezeDoa = logical(seedCandidate.freezeDoa);
dynOpt.disableUnknownDoaReleaseFloor = true;
dynOpt.unknownDoaReleaseHalfWidth = seedCandidate.doaHalfWidthDeg(:);

initParamSeed = buildDynamicInitParamFromCase(selectedSubsetCase, false, selectedSubsetCase.estResult.fdRateEst);
if ~isempty(initParamSeed) && numel(initParamSeed) >= numel(seedCandidate.doaInitParam)
  initParamSeed(1:numel(seedCandidate.doaInitParam)) = seedCandidate.doaInitParam(:);
end
initCandidate = struct( ...
  'startTag', string(seedCandidate.startTag), ...
  'initParam', initParamSeed, ...
  'initDoaParam', seedCandidate.doaInitParam(:), ...
  'initDoaHalfWidth', seedCandidate.doaHalfWidthDeg(:), ...
  'freezeDoa', logical(seedCandidate.freezeDoa));

caseUse = runDynamicDoaDopplerCase(displayName, localInferSatMode(viewUse), viewUse, truthUse, pilotWave, ...
  carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, ...
  optVerbose, dynOpt, false, debugTruth, initCandidate);
summaryUse = buildDynamicUnknownCaseSummary(caseUse, toothStepHz, truthUse);
summaryUse.seedSource = string(seedCandidate.seedSource);
summaryUse.startTag = string(seedCandidate.startTag);
summaryUse.subsetDriftFromStaticDeg = seedCandidate.subsetDriftFromStaticDeg;
end

function satMode = localInferSatMode(viewUse)
satMode = "multi";
if isstruct(viewUse) && isfield(viewUse, 'numSat')
  if isequal(viewUse.numSat, 1)
    satMode = "single";
  end
end
end

function [viewUse, truthUse] = localResolveViewAndTruth(satMode, fixtureUse)
if satMode == "single"
  viewUse = fixtureUse.viewRefOnly;
  truthUse = fixtureUse.truth;
else
  viewUse = fixtureUse.viewMs;
  truthUse = fixtureUse.truth;
end
end
