function regressionMfDoaBasinEntrySplit(varargin)
% Regression check for MF DoA basin-entry branch wiring.
% This regression locks the branch-level contract around the SS-MF DoA
% basin-entry split: the default single-sat continuous MF path must evaluate
% the entry candidates, while disableDoaBasinEntry must bypass them.  It does
% not assert RMSE/CRB performance or candidate superiority.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fprintf('Running regressionMfDoaBasinEntrySplit ...\n');

fixture = localBuildRegressionFixture();
baseModelOpt = localBuildBaseModelOpt(fixture, verbose);

[solveEnabled, optimInfoEnabled] = localRunBranchCase(fixture, baseModelOpt);
entryDiag = localRequireDoaBasinEntryDiag(optimInfoEnabled, 'enabled');
localCheckEnabledDiag(entryDiag, baseModelOpt);
localCheckDoaInsideParent(solveEnabled.optVar(1:2), baseModelOpt, 'enabled-final');
localCheckReferenceDopplerInvariant(solveEnabled.finalEvalDiag, 'enabled');

modelOptDisabled = baseModelOpt;
modelOptDisabled.disableDoaBasinEntry = true;
[solveDisabled, optimInfoDisabled] = localRunBranchCase(fixture, modelOptDisabled);
disableDiag = localRequireDoaBasinEntryDiag(optimInfoDisabled, 'disabled');
localCheckDisabledDiag(disableDiag);
localCheckDoaInsideParent(solveDisabled.optVar(1:2), baseModelOpt, 'disabled-final');
localCheckReferenceDopplerInvariant(solveDisabled.finalEvalDiag, 'disabled');

fprintf('  enabled reason      : %s\n', string(entryDiag.reason));
fprintf('  enabled candidates  : %d\n', entryDiag.numCandidate);
fprintf('  disabled reason     : %s\n', string(disableDiag.reason));
fprintf('  selected variant    : %s\n', string(entryDiag.selectedVariant));
fprintf('PASS: regressionMfDoaBasinEntrySplit\n');
end

function [solveUse, optimInfo] = localRunBranchCase(fixture, modelOpt)
%LOCALRUNBRANCHCASE Build and solve one branch case.

[model, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, modelOpt);
[initParam, ~] = buildDoaDopplerMfInit(model, []);
[solveUse, optimInfo, ~] = solveDoaDopplerMfBranches(model, initParam, modelOpt.verbose);
end

function modelOpt = localBuildBaseModelOpt(fixture, verbose)
%LOCALBUILDBASEMODELOPT Build one lightweight SS-MF branch setup.

modelOpt = struct();
modelOpt.useLogObjective = false;
modelOpt.phaseMode = 'continuous';
modelOpt.fdRateMode = 'unknown';
modelOpt.steeringMode = 'framewise';
modelOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
modelOpt.initDoaParam = fixture.truth.latlonTrueDeg(:) + [0.018; -0.014];
modelOpt.initDoaHalfWidth = [0.003; 0.003];
modelOpt.doaBasinEntryHalfWidthList = [0.012, 0.024, 0.048];
modelOpt.doaBasinEntryAdoptionMode = 'diagnostic-only';
modelOpt.disableUnknownWarmAnchor = true;
modelOpt.initFdCount = 17;
modelOpt.phaseGridCount = 16;
modelOpt.verbose = verbose;
modelOpt.optimOpt = struct( ...
  'MaxIterations', 25, ...
  'StepTolerance', 1e-9, ...
  'OptimalityTolerance', 1e-8, ...
  'ConstraintTolerance', 1e-9);
end

function entryDiag = localRequireDoaBasinEntryDiag(optimInfo, caseTag)
%LOCALREQUIREDOABASINENTRYDIAG Fetch the basin-entry diagnostic.

if ~isstruct(optimInfo) || ~isfield(optimInfo, 'doaBasinEntry')
  error('regressionMfDoaBasinEntrySplit:MissingDoaBasinEntryDiag', ...
    '[%s] optimInfo must expose doaBasinEntry after branch solving.', caseTag);
end
entryDiag = optimInfo.doaBasinEntry;
end

function localCheckEnabledDiag(entryDiag, modelOpt)
%LOCALCHECKENABLEDDIAG Verify the evaluated SS-MF entry diagnostic.

if ~isstruct(entryDiag) || ~isfield(entryDiag, 'enabled') || ~logical(entryDiag.enabled)
  error('regressionMfDoaBasinEntrySplit:EntryNotEnabled', ...
    'The default SS-MF continuous branch must evaluate DoA basin-entry.');
end
if ~isfield(entryDiag, 'reason') || string(entryDiag.reason) ~= "evaluated"
  error('regressionMfDoaBasinEntrySplit:UnexpectedEnabledReason', ...
    'Enabled DoA basin-entry must report reason="evaluated".');
end
if ~isfield(entryDiag, 'adoptionMode') || string(entryDiag.adoptionMode) ~= "diagnostic-only"
  error('regressionMfDoaBasinEntrySplit:UnexpectedAdoptionMode', ...
    'This regression expects diagnostic-only basin-entry adoption.');
end
if ~isfield(entryDiag, 'numCandidate') || entryDiag.numCandidate ~= 1
  error('regressionMfDoaBasinEntrySplit:UnexpectedCandidateCount', ...
    'The collapsed actual DoA entry paths must evaluate exactly one candidate.');
end
if ~isfield(entryDiag, 'numCandidateProposed') || entryDiag.numCandidateProposed ~= 3
  error('regressionMfDoaBasinEntrySplit:UnexpectedProposedCandidateCount', ...
    'The regression must propose three distinct wide DoA boxes.');
end
if ~isfield(entryDiag, 'numCandidateAfterDedup') || entryDiag.numCandidateAfterDedup ~= 1
  error('regressionMfDoaBasinEntrySplit:UnexpectedDedupCandidateCount', ...
    'Parent-envelope clipping should collapse the three wide boxes to one actual path.');
end
if ~isfield(entryDiag, 'numCandidateDeduped') || entryDiag.numCandidateDeduped ~= 2
  error('regressionMfDoaBasinEntrySplit:UnexpectedDedupCount', ...
    'Exactly two duplicated actual DoA entry paths should be skipped.');
end
if ~isfield(entryDiag, 'numCandidateSkippedNoOverlap') || entryDiag.numCandidateSkippedNoOverlap ~= 0
  error('regressionMfDoaBasinEntrySplit:UnexpectedNoOverlapCount', ...
    'The default regression boxes should overlap the parent envelope.');
end
if ~isfield(entryDiag, 'entryTable') || numel(entryDiag.entryTable) ~= entryDiag.numCandidate
  error('regressionMfDoaBasinEntrySplit:EntryTableSizeMismatch', ...
    'The basin-entry table size must match numCandidate.');
end
row = entryDiag.entryTable(1);
requiredField = ["tag", "halfWidth1", "halfWidth2", "centerSource", ...
  "adoptionMode", "baselineFval", "entryFval", "selectedFval", ...
  "selectedVariant", "adoptedOverPreviousBest"];
for iField = 1:numel(requiredField)
  fieldName = char(requiredField(iField));
  if ~isfield(row, fieldName)
    error('regressionMfDoaBasinEntrySplit:MissingEntryRowField', ...
      'Basin-entry row is missing field: %s.', fieldName);
  end
end
if string(row.adoptionMode) ~= "diagnostic-only"
  error('regressionMfDoaBasinEntrySplit:EntryRowAdoptionModeMismatch', ...
    'Entry row adoptionMode must follow the model adoption mode.');
end
if logical(row.adoptedOverPreviousBest)
  error('regressionMfDoaBasinEntrySplit:DiagnosticOnlyAdopted', ...
    'Diagnostic-only basin-entry must not adopt an entry candidate.');
end
if ~isfinite(row.baselineFval) || ~isfinite(row.entryFval) || ~isfinite(row.selectedFval)
  error('regressionMfDoaBasinEntrySplit:InvalidEntryObjective', ...
    'Basin-entry baseline, entry, and selected objectives must be finite.');
end
localCheckDoaInsideParent([row.entryLatDeg; row.entryLonDeg], modelOpt, 'entry-candidate');
localCheckDoaInsideParent([row.selectedLatDeg; row.selectedLonDeg], modelOpt, 'selected-candidate');
if ~isfield(entryDiag, 'selectedVariant') || strlength(string(entryDiag.selectedVariant)) == 0
  error('regressionMfDoaBasinEntrySplit:MissingSelectedVariant', ...
    'Basin-entry diagnostic must expose selectedVariant.');
end
end


function localCheckDoaInsideParent(doaParam, modelOpt, caseTag)
%LOCALCHECKDOAINSIDEPARENT Guard child DoA boxes against parent-envelope escape.

center = reshape(modelOpt.initDoaParam, [], 1);
halfWidth = reshape(modelOpt.initDoaHalfWidth, [], 1);
doaParam = reshape(doaParam, [], 1);
if numel(center) ~= 2 || numel(halfWidth) ~= 2 || numel(doaParam) < 2 || ...
    any(~isfinite(doaParam(1:2)))
  error('regressionMfDoaBasinEntrySplit:InvalidDoaParentCheckInput', ...
    '[%s] DoA parent-envelope check received invalid input.', caseTag);
end
lb = center - halfWidth;
ub = center + halfWidth;
tol = 1e-9;
if any(doaParam(1:2) < lb - tol) || any(doaParam(1:2) > ub + tol)
  error('regressionMfDoaBasinEntrySplit:DoaEscapedParentEnvelope', ...
    ['[%s] MF branch returned DoA outside the upper-level initDoaParam ', ...
     '+/- initDoaHalfWidth envelope.'], caseTag);
end
end

function localCheckDisabledDiag(entryDiag)
%LOCALCHECKDISABLEDDIAG Verify the disableDoaBasinEntry bypass.

if ~isstruct(entryDiag) || ~isfield(entryDiag, 'enabled') || logical(entryDiag.enabled)
  error('regressionMfDoaBasinEntrySplit:EntryNotBypassed', ...
    'disableDoaBasinEntry=true must bypass DoA basin-entry.');
end
if ~isfield(entryDiag, 'reason') || string(entryDiag.reason) ~= "disabled-by-option"
  error('regressionMfDoaBasinEntrySplit:UnexpectedDisabledReason', ...
    'Disabled DoA basin-entry must report reason="disabled-by-option".');
end
if ~isfield(entryDiag, 'numCandidate') || entryDiag.numCandidate ~= 0
  error('regressionMfDoaBasinEntrySplit:DisabledCandidateCount', ...
    'Disabled DoA basin-entry must not evaluate candidates.');
end
if isfield(entryDiag, 'entryTable') && ~isempty(entryDiag.entryTable)
  error('regressionMfDoaBasinEntrySplit:DisabledEntryTableNotEmpty', ...
    'Disabled DoA basin-entry must leave an empty entry table.');
end
end

function localCheckReferenceDopplerInvariant(evalDiag, caseTag)
%LOCALCHECKREFERENCEDOPPLERINVARIANT Guard reference Doppler composition.

fdRef = evalDiag.fdRef;
fdSat = reshape(evalDiag.fdSat, [], 1);
deltaFdRef = reshape(evalDiag.deltaFdRef, [], 1);
refSatIdxLocal = evalDiag.refSatIdxLocal;
if isempty(fdSat) || isempty(deltaFdRef) || ~isfinite(refSatIdxLocal)
  error('regressionMfDoaBasinEntrySplit:MissingDopplerDiag', ...
    '[%s] finalEvalDiag must expose fdRef, fdSat, deltaFdRef, and refSatIdxLocal.', caseTag);
end
composeVec = fdRef + deltaFdRef;
tolHz = 1e-9 * max([1; abs(fdSat(:)); abs(composeVec(:)); abs(fdRef)]);
if abs(deltaFdRef(refSatIdxLocal)) > tolHz
  error('regressionMfDoaBasinEntrySplit:RefDeltaFdNonzero', ...
    '[%s] deltaFdRef(refSat) must be zero.', caseTag);
end
if abs(fdSat(refSatIdxLocal) - fdRef) > tolHz
  error('regressionMfDoaBasinEntrySplit:RefFdMismatch', ...
    '[%s] fdSat(refSat) must equal fdRef.', caseTag);
end
if max(abs(fdSat(:) - composeVec(:))) > tolHz
  error('regressionMfDoaBasinEntrySplit:FdCompositionMismatch', ...
    '[%s] fdSat must equal fdRef + deltaFdRef.', caseTag);
end
end

function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact single-sat MF case.

rng(253);

numUsr = 1;
numFrame = 5;
refFrameIdx = ceil(numFrame / 2);
frameIntvlSec = 1 / 750;

sampleRate = 122.88e6;
symbolRate = 3.84e6;
osf = sampleRate / symbolRate;
numSym = 16;
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
arrUpa = createUpa([2, 2], wavelen / 2);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, 1, [], arrUpa, ...
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

pathGainCell = repmat({ones(sceneSeq.numSat, sceneSeq.numUser)}, 1, sceneSeq.numFrame);
noiseVar = 1e-6;
[rxSigCell, ~, ~, ~, ~] = genMultiFrameSnapshots(sceneSeq, pilotWave, carrierFreq, ...
  waveInfo.sampleRate, noiseVar, pathGainCell, snapOpt);

searchRange = [usrLla(1, 1) - 0.08, usrLla(1, 1) + 0.08; ...
               usrLla(2, 1) - 0.08, usrLla(2, 1) + 0.08];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, [15, 15], searchRange, E, ...
  struct('sceneSeq', sceneSeq));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateRange = expandRangeToTruth([-1e4, 0], ...
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
