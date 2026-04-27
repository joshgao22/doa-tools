function regressionMfCpSupportCollapsePenalty(varargin)
% Regression check for MF continuous-phase support-collapse penalties.
% The wrong-tooth candidate should expose worse support-collapse
% diagnostics and a smaller continuous-phase consistency score in the
% shared-phase summary.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();

  fprintf('Running regressionMfCpSupportCollapsePenalty ...\n');

modelOpt = struct();
modelOpt.useLogObjective = false;
modelOpt.phaseMode = 'continuous';
modelOpt.steeringMode = 'framewise';
modelOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
modelOpt.fdRateMode = 'known';
modelOpt.fdRateKnown = fixture.truth.fdRateFit;
modelOpt.continuousPhaseConsistencyWeight = 0.05;
modelOpt.continuousPhaseCollapsePenaltyWeight = 0.10;
modelOpt.continuousPhaseNegativeProjectionPenaltyWeight = 0.10;

[model, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, fixture.fdRateRange, modelOpt);

toothStepHz = 1 / fixture.frameIntvlSec;
hypTruth = struct();
hypTruth.localDoaArr = fixture.sceneSeq.localDoa;
hypTruth.fdSat = fixture.truth.fdSatFitHz(:);
hypTruth.fdRateSat = fixture.truth.fdRateSatTrueHzPerSec(:);

hypWrong = hypTruth;
hypWrong.fdSat = hypWrong.fdSat + toothStepHz;

[objTruth, profTruth, auxTruth] = evalDoaDopplerDynProfileLike(model, hypTruth);
[objWrong, profWrong, auxWrong] = evalDoaDopplerDynProfileLike(model, hypWrong);

truthSupport = localGetPenaltyMetric(profTruth, auxTruth, 'effectiveFrameSupportRatioSat', @min, NaN);
wrongSupport = localGetPenaltyMetric(profWrong, auxWrong, 'effectiveFrameSupportRatioSat', @min, NaN);
truthNegRatio = localGetPenaltyMetric(profTruth, auxTruth, 'negativeProjectionRatioSat', @max, NaN);
wrongNegRatio = localGetPenaltyMetric(profWrong, auxWrong, 'negativeProjectionRatioSat', @max, NaN);
truthConsistency = localGetPenaltyMetric(profTruth, auxTruth, 'consistencyScoreSat', @min, NaN);
wrongConsistency = localGetPenaltyMetric(profWrong, auxWrong, 'consistencyScoreSat', @min, NaN);

if ~(localHasPenaltyMetric(profTruth, auxTruth, 'effectiveFrameSupportRatioSat') && ...
    localHasPenaltyMetric(profWrong, auxWrong, 'negativeProjectionRatioSat') && ...
    localHasPenaltyMetric(profTruth, auxTruth, 'consistencyScoreSat') && ...
    localHasPenaltyMetric(profWrong, auxWrong, 'consistencyScoreSat'))
  error('regressionMfCpSupportCollapsePenalty:MissingConsistencyDiagnostics', ...
    'The evaluator must expose the CP support-collapse diagnostics and scores.');
end

metricTol = 1e-9;
supportCollapsed = isfinite(truthSupport) && isfinite(wrongSupport) && (wrongSupport < truthSupport - metricTol);
negativeIncreased = isfinite(truthNegRatio) && isfinite(wrongNegRatio) && (wrongNegRatio > truthNegRatio + metricTol);
consistencyDegraded = isfinite(truthConsistency) && isfinite(wrongConsistency) && (wrongConsistency < truthConsistency - metricTol);
objectivePenalized = isfinite(objTruth) && isfinite(objWrong) && (objWrong > objTruth + 1e-6);

if ~(supportCollapsed || negativeIncreased || consistencyDegraded || objectivePenalized)
  error('regressionMfCpSupportCollapsePenalty:WrongToothPenaltyMissing', ...
    ['The wrong-tooth candidate should trigger one CP penalty signal ', ...
     '(support collapse, negative projection, consistency score, or objective loss).']);
end

  fprintf('  truth obj                 : %.6f\n', objTruth);
  fprintf('  wrong-tooth obj           : %.6f\n', objWrong);
  fprintf('  truth min support ratio   : %.6f\n', truthSupport);
  fprintf('  wrong min support ratio   : %.6f\n', wrongSupport);
  fprintf('  truth max neg ratio       : %.6f\n', truthNegRatio);
  fprintf('  wrong max neg ratio       : %.6f\n', wrongNegRatio);
  fprintf('  truth min consistency     : %.6f\n', truthConsistency);
  fprintf('  wrong min consistency     : %.6f\n', wrongConsistency);
  fprintf('PASS: regressionMfCpSupportCollapsePenalty\n');




end

function tf = localHasPenaltyMetric(profStruct, auxStruct, fieldName)
%LOCALHASPENALTYMETRIC Check whether one CP penalty metric is exposed.

tf = false;
if isstruct(profStruct) && isfield(profStruct, fieldName) && ~isempty(profStruct.(fieldName))
  tf = true;
  return;
end
if isstruct(auxStruct) && isfield(auxStruct, fieldName) && ~isempty(auxStruct.(fieldName))
  tf = true;
end
end


function value = localGetPenaltyMetric(profStruct, auxStruct, fieldName, reducer, defaultValue)
%LOCALGETPENALTYMETRIC Read one CP penalty metric from prof/aux structs.

value = defaultValue;
rawValue = [];
if isstruct(profStruct) && isfield(profStruct, fieldName)
  rawValue = profStruct.(fieldName);
elseif isstruct(auxStruct) && isfield(auxStruct, fieldName)
  rawValue = auxStruct.(fieldName);
end
if isempty(rawValue)
  return;
end
rawValue = rawValue(isfinite(rawValue));
if isempty(rawValue)
  return;
end
value = reducer(rawValue);
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
fixture.frameIntvlSec = frameIntvlSec;
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
