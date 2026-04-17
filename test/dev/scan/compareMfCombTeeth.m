% Diagnostic compare of the 1/T_f comb teeth in the MF objective.
% This script does not run the optimizer. Instead, it evaluates the dynamic
% profile objective on several fdRef teeth around truth:
%   fdRef = fdRefTruth + k / T_f,  k in {-2,-1,0,1,2}.
% It compares one periodic frame schedule against one non-periodic subset
% schedule while keeping the same current CP model and truth local DoA.
clear(); close all; clc;

localAddProjectPath();

periodicOffsetIdx = -4:5;
masterOffsetIdx = -9:10;
rng(253);
masterNoZero = masterOffsetIdx(masterOffsetIdx ~= 0);
nonPeriodicOffsetIdx = sort([0, masterNoZero(randperm(numel(masterNoZero), 9))]);

fixturePeriodic = localBuildFixture(periodicOffsetIdx);
fixtureNonPeriodic = localBuildFixture(nonPeriodicOffsetIdx);

toothStepHz = 1 / fixturePeriodic.frameIntvlSec;
toothIdxList = -2:2;

fprintf('Running compareMfCombTeeth ...\n');

curvePeriodic = localEvalToothCurve(fixturePeriodic, toothIdxList, toothStepHz);
curveNonPeriodic = localEvalToothCurve(fixtureNonPeriodic, toothIdxList, toothStepHz);

if ~all(isfinite(curvePeriodic.objList)) || ~all(isfinite(curveNonPeriodic.objList))
  error('compareMfCombTeeth:InvalidObjective', ...
    'All evaluated comb-tooth objectives must stay finite.');
end

[bestObjPeriodic, bestIdxPeriodic] = min(curvePeriodic.objList);
[bestObjNonPeriodic, bestIdxNonPeriodic] = min(curveNonPeriodic.objList);
offMask = toothIdxList ~= 0;
minOffGapPeriodic = min(curvePeriodic.gapToTruth(offMask));
minOffGapNonPeriodic = min(curveNonPeriodic.gapToTruth(offMask));

fprintf('  tooth step (Hz)                    : %.6f\n', toothStepHz);
fprintf('  tooth index list                   : %s\n', localFormatIntegerRow(toothIdxList));
fprintf('  periodic offsets                   : %s\n', localFormatIntegerRow(periodicOffsetIdx));
fprintf('  non-periodic offsets               : %s\n', localFormatIntegerRow(nonPeriodicOffsetIdx));
fprintf('  periodic objective                 : %s\n', localFormatNumericRow(curvePeriodic.objList));
fprintf('  periodic gap to truth tooth        : %s\n', localFormatNumericRow(curvePeriodic.gapToTruth));
fprintf('  non-periodic objective             : %s\n', localFormatNumericRow(curveNonPeriodic.objList));
fprintf('  non-periodic gap to truth tooth    : %s\n', localFormatNumericRow(curveNonPeriodic.gapToTruth));
fprintf('  periodic best tooth index          : %d\n', toothIdxList(bestIdxPeriodic));
fprintf('  non-periodic best tooth index      : %d\n', toothIdxList(bestIdxNonPeriodic));
fprintf('  periodic min off-tooth gap         : %.6e\n', minOffGapPeriodic);
fprintf('  non-periodic min off-tooth gap     : %.6e\n', minOffGapNonPeriodic);
fprintf('  off-tooth gap ratio (non/periodic) : %.6f\n', minOffGapNonPeriodic / max(minOffGapPeriodic, eps));

if minOffGapNonPeriodic <= minOffGapPeriodic
  fprintf('  note                               : non-periodic schedule did not increase the off-tooth gap in this realization.\n');
end
if toothIdxList(bestIdxPeriodic) ~= 0
  fprintf('  note                               : periodic schedule prefers a non-truth tooth in this realization.\n');
end
if toothIdxList(bestIdxNonPeriodic) ~= 0
  fprintf('  note                               : non-periodic schedule still prefers a non-truth tooth in this realization.\n');
end
fprintf('PASS: compareMfCombTeeth\n');


function curve = localEvalToothCurve(fixture, toothIdxList, toothStepHz)
modelOpt = struct();
modelOpt.useLogObjective = false;
modelOpt.phaseMode = 'continuous';
modelOpt.fdRateMode = 'unknown';
modelOpt.steeringMode = 'framewise';
modelOpt.steeringRefFrameIdx = fixture.sceneSeq.refFrameIdx;
modelOpt.debugEnable = false;

[modelUse, ~, ~] = buildDoaDopplerMfModel( ...
  fixture.sceneSeq, fixture.rxSigCell, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, fixture.fdRateRange, modelOpt);

truthHyp = struct();
truthHyp.localDoaArr = localExtractSceneLocalDoa(fixture.sceneSeq);
truthHyp.fdSat = fixture.truth.fdSatFitHz(:);
truthHyp.fdRateSat = fixture.truth.fdRateSatTrueHzPerSec(:);

objList = nan(size(toothIdxList));
for iTooth = 1:numel(toothIdxList)
  hypUse = truthHyp;
  hypUse.fdSat = truthHyp.fdSat + toothIdxList(iTooth) * toothStepHz;
  [objList(iTooth), ~] = evalDoaDopplerDynProfileLike(modelUse, hypUse);
end

truthPos = find(toothIdxList == 0, 1);
curve = struct();
curve.objList = reshape(objList, 1, []);
curve.gapToTruth = curve.objList - curve.objList(truthPos);
end


function fixture = localBuildFixture(offsetIdxVec)
rng(253);

numUsr = 1;
refFrameIdx = find(offsetIdxVec == 0, 1);
if isempty(refFrameIdx)
  error('compareMfCombTeeth:MissingReferenceFrame', ...
    'offsetIdxVec must contain one zero-offset reference frame.');
end

frameIntvlSec = 1 / 750;
sampleRate = 512e6;
symbolRate = 128e6;
osf = sampleRate / symbolRate;
numSym = 512;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;
E = referenceEllipsoid('sphere');

snrDb = 10;
pwrSource = 1;
pwrNoise = pwrSource / (10^(snrDb / 10));

usrLla = [37.78; 36.59; 0];
utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(offsetIdxVec(:).' * frameIntvlSec);

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[satIdx, ~] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = satIdx(1);
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, satIdx, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, refFrameIdx);
sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};

linkParamCell = cell(1, sceneSeq.numFrame);
for iFrame = 1:sceneSeq.numFrame
  linkParamCell{iFrame} = getLinkParam(sceneSeq.sceneCell{iFrame}, wavelen);
end
truth = buildDynTruthFromLinkParam(linkParamCell, sceneSeq.timeOffsetSec, sampleRate, 1);
truth.latlonTrueDeg = usrLla(1:2, 1);
truth.selectedSatIdxGlobal = satIdx(:).';

[pilotSym, ~] = genPilotSymbol(numUsr, numSym, 'pn', pwrSource);
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
  sceneSeq, pilotWave, carrierFreq, waveInfo.sampleRate, ...
  pwrNoise, pathGainCell, snapOpt);

searchRange = [usrLla(1, 1) - 5, usrLla(1, 1) + 5; ...
               usrLla(2, 1) - 5, usrLla(2, 1) + 5];
viewMs = buildDoaDopplerEstView(sceneRef, rxSigCell, [50, 50], searchRange, E, struct('sceneSeq', sceneSeq));

fdRange = expandRangeToTruth([-1e5, 0], [truth.fdRefFit; truth.fdSatFitHz(:)], 0.1, 2e4);
fdRateTruthCand = [truth.fdRateFit; truth.fdRateFit + reshape(localGetFieldOrDefault(truth, 'deltaFdRate', []), [], 1)];
fdRateRange = expandRangeToTruth([-1e4, 0], fdRateTruthCand, 0.1, 5e2);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.rxSigCell = rxSigCell;
fixture.viewMs = viewMs;
fixture.truth = truth;
fixture.pilotWave = pilotWave;
fixture.sampleRate = waveInfo.sampleRate;
fixture.carrierFreq = carrierFreq;
fixture.fdRange = fdRange;
fixture.fdRateRange = fdRateRange;
fixture.frameIntvlSec = frameIntvlSec;
end


function localDoa = localExtractSceneLocalDoa(sceneSeq)
localDoaRaw = sceneSeq.localDoa;
if ndims(localDoaRaw) == 3
  localDoa = localDoaRaw;
elseif ndims(localDoaRaw) == 4
  localDoa = localDoaRaw(:, :, 1, :);
else
  error('compareMfCombTeeth:InvalidLocalDoaSize', ...
    'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isempty(dataStruct)
  return;
end
if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    value = dataStruct.(fieldName);
  end
  return;
end
if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    value = dataStruct.(fieldName);
  end
  return;
end
try
  value = dataStruct.(fieldName);
catch
end
end


function textOut = localFormatIntegerRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueVec, 'UniformOutput', false), ', '), ']'];
end


function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false), ', '), ']'];
end


function localAddProjectPath()
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
