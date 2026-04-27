function regressionMfCpVsIpRepresentative(varargin)
%REGRESSIONMFCPVSIPREPRESENTATIVE Guard one representative CP-vs-IP MF case.
% The goal is not to over-encode the full paper claim into one seed. Instead
% this regression keeps a narrower, behavior-stable guardrail on the main
% phase-mode distinction:
%   1) both CP and IP MF estimators must run from the same representative
%      static MS seed without interface drift;
%   2) the truth-centered fdRef tooth scan must stay finite and comparable
%      across phase modes;
%   3) on this representative case, CP should show at least one meaningful
%      advantage (angle / fdRef / off-tooth gap), while not regressing
%      materially in both angle and fdRef at the same time.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

context = buildDynamicDualSatEciContext(struct( ...
  'baseSeed', 253, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', struct('enableSubsetEvalParfor', false, 'minSubsetEvalParfor', 4)));
repeatData = buildDynamicRepeatData(context, 15, 253);
periodicFixture = repeatData.periodicFixture;
truth = periodicFixture.truth;

flowOpt = buildSimpleDynamicFlowOpt(struct('parallelOpt', context.parallelOpt));
flowOpt.verbose = verbose;
staticBundle = buildDoaDopplerStaticTransitionBundle( ...
  periodicFixture.viewRefOnly, periodicFixture.viewOtherOnly, periodicFixture.viewMs, ...
  context.wavelen, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  periodicFixture.fdRange, truth, periodicFixture.otherSatIdxGlobal, ...
  false, flowOpt.doaOnlyOpt, flowOpt.staticBaseOpt, zeros(0, 1), flowOpt.staticMsHalfWidth);
staticSeedCase = staticBundle.bestStaticMsCase;
initParamKnown = buildDynamicInitParamFromCase(staticSeedCase, true, truth.fdRateFit);
if isempty(initParamKnown)
  error('regressionMfCpVsIpRepresentative:MissingStaticSeed', ...
    'The representative CP/IP regression requires one resolved static MS seed.');
end

baseModelOpt = flowOpt.dynBaseOpt;
baseModelOpt.fdRateMode = 'known';
baseModelOpt.fdRateKnown = truth.fdRateFit;
baseModelOpt.steeringMode = 'framewise';
baseModelOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
baseModelOpt.enableFdAliasUnwrap = false;
baseModelOpt.initDoaParam = staticSeedCase.estResult.doaParamEst(:);
baseModelOpt.initDoaHalfWidth = [0.003; 0.003];
baseModelOpt.debugEnable = false;
baseModelOpt.debugStoreEvalTrace = false;
baseModelOpt.debugMaxEvalTrace = 1;

cpOpt = baseModelOpt;
cpOpt.phaseMode = 'continuous';
ipOpt = baseModelOpt;
ipOpt.phaseMode = 'independent';

[cpResult, ~, ~] = estimatorDoaDopplerMlePilotMfOpt( ...
  periodicFixture.sceneSeq, periodicFixture.rxSigCell, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, periodicFixture.viewMs.doaGrid, ...
  periodicFixture.fdRange, [], initParamKnown, false, cpOpt);
[ipResult, ~, ~] = estimatorDoaDopplerMlePilotMfOpt( ...
  periodicFixture.sceneSeq, periodicFixture.rxSigCell, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, periodicFixture.viewMs.doaGrid, ...
  periodicFixture.fdRange, [], initParamKnown, false, ipOpt);

localAssertResolved(cpResult, 'CP');
localAssertResolved(ipResult, 'IP');
localAssertPhaseMode(cpResult, 'continuous');
localAssertPhaseMode(ipResult, 'independent');

cpAngleErrDeg = calcLatlonAngleError(cpResult.doaParamEst(:), truth.latlonTrueDeg(:));
ipAngleErrDeg = calcLatlonAngleError(ipResult.doaParamEst(:), truth.latlonTrueDeg(:));
cpFdErrHz = abs(cpResult.fdRefEst - truth.fdRefFit);
ipFdErrHz = abs(ipResult.fdRefEst - truth.fdRefFit);

[cpModel, ~, ~, ~] = buildDoaDopplerMfModel( ...
  periodicFixture.sceneSeq, periodicFixture.rxSigCell, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, periodicFixture.viewMs.doaGrid, ...
  periodicFixture.fdRange, [], cpOpt);
[ipModel, ~, ~, ~] = buildDoaDopplerMfModel( ...
  periodicFixture.sceneSeq, periodicFixture.rxSigCell, context.pilotWave, ...
  context.carrierFreq, context.waveInfo.sampleRate, periodicFixture.viewMs.doaGrid, ...
  periodicFixture.fdRange, [], ipOpt);

aliasStepHz = localResolveAliasStep(periodicFixture.sceneSeq.timeOffsetSec);
cpLine = localScanAliasToothLine(cpModel, truth, aliasStepHz);
ipLine = localScanAliasToothLine(ipModel, truth, aliasStepHz);

  fprintf('Running regressionMfCpVsIpRepresentative ...\n');
  fprintf('  representative seed           : %d\n', repeatData.taskSeed);
  fprintf('  alias step (Hz)               : %.6f\n', aliasStepHz);
  fprintf('  CP angle err (deg)            : %.6g\n', cpAngleErrDeg);
  fprintf('  IP angle err (deg)            : %.6g\n', ipAngleErrDeg);
  fprintf('  CP |fdRef err| (Hz)           : %.6g\n', cpFdErrHz);
  fprintf('  IP |fdRef err| (Hz)           : %.6g\n', ipFdErrHz);
  fprintf('  CP first-tooth gap            : %.6g\n', cpLine.firstAliasGap);
  fprintf('  IP first-tooth gap            : %.6g\n', ipLine.firstAliasGap);
  fprintf('  CP second-tooth gap           : %.6g\n', cpLine.secondAliasGap);
  fprintf('  IP second-tooth gap           : %.6g\n', ipLine.secondAliasGap);

if ~(isfinite(cpLine.truthObj) && isfinite(ipLine.truthObj) && ...
     isfinite(cpLine.firstAliasGap) && isfinite(ipLine.firstAliasGap) && ...
     isfinite(cpLine.secondAliasGap) && isfinite(ipLine.secondAliasGap))
  error('regressionMfCpVsIpRepresentative:InvalidToothScan', ...
    'The representative CP/IP tooth scan must stay finite for both phase modes.');
end

materialAngleTolDeg = max(5e-4, 0.25 * max(ipAngleErrDeg, eps));
materialFdTolHz = max(1, 0.25 * max(ipFdErrHz, eps));
if (cpAngleErrDeg > ipAngleErrDeg + materialAngleTolDeg) && ...
    (cpFdErrHz > ipFdErrHz + materialFdTolHz)
  error('regressionMfCpVsIpRepresentative:RepresentativeEstimateRegressed', ...
    ['CP is materially worse than IP in both angle and fdRef on the ', ...
     'representative MF case.']);
end

gapAdvTol = max(50, 0.05 * max(abs([cpLine.firstAliasGap; ipLine.firstAliasGap; ...
  cpLine.secondAliasGap; ipLine.secondAliasGap; 1])));
cpBetterAngle = cpAngleErrDeg + materialAngleTolDeg < ipAngleErrDeg;
cpBetterFd = cpFdErrHz + materialFdTolHz < ipFdErrHz;
cpBetterFirstGap = cpLine.firstAliasGap > ipLine.firstAliasGap + gapAdvTol;
cpBetterSecondGap = cpLine.secondAliasGap > ipLine.secondAliasGap + gapAdvTol;

if ~(cpBetterAngle || cpBetterFd || cpBetterFirstGap || cpBetterSecondGap)
  error('regressionMfCpVsIpRepresentative:NoRepresentativeCpAdvantage', ...
    ['On the representative MF case, CP must keep at least one meaningful ', ...
     'advantage over IP (angle / fdRef / off-tooth gap).']);
end

advTagList = strings(0, 1);
if cpBetterAngle
  advTagList(end + 1, 1) = "angle"; %#ok<AGROW>
end
if cpBetterFd
  advTagList(end + 1, 1) = "fdRef"; %#ok<AGROW>
end
if cpBetterFirstGap
  advTagList(end + 1, 1) = "firstGap"; %#ok<AGROW>
end
if cpBetterSecondGap
  advTagList(end + 1, 1) = "secondGap"; %#ok<AGROW>
end
  fprintf('  CP representative advantage   : [%s]\n', strjoin(cellstr(advTagList), ', '));
  fprintf('PASS: regressionMfCpVsIpRepresentative\n');
end

function localAssertResolved(estResult, modeTag)
if ~isstruct(estResult) || ~logical(getDoaDopplerFieldOrDefault(estResult, 'isResolved', false))
  error('regressionMfCpVsIpRepresentative:UnresolvedEstimate', ...
    '%s representative estimate must resolve.', char(string(modeTag)));
end
end

function localAssertPhaseMode(estResult, expectedPhaseMode)
phaseMode = string(getDoaDopplerFieldOrDefault(estResult, 'phaseMode', ""));
if phaseMode ~= string(expectedPhaseMode)
  error('regressionMfCpVsIpRepresentative:InvalidPhaseMode', ...
    'Expected phaseMode=%s, got %s.', char(string(expectedPhaseMode)), char(phaseMode));
end
end

function aliasStepHz = localResolveAliasStep(timeOffsetSec)
timeOffsetSec = reshape(timeOffsetSec, [], 1);
if numel(timeOffsetSec) <= 1
  error('regressionMfCpVsIpRepresentative:InvalidTimeOffsets', ...
    'At least two frame offsets are required to resolve the alias step.');
end

deltaSec = diff(timeOffsetSec);
scaleUse = max(abs(timeOffsetSec));
if isempty(scaleUse) || ~isfinite(scaleUse) || scaleUse <= 0
  scaleUse = 1;
end
deltaSec = deltaSec(abs(deltaSec) > eps(scaleUse));
if isempty(deltaSec)
  error('regressionMfCpVsIpRepresentative:DegenerateTimeOffsets', ...
    'The representative frame offsets collapse to one repeated time.');
end
aliasStepHz = 1 / median(deltaSec);
end

function lineInfo = localScanAliasToothLine(model, truth, aliasStepHz)
numAliasSide = 2;
scanHalfWidthHz = min(60, 0.1 * aliasStepHz);
numGridPerTooth = 121;
probeOpt = struct('debugEnable', false);
basePoint = struct();
basePoint.latlon = truth.latlonTrueDeg(:);
basePoint.fdRate = truth.fdRateFit;
basePoint.fdRef = truth.fdRefFit;

aliasIdxList = (-numAliasSide:numAliasSide).';
minObjVec = nan(numel(aliasIdxList), 1);
minFdRefVec = nan(numel(aliasIdxList), 1);

for iAlias = 1:numel(aliasIdxList)
  aliasIdx = aliasIdxList(iAlias);
  centerFd = truth.fdRefFit + aliasIdx * aliasStepHz;
  fdGrid = linspace(centerFd - scanHalfWidthHz, centerFd + scanHalfWidthHz, numGridPerTooth);
  objVec = nan(numGridPerTooth, 1);
  for iGrid = 1:numGridPerTooth
    pointCur = basePoint;
    pointCur.fdRef = fdGrid(iGrid);
    probeEval = evalDoaDopplerMfProbePoint(model, pointCur, probeOpt);
    objVec(iGrid) = getDoaDopplerFieldOrDefault(probeEval, 'obj', NaN);
  end
  [minObjVec(iAlias), idxMin] = min(objVec);
  minFdRefVec(iAlias) = fdGrid(idxMin);
end

idxTruth = find(aliasIdxList == 0, 1, 'first');
truthObj = minObjVec(idxTruth);
firstAliasObj = min(minObjVec(abs(aliasIdxList) == 1));
secondAliasObj = min(minObjVec(abs(aliasIdxList) == 2));

lineInfo = struct();
lineInfo.aliasIdxList = aliasIdxList;
lineInfo.minObjVec = minObjVec;
lineInfo.minFdRefVec = minFdRefVec;
lineInfo.truthObj = truthObj;
lineInfo.firstAliasGap = firstAliasObj - truthObj;
lineInfo.secondAliasGap = secondAliasObj - truthObj;
end
