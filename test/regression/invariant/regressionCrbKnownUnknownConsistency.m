function regressionCrbKnownUnknownConsistency(varargin)
% Regression check for MF dynamic CRB known-vs-unknown consistency.
% Keep one direct guardrail on the paper-facing nuisance-rate claim:
%   1) the unified MF CRB entry must return stable parameter ordering for
%      known-rate and unknown-rate modes;
%   2) the unknown-rate CRB must not improve the DoA or fdRef lower bound
%      relative to the corresponding known-rate CRB.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
warnCleanup = localSuppressExpectedCrbWarnings(verbose); %#ok<NASGU>

context = buildDynamicDualSatEciContext(struct( ...
  'baseSeed', 253, ...
  'numSubsetRandomTrial', 0, ...
  'parallelOpt', struct('enableSubsetEvalParfor', false, 'minSubsetEvalParfor', 4)));
repeatData = buildDynamicRepeatData(context, 10, 253);
periodicFixture = repeatData.periodicFixture;
truth = periodicFixture.truth;

numSat = periodicFixture.sceneSeq.numSat;
numFrame = periodicFixture.sceneSeq.numFrame;
pathGain = ones(numSat, numFrame);
noiseVar = 0.1 * ones(numSat, numFrame);

knownOpt = struct();
knownOpt.doaType = 'latlon';
knownOpt.phaseMode = 'continuous';
knownOpt.fdRateMode = 'known';
knownOpt.steeringMode = 'framewise';
knownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;

unknownOpt = knownOpt;
unknownOpt.fdRateMode = 'unknown';

[crbKnown, auxKnown] = crbPilotMfDoaDoppler( ...
  periodicFixture.sceneSeq, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGain, noiseVar, knownOpt);
[crbUnknown, auxUnknown] = crbPilotMfDoaDoppler( ...
  periodicFixture.sceneSeq, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGain, noiseVar, unknownOpt);

localAssertCrbShape(crbKnown, auxKnown, 'known');
localAssertCrbShape(crbUnknown, auxUnknown, 'unknown');

paramNameKnown = string(auxKnown.paramName(:));
paramNameUnknown = string(auxUnknown.paramName(:));
if ~isequal(paramNameKnown, paramNameUnknown)
  error('regressionCrbKnownUnknownConsistency:ParamOrderMismatch', ...
    'Known/unknown MF CRB must keep the same interest-parameter ordering.');
end

if string(auxKnown.phaseMode) ~= "continuous" || string(auxUnknown.phaseMode) ~= "continuous"
  error('regressionCrbKnownUnknownConsistency:InvalidPhaseMode', ...
    'This regression expects continuous-phase CRBs for both known and unknown modes.');
end
if string(auxKnown.fdRateMode) ~= "known" || string(auxUnknown.fdRateMode) ~= "unknown"
  error('regressionCrbKnownUnknownConsistency:InvalidFdRateMode', ...
    'Known/unknown CRB outputs returned unexpected fdRateMode tags.');
end

knownDiag = max(real(diag(crbKnown)), 0);
unknownDiag = max(real(diag(crbUnknown)), 0);
knownAngleStdDeg = projectCrbToAngleMetric(crbKnown(1:2, 1:2), truth.latlonTrueDeg, 'latlon');
unknownAngleStdDeg = projectCrbToAngleMetric(crbUnknown(1:2, 1:2), truth.latlonTrueDeg, 'latlon');
knownFdStdHz = sqrt(knownDiag(3));
unknownFdStdHz = sqrt(unknownDiag(3));

  fprintf('Running regressionCrbKnownUnknownConsistency ...\n');
  fprintf('  known angle std (deg)         : %.6g\n', knownAngleStdDeg);
  fprintf('  unknown angle std (deg)       : %.6g\n', unknownAngleStdDeg);
  fprintf('  known fdRef std (Hz)          : %.6g\n', knownFdStdHz);
  fprintf('  unknown fdRef std (Hz)        : %.6g\n', unknownFdStdHz);
  fprintf('  known diag                    : [%g, %g, %g]\n', knownDiag(1), knownDiag(2), knownDiag(3));
  fprintf('  unknown diag                  : [%g, %g, %g]\n', unknownDiag(1), unknownDiag(2), unknownDiag(3));

tolAngle = max(1e-12, 1e-9 * max(unknownAngleStdDeg, 1));
tolFd = max(1e-12, 1e-9 * max(unknownFdStdHz, 1));
tolDiag = max(1e-12, 1e-9 * max(max(unknownDiag), 1));

if unknownAngleStdDeg + tolAngle < knownAngleStdDeg
  error('regressionCrbKnownUnknownConsistency:AngleBoundImprovedByUnknownRate', ...
    'Unknown-rate angle CRB must not improve over the known-rate CRB.');
end
if unknownFdStdHz + tolFd < knownFdStdHz
  error('regressionCrbKnownUnknownConsistency:FdBoundImprovedByUnknownRate', ...
    'Unknown-rate fdRef CRB must not improve over the known-rate CRB.');
end
if any(unknownDiag + tolDiag < knownDiag)
  error('regressionCrbKnownUnknownConsistency:DiagonalImprovedByUnknownRate', ...
    'Unknown-rate CRB must not improve any interest-parameter variance.');
end

  fprintf('PASS: regressionCrbKnownUnknownConsistency\n');


end


function cleanupObj = localSuppressExpectedCrbWarnings(verbose)
%LOCALSUPPRESSEXPECTEDCRBWARNINGS Suppress known full-FIM warnings in default regression output.

cleanupObj = [];
if verbose
  return;
end

warnState(1) = warning('off', 'crbPilotMfDoaDoppler:IllConditionedFim');
warnState(2) = warning('off', 'crbPilotMfDoaDoppler:IllConditionedFullFim');
cleanupObj = onCleanup(@() warning(warnState));
end

function localAssertCrbShape(crb, aux, modeTag)
if ~isequal(size(crb), [3, 3])
  error('regressionCrbKnownUnknownConsistency:InvalidCrbSize', ...
    '%s MF CRB must be 3x3.', char(string(modeTag)));
end
if any(~isfinite(crb(:)))
  error('regressionCrbKnownUnknownConsistency:NonFiniteCrb', ...
    '%s MF CRB must be finite.', char(string(modeTag)));
end
paramName = string(getDoaDopplerFieldOrDefault(aux, 'paramName', strings(0, 1)));
if numel(paramName) ~= 3
  error('regressionCrbKnownUnknownConsistency:InvalidParamNameList', ...
    '%s MF CRB must expose exactly three interest-parameter names.', ...
    char(string(modeTag)));
end
end
