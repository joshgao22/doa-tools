function regressionSfStaticJointCouplingCaseLocked(varargin)
% Regression check for the locked SF static SS/MS joint-coupling case.
% This script mirrors the current doaDopplerStatDualSatUraEci setup and
% keeps the contract intentionally narrow:
%   1) the locked SS/MS static cases must remain reproducible and resolved;
%   2) free, truth, and mixed-point objective re-evaluations must be
%      internally consistent;
%   3) the printed mixed-point diagnostics should help decide whether the
%      remaining MS static gap is dominated by DoA coupling, fdRef coupling,
%      or both.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
localPrint = @(varargin) fprintf(varargin{:});
fixture = buildSfStaticJointCouplingFixture();

localPrint('Running regressionSfStaticJointCouplingCaseLocked ...\n');

caseBundle = fixture.caseBundle;
caseSsStatic = caseBundle.caseStaticRefOnly;
caseMsStatic = caseBundle.caseStaticMs;

if ~localCaseResolved(caseSsStatic) || ~localCaseResolved(caseMsStatic)
  error('regressionSfStaticJointCouplingCaseLocked:UnresolvedStaticCase', ...
    'The locked SS/MS static cases must both resolve before coupling checks.');
end

truthDoa = fixture.truth.latlonTrueDeg(:);
truthFdRef = fixture.truth.fdRefTrueHz;
freeSsDoa = caseSsStatic.estResult.doaParamEst(:);
freeMsDoa = caseMsStatic.estResult.doaParamEst(:);
freeSsFdRef = caseSsStatic.estResult.fdRefEst;
freeMsFdRef = caseMsStatic.estResult.fdRefEst;

freeSsEval = localEvalStaticPoint(fixture.viewRefOnly, fixture, freeSsDoa, freeSsFdRef, 1);
freeMsEval = localEvalStaticPoint(fixture.viewMs, fixture, freeMsDoa, freeMsFdRef, []);
truthSsEval = localEvalStaticPoint(fixture.viewRefOnly, fixture, truthDoa, truthFdRef, 1);
truthMsEval = localEvalStaticPoint(fixture.viewMs, fixture, truthDoa, truthFdRef, []);

mixMsTruthDoaEval = localEvalStaticPoint(fixture.viewMs, fixture, truthDoa, freeMsFdRef, []);
mixMsTruthFdEval = localEvalStaticPoint(fixture.viewMs, fixture, freeMsDoa, truthFdRef, []);
mixSsTruthDoaEval = localEvalStaticPoint(fixture.viewRefOnly, fixture, truthDoa, freeSsFdRef, 1);
mixSsTruthFdEval = localEvalStaticPoint(fixture.viewRefOnly, fixture, freeSsDoa, truthFdRef, 1);

localAssertConsistentFval(caseSsStatic, freeSsEval, 'SS');
localAssertConsistentFval(caseMsStatic, freeMsEval, 'MS');
localAssertFiniteEval(truthSsEval, 'truthSs');
localAssertFiniteEval(truthMsEval, 'truthMs');
localAssertFiniteEval(mixMsTruthDoaEval, 'mixMsTruthDoa');
localAssertFiniteEval(mixMsTruthFdEval, 'mixMsTruthFd');
localAssertFiniteEval(mixSsTruthDoaEval, 'mixSsTruthDoa');
localAssertFiniteEval(mixSsTruthFdEval, 'mixSsTruthFd');

freeSsAngleErrDeg = calcLatlonAngleError(freeSsDoa, truthDoa);
freeMsAngleErrDeg = calcLatlonAngleError(freeMsDoa, truthDoa);
freeSsFdErrHz = freeSsFdRef - truthFdRef;
freeMsFdErrHz = freeMsFdRef - truthFdRef;

msTruthDoaImprove = freeMsEval.fval - mixMsTruthDoaEval.fval;
msTruthFdImprove = freeMsEval.fval - mixMsTruthFdEval.fval;
ssTruthDoaImprove = freeSsEval.fval - mixSsTruthDoaEval.fval;
ssTruthFdImprove = freeSsEval.fval - mixSsTruthFdEval.fval;

msDominantAxis = localClassifyDominantAxis(msTruthDoaImprove, msTruthFdImprove);
ssDominantAxis = localClassifyDominantAxis(ssTruthDoaImprove, ssTruthFdImprove);

localPrint('  free SS angle err (deg)        : %.6g\n', freeSsAngleErrDeg);
localPrint('  free MS angle err (deg)        : %.6g\n', freeMsAngleErrDeg);
localPrint('  free SS fdRef err (Hz)         : %.6f\n', freeSsFdErrHz);
localPrint('  free MS fdRef err (Hz)         : %.6f\n', freeMsFdErrHz);
localPrint('  free SS objective              : %.6g\n', freeSsEval.fval);
localPrint('  free MS objective              : %.6g\n', freeMsEval.fval);
localPrint('  truth SS objective             : %.6g\n', truthSsEval.fval);
localPrint('  truth MS objective             : %.6g\n', truthMsEval.fval);
localPrint('  mix MS truthDoA/freeFd obj     : %.6g\n', mixMsTruthDoaEval.fval);
localPrint('  mix MS freeDoA/truthFd obj     : %.6g\n', mixMsTruthFdEval.fval);
localPrint('  mix SS truthDoA/freeFd obj     : %.6g\n', mixSsTruthDoaEval.fval);
localPrint('  mix SS freeDoA/truthFd obj     : %.6g\n', mixSsTruthFdEval.fval);
localPrint('  MS truthDoA improve            : %.6g\n', msTruthDoaImprove);
localPrint('  MS truthFd improve             : %.6g\n', msTruthFdImprove);
localPrint('  SS truthDoA improve            : %.6g\n', ssTruthDoaImprove);
localPrint('  SS truthFd improve             : %.6g\n', ssTruthFdImprove);
localPrint('  MS dominant coupling axis      : %s\n', msDominantAxis);
localPrint('  SS dominant coupling axis      : %s\n', ssDominantAxis);

localPrint('  free MS objectiveSat           : %s\n', localFormatNumericRow(freeMsEval.objectiveSat));
localPrint('  truth MS objectiveSat          : %s\n', localFormatNumericRow(truthMsEval.objectiveSat));
localPrint('  mix MS truthDoA objectiveSat   : %s\n', localFormatNumericRow(mixMsTruthDoaEval.objectiveSat));
localPrint('  mix MS truthFd objectiveSat    : %s\n', localFormatNumericRow(mixMsTruthFdEval.objectiveSat));

localPrint('PASS: regressionSfStaticJointCouplingCaseLocked\n');


end

function evalInfo = localEvalStaticPoint(view, fixture, doaParam, fdRefHz, satWeight)
%LOCALEVALSTATICPOINT Evaluate one locked SF static point through the estimator.

modelOpt = fixture.staticBaseOpt;
modelOpt.evalOnly = true;
if nargin >= 5 && ~isempty(satWeight)
  modelOpt.satWeight = satWeight;
end
initParam = [reshape(doaParam, [], 1); fdRefHz];
[estResult, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  view.sceneRef, view.rxSigSf, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, view.doaGrid, fixture.fdRange, 1, initParam, false, modelOpt);

evalInfo = struct();
evalInfo.fval = estResult.fval;
evalInfo.objectiveSat = reshape(getDoaDopplerFieldOrDefault(estResult.aux, 'objectiveSat', NaN), [], 1);
evalInfo.residualNormSat = reshape(getDoaDopplerFieldOrDefault(estResult.aux, 'residualNormSat', NaN), [], 1);
evalInfo.fdSatEstHz = reshape(getDoaDopplerFieldOrDefault(estResult.aux, 'fdSatEst', NaN), [], 1);
evalInfo.deltaFdRefEstHz = reshape(getDoaDopplerFieldOrDefault(estResult.aux, 'deltaFdRefEst', NaN), [], 1);
evalInfo.evalResult = estResult;
end


function localAssertConsistentFval(caseInfo, evalInfo, caseTag)
%LOCALASSERTCONSISTENTFVAL Require eval-only replay to match stored fval.

storedFval = getDoaDopplerFieldOrDefault(caseInfo.estResult, 'fval', NaN);
if ~isfinite(storedFval) || ~isfinite(evalInfo.fval)
  error('regressionSfStaticJointCouplingCaseLocked:NonFiniteReplay', ...
    '%s static replay returned a non-finite objective.', caseTag);
end

fvalTol = 1e-8 * max(1, abs(storedFval));
if abs(evalInfo.fval - storedFval) > fvalTol
  error('regressionSfStaticJointCouplingCaseLocked:ReplayMismatch', ...
    '%s static replay must match the stored free-case objective.', caseTag);
end
end


function localAssertFiniteEval(evalInfo, evalTag)
%LOCALASSERTFINITEEVAL Require one mixed-point evaluation to stay finite.

if ~isfinite(evalInfo.fval)
  error('regressionSfStaticJointCouplingCaseLocked:NonFiniteObjective', ...
    'The %s objective must remain finite.', evalTag);
end
if isempty(evalInfo.objectiveSat) || any(~isfinite(evalInfo.objectiveSat))
  error('regressionSfStaticJointCouplingCaseLocked:NonFiniteObjectiveSat', ...
    'The %s per-satellite objective slices must remain finite.', evalTag);
end
end


function label = localClassifyDominantAxis(doaImprove, fdImprove)
%LOCALCLASSIFYDOMINANTAXIS Label the stronger mixed-point improvement axis.

tol = 1e-10 * max(1, abs(doaImprove) + abs(fdImprove));
if doaImprove > fdImprove + tol
  label = "doa-dominant";
elseif fdImprove > doaImprove + tol
  label = "fd-dominant";
else
  label = "balanced";
end
end


function tf = localCaseResolved(caseInfo)
%LOCALCASERESOLVED Return true when one case contains a usable estimate.

tf = false;
if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end

estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  tf = logical(estResult.isResolved);
else
  tf = true;
end
if ~tf
  return;
end

tf = isfield(estResult, 'doaParamEst') && numel(estResult.doaParamEst) == 2 && ...
  all(isfinite(estResult.doaParamEst(:))) && isfield(estResult, 'fdRefEst') && ...
  isscalar(estResult.fdRefEst) && isfinite(estResult.fdRefEst);
end


function text = localFormatNumericRow(value)
%LOCALFORMATNUMERICROW Format one numeric vector for compact printing.

value = reshape(value, 1, []);
text = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%.6g', x), value, ...
  'UniformOutput', false), ', '));
end
