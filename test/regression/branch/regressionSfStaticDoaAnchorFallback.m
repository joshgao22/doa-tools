function regressionSfStaticDoaAnchorFallback(varargin)
% Regression check for the SF static DoA-anchor fallback candidate.
% This script verifies one branch-level contract only:
%   1) the multi-sat static estimator may compare one fixed-DoA anchor
%      candidate against the free joint solve;
%   2) when the objective gap stays within the configured tolerance, the
%      selected candidate should not degrade angle or fdRef relative to the
%      pure free solve.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;
localPrint = @(varargin) fprintf(varargin{:});
fixture = buildSfStaticRegressionFixture();

localPrint('Running regressionSfStaticDoaAnchorFallback ...\n');

optFree = fixture.caseBundle.staticMsOpt;
optFree.initDoaHalfWidth = [];
optFree.enableDoaAnchorFallback = false;
[estFree, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  fixture.viewMs.sceneRef, fixture.viewMs.rxSigSf, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, 1, [], false, optFree);

optFallback = fixture.caseBundle.staticMsOpt;
optFallback.initDoaHalfWidth = [];
optFallback.enableDoaAnchorFallback = true;
[estFallback, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  fixture.viewMs.sceneRef, fixture.viewMs.rxSigSf, fixture.pilotWave, ...
  fixture.carrierFreq, fixture.sampleRate, fixture.viewMs.doaGrid, ...
  fixture.fdRange, 1, [], false, optFallback);

freeAngleErr = calcLatlonAngleError(estFree.doaParamEst(:), fixture.truth.latlonTrueDeg(:));
fallbackAngleErr = calcLatlonAngleError(estFallback.doaParamEst(:), fixture.truth.latlonTrueDeg(:));
freeFdErr = estFree.fdRefEst - fixture.truth.fdRefTrueHz;
fallbackFdErr = estFallback.fdRefEst - fixture.truth.fdRefTrueHz;
selectedTag = string(getDoaDopplerFieldOrDefault(estFallback.optimInfo, 'selectedCandidateTag', ""));
objGap = getDoaDopplerFieldOrDefault(estFallback.optimInfo, 'doaAnchorObjGap', NaN);
objTol = getDoaDopplerFieldOrDefault(estFallback.optimInfo, 'doaAnchorObjTol', NaN);

if selectedTag == "fixedDoaAnchor"
  if fallbackAngleErr > freeAngleErr + 1e-12
    error('regressionSfStaticDoaAnchorFallback:AngleDegraded', ...
      'The selected fixed-DoA anchor candidate must not degrade the angle error relative to the free solve.');
  end
  if abs(fallbackFdErr) > abs(freeFdErr) + 1e-9
    error('regressionSfStaticDoaAnchorFallback:FdRefDegraded', ...
      'The selected fixed-DoA anchor candidate must not degrade the fdRef error relative to the free solve.');
  end
end

if ~all(isfinite([freeAngleErr, fallbackAngleErr, freeFdErr, fallbackFdErr]))
  error('regressionSfStaticDoaAnchorFallback:NonFiniteMetric', ...
    'The fallback comparison metrics must stay finite.');
end

localPrint('  selected candidate tag        : %s\n', selectedTag);
localPrint('  free angle err (deg)          : %.6g\n', freeAngleErr);
localPrint('  fallback angle err (deg)      : %.6g\n', fallbackAngleErr);
localPrint('  free fdRef err (Hz)           : %.6f\n', freeFdErr);
localPrint('  fallback fdRef err (Hz)       : %.6f\n', fallbackFdErr);
localPrint('  fallback obj gap              : %.6g\n', objGap);
localPrint('  fallback obj tol              : %.6g\n', objTol);
localPrint('PASS: regressionSfStaticDoaAnchorFallback\n');
end
