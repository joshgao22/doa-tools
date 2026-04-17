% Regression check for SF static multi-sat objective decomposition.
% Contract:
%   1) the full SF static objective at weight [1; alpha] must equal the sum
%      of the ref-only full-scene objective and alpha times the non-ref-only
%      full-scene objective at the same fixed probe point;
%   2) the same decomposition must hold per-satellite in aux.objectiveSat.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticRegressionFixture();

fprintf('Running regressionMsSfStaticObjectiveDecomposition ...\n');

alphaList = reshape(fixture.weightSweepAlpha, [], 1);
probePoint = struct();
probePoint.doaParam = fixture.caseBundle.caseStaticMs.estResult.doaParamEst(:);
probePoint.fdRef = fixture.caseBundle.caseStaticMs.estResult.fdRefEst;

probeRefOnly = localEvalWeightCase(fixture, [1; 0], probePoint);
probeOtherOnly = localEvalWeightCase(fixture, [0; 1], probePoint);

numAlpha = numel(alphaList);
objErr = nan(numAlpha, 1);
perSatErr = nan(numAlpha, 2);
objFull = nan(numAlpha, 1);
objPred = nan(numAlpha, 1);

for iAlpha = 1:numAlpha
  alpha = alphaList(iAlpha);
  probeFull = localEvalWeightCase(fixture, [1; alpha], probePoint);

  objFull(iAlpha) = probeFull.obj;
  objPred(iAlpha) = probeRefOnly.obj + alpha * probeOtherOnly.obj;
  objErr(iAlpha) = objFull(iAlpha) - objPred(iAlpha);
  perSatErr(iAlpha, :) = probeFull.aux.objectiveSat(:).' - ...
    [probeRefOnly.aux.objectiveSat(1), alpha * probeOtherOnly.aux.objectiveSat(2)];
end

objTol = 1e-9 * max(1, max(abs([objFull; objPred])));
perSatScale = [perSatErr(:); probeRefOnly.aux.objectiveSat(:); probeOtherOnly.aux.objectiveSat(:)];
perSatTol = 1e-9 * max(1, max(abs(perSatScale)));
if max(abs(objErr)) > objTol
  error('regressionMsSfStaticObjectiveDecomposition:ObjectiveMismatch', ...
    'The full SF static objective must decompose into ref-only and non-ref-only weighted pieces.');
end
if max(abs(perSatErr(:))) > perSatTol
  error('regressionMsSfStaticObjectiveDecomposition:PerSatMismatch', ...
    'The per-satellite SF static objective pieces must follow the same weighted decomposition.');
end

fprintf('  alpha list                     : %s\n', localFormatNumericRow(alphaList));
fprintf('  full objective sweep           : %s\n', localFormatNumericRow(objFull));
fprintf('  predicted objective sweep      : %s\n', localFormatNumericRow(objPred));
fprintf('  objective decomposition error  : %s\n', localFormatNumericRow(objErr));
fprintf('  max per-sat error              : %.6g\n', max(abs(perSatErr(:))));
fprintf('PASS: regressionMsSfStaticObjectiveDecomposition\n');


function probe = localEvalWeightCase(fixture, satWeight, probePoint)
%LOCALEVALWEIGHTCASE Evaluate one SF static full-scene weight case.

modelOpt = fixture.staticBaseOpt;
modelOpt.satWeight = satWeight(:);
probe = evalDoaDopplerSfProbePoint( ...
  fixture.viewMs, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
  fixture.fdRange, modelOpt, probePoint);
end

function text = localFormatNumericRow(value)
%LOCALFORMATNUMERICROW Format one numeric vector for console output.

value = reshape(value, 1, []);
text = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%.6g', x), value, ...
  'UniformOutput', false), ', '));
end

function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(scriptDir)));
addpath(genpath(projectRoot));
end
