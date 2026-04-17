% Regression check for SF static sat-weight injection at one fixed probe point.
% Contract:
%   1) the per-satellite residual/noise terms must stay unchanged when only
%      satWeight changes at a fixed [DoA; fdRef] probe point;
%   2) the ref-sat objective contribution must stay unchanged;
%   3) the non-ref-sat objective contribution must scale linearly with alpha.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticRegressionFixture();

fprintf('Running regressionMsSfStaticWeightInjection ...\n');

alphaList = reshape(fixture.weightSweepAlpha, [], 1);
probePoint = struct();
probePoint.doaParam = fixture.caseBundle.caseStaticMs.estResult.doaParamEst(:);
probePoint.fdRef = fixture.caseBundle.caseStaticMs.estResult.fdRefEst;

numAlpha = numel(alphaList);
probeList = cell(numAlpha, 1);
residualEnergySat = nan(2, numAlpha);
noiseVarSat = nan(2, numAlpha);
objectiveSat = nan(2, numAlpha);

for iAlpha = 1:numAlpha
  alpha = alphaList(iAlpha);
  modelOpt = fixture.staticBaseOpt;
  modelOpt.satWeight = [1; alpha];
  probeUse = evalDoaDopplerSfProbePoint( ...
    fixture.viewMs, fixture.pilotWave, fixture.carrierFreq, fixture.sampleRate, ...
    fixture.fdRange, modelOpt, probePoint);
  probeList{iAlpha} = probeUse;
  residualEnergySat(:, iAlpha) = probeUse.aux.residualNormSat(:);
  noiseVarSat(:, iAlpha) = probeUse.aux.noiseVarSat(:);
  objectiveSat(:, iAlpha) = probeUse.aux.objectiveSat(:);
end

baseResidualEnergySat = residualEnergySat(:, 1);
baseNoiseVarSat = noiseVarSat(:, 1);
refObjectiveBase = objectiveSat(1, 1);

residualTol = 1e-9 * max(1, max(abs(baseResidualEnergySat)));
noiseTol = 1e-12 * max(1, max(abs(baseNoiseVarSat)));
objectiveTol = 1e-9 * max(1, max(abs(objectiveSat(:))));

residualEnergyDrift = residualEnergySat - repmat(baseResidualEnergySat, 1, numAlpha);
noiseVarDrift = noiseVarSat - repmat(baseNoiseVarSat, 1, numAlpha);
if max(abs(residualEnergyDrift(:))) > residualTol
  error('regressionMsSfStaticWeightInjection:ResidualDrift', ...
    'Changing satWeight at a fixed SF static probe point must not change per-satellite residual energy.');
end
if max(abs(noiseVarDrift(:))) > noiseTol
  error('regressionMsSfStaticWeightInjection:NoiseDrift', ...
    'Changing satWeight at a fixed SF static probe point must not change per-satellite noise estimates.');
end
if max(abs(objectiveSat(1, :) - refObjectiveBase)) > objectiveTol
  error('regressionMsSfStaticWeightInjection:RefObjectiveDrift', ...
    'The ref-sat SF static objective contribution must stay invariant when only sat2 weight changes.');
end

positiveIdx = find(alphaList > 0);
if isempty(positiveIdx)
  error('regressionMsSfStaticWeightInjection:MissingPositiveWeight', ...
    'The SF static weight probe requires at least one positive sat2 weight.');
end

otherObjectivePerUnit = objectiveSat(2, positiveIdx) ./ alphaList(positiveIdx).';
if max(abs(otherObjectivePerUnit - otherObjectivePerUnit(1))) > objectiveTol
  error('regressionMsSfStaticWeightInjection:OtherObjectiveScaling', ...
    'The non-ref-sat SF static objective contribution must scale linearly with sat2 weight.');
end
if abs(objectiveSat(2, 1)) > objectiveTol
  error('regressionMsSfStaticWeightInjection:ZeroWeightLeak', ...
    'The non-ref-sat SF static objective contribution must vanish at alpha = 0.');
end

fprintf('  alpha list                     : %s\n', localFormatNumericRow(alphaList));
fprintf('  residual energy sat1 sweep     : %s\n', localFormatNumericRow(residualEnergySat(1, :).')); 
fprintf('  residual energy sat2 sweep     : %s\n', localFormatNumericRow(residualEnergySat(2, :).'));
fprintf('  objective sat1 sweep           : %s\n', localFormatNumericRow(objectiveSat(1, :).'));
fprintf('  objective sat2 sweep           : %s\n', localFormatNumericRow(objectiveSat(2, :).'));
fprintf('  sat2 objective / alpha         : %s\n', ...
  localFormatNumericRow(otherObjectivePerUnit(:)));
fprintf('PASS: regressionMsSfStaticWeightInjection\n');


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
