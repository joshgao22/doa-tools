% Diagnostic regression for SF static nonzero sat2-weight DoA sensitivity.
% Contract:
%   Positive sat2 weights must not collapse to an essentially identical DoA
%   solution while fdRef continues to move materially across the same weight
%   sweep. If that pattern appears, it is a strong sign that sat2 weight is
%   not entering the DoA optimization path correctly.
%
% Note:
%   This regression is intentionally diagnostic. It is expected to fail on
%   the current bugged implementation and should pass only after the
%   nonzero-weight DoA coupling is repaired.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticJointCouplingFixture();

fprintf('Running regressionMsSfStaticNonzeroWeightDoaSensitivity ...\n');

alphaList = reshape(fixture.weightSweepAlpha, [], 1);
weightCase = reshape(fixture.caseBundle.weightCase, [], 1);
positiveIdx = find(alphaList > 0);
if numel(positiveIdx) < 2
  error('regressionMsSfStaticNonzeroWeightDoaSensitivity:NeedMultiplePositiveWeights', ...
    'The diagnostic regression needs at least two positive sat2 weights.');
end

numPositive = numel(positiveIdx);
doaParamMat = nan(2, numPositive);
fdRefList = nan(numPositive, 1);
angleErrList = nan(numPositive, 1);
truthLatlon = fixture.truth.latlonTrueDeg(:);

for iPos = 1:numPositive
  caseInfo = weightCase(positiveIdx(iPos));
  localAssertResolved(caseInfo, sprintf('alpha=%.2f', alphaList(positiveIdx(iPos))));
  doaParamMat(:, iPos) = caseInfo.estResult.doaParamEst(:);
  fdRefList(iPos) = caseInfo.estResult.fdRefEst;
  angleErrList(iPos) = calcLatlonAngleError(doaParamMat(:, iPos), truthLatlon);
end

pairIdx = nchoosek(1:numPositive, 2);
numPair = size(pairIdx, 1);
doaDiffDeg = nan(numPair, 1);
fdSpreadHz = nan(numPair, 1);
alphaPairMat = nan(numPair, 2);
for iPair = 1:numPair
  idxA = pairIdx(iPair, 1);
  idxB = pairIdx(iPair, 2);
  doaDiffDeg(iPair) = calcLatlonAngleError(doaParamMat(:, idxA), doaParamMat(:, idxB));
  fdSpreadHz(iPair) = abs(fdRefList(idxA) - fdRefList(idxB));
  alphaPairMat(iPair, :) = alphaList(positiveIdx([idxA, idxB])).';
end

maxDoaDiffDeg = max(doaDiffDeg);
maxFdSpreadHz = max(fdSpreadHz);
doaCollapseTolDeg = 1e-6;
fdSpreadTolHz = 5;

if (maxDoaDiffDeg <= doaCollapseTolDeg) && (maxFdSpreadHz > fdSpreadTolHz)
  [~, worstIdx] = max(fdSpreadHz);
  error('regressionMsSfStaticNonzeroWeightDoaSensitivity:CollapsedDoaPath', ...
    ['Positive sat2 weights collapsed to the same DoA (max pairwise diff %.6g deg) ', ...
     'while fdRef still moved by %.6g Hz between alpha=%.2f and alpha=%.2f.'], ...
    maxDoaDiffDeg, maxFdSpreadHz, alphaPairMat(worstIdx, 1), alphaPairMat(worstIdx, 2));
end

fprintf('  positive alpha list            : %s\n', localFormatNumericRow(alphaList(positiveIdx)));
fprintf('  positive-weight angle err (deg): %s\n', localFormatNumericRow(angleErrList));
fprintf('  positive-weight fdRef (Hz)     : %s\n', localFormatNumericRow(fdRefList));
fprintf('  max pairwise DoA diff (deg)    : %.6g\n', maxDoaDiffDeg);
fprintf('  max pairwise fdRef spread (Hz) : %.6g\n', maxFdSpreadHz);
fprintf('PASS: regressionMsSfStaticNonzeroWeightDoaSensitivity\n');


function localAssertResolved(caseInfo, displayName)
%LOCALASSERTRESOLVED Ensure one regression case resolved successfully.

if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  error('regressionMsSfStaticNonzeroWeightDoaSensitivity:MissingCase', ...
    'Missing regression case: %s.', displayName);
end
estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  if ~logical(estResult.isResolved)
    error('regressionMsSfStaticNonzeroWeightDoaSensitivity:UnresolvedCase', ...
      'Regression case did not resolve: %s.', displayName);
  end
end
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
