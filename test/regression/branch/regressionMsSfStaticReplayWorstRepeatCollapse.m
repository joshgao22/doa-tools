% Diagnostic regression for the worst fixed-SNR MS static angle-loss repeat.
% Contract:
%   On the locked replay of repeat 25 from doaDopplerStatDualSatUraEci,
%   positive sat2 weights must not collapse to the same DoA while fdRef
%   still spreads materially across the sweep. If that collapse appears,
%   the solve path is still stuck on the known MS static branch bug.
clear(); close all;

localAddProjectPath();
fixture = buildSfStaticSingleSnrRepeatFixture(25, 10);

fprintf('Running regressionMsSfStaticReplayWorstRepeatCollapse ...\n');
fprintf('  replay repeat / taskSeed        : %d / %d\n', ...
  fixture.repeatIdx, fixture.taskSeed);

alphaList = reshape(fixture.weightSweepAlpha, [], 1);
weightCase = reshape(fixture.caseBundle.weightCase, [], 1);
positiveIdx = find(alphaList > 0);
if numel(positiveIdx) < 2
  error('regressionMsSfStaticReplayWorstRepeatCollapse:NeedMultiplePositiveWeights', ...
    'The replay regression needs at least two positive sat2 weights.');
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
doaCollapseTolDeg = 1e-8;
fdSpreadTolHz = 5;

fprintf('  positive alpha list            : %s\n', localFormatNumericRow(alphaList(positiveIdx)));
fprintf('  positive-weight angle err (deg): %s\n', localFormatNumericRow(angleErrList));
fprintf('  positive-weight fdRef (Hz)     : %s\n', localFormatNumericRow(fdRefList));
fprintf('  max pairwise DoA diff (deg)    : %.6g\n', maxDoaDiffDeg);
fprintf('  max pairwise fdRef spread (Hz) : %.6g\n', maxFdSpreadHz);

if (maxDoaDiffDeg <= doaCollapseTolDeg) && (maxFdSpreadHz > fdSpreadTolHz)
  [~, worstIdx] = max(fdSpreadHz);
  error('regressionMsSfStaticReplayWorstRepeatCollapse:CollapsedPositiveWeightBranch', ...
    ['Replay repeat %d collapsed all positive sat2 weights to the same DoA ', ...
     '(max pairwise diff %.6g deg) while fdRef still moved by %.6g Hz ', ...
     'between alpha=%.2f and alpha=%.2f.'], ...
    fixture.repeatIdx, maxDoaDiffDeg, maxFdSpreadHz, ...
    alphaPairMat(worstIdx, 1), alphaPairMat(worstIdx, 2));
end

fprintf('PASS: regressionMsSfStaticReplayWorstRepeatCollapse\n');


function localAssertResolved(caseInfo, displayName)
%LOCALASSERTRESOLVED Ensure one regression case resolved successfully.

if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  error('regressionMsSfStaticReplayWorstRepeatCollapse:MissingCase', ...
    'Missing regression case: %s.', displayName);
end
estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  if ~logical(estResult.isResolved)
    error('regressionMsSfStaticReplayWorstRepeatCollapse:UnresolvedCase', ...
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
