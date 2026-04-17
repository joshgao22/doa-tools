% Diagnostic regression for SF static MS weight replay solve paths.
% Contract:
%   On locked replay repeats that currently expose the MS static angle-loss
%   symptom, positive sat2 weights must not stay frozen on the MS-SF-DoA
%   anchor while fdRef still moves materially away from the W0 branch. If
%   that frozen-anchor pattern appears, the remaining bug is in the solve
%   path rather than in evaluator weight injection.
clear(); close all;

localAddProjectPath();
repeatIdxList = [22; 25];
snrDb = 10;
anchorTolDeg = 1e-8;
fdShiftTolHz = 5;
angleAwayFromW0TolDeg = 1e-4;

fprintf('Running regressionMsSfStaticSolverPathWeightReplay ...\n');
resultCell = cell(numel(repeatIdxList), 1);
flaggedRepeat = false(numel(repeatIdxList), 1);

for iRepeat = 1:numel(repeatIdxList)
  fixture = buildSfStaticSingleSnrRepeatFixture(repeatIdxList(iRepeat), snrDb);
  resultCell{iRepeat} = localAnalyzeRepeat(fixture, anchorTolDeg, fdShiftTolHz, angleAwayFromW0TolDeg);
  flaggedRepeat(iRepeat) = resultCell{iRepeat}.isFrozenAnchor;
end

pathDiagCell = cellfun(@(s) s.pathDiagTable, resultCell, 'UniformOutput', false);
pathDiagTable = vertcat(pathDiagCell{:});
summaryTable = table(repeatIdxList, cellfun(@(s) s.taskSeed, resultCell), ...
  cellfun(@(s) s.maxAngleToAnchorDeg, resultCell), ...
  cellfun(@(s) s.maxFdShiftVsW0Hz, resultCell), ...
  cellfun(@(s) s.minAngleToW0Deg, resultCell), flaggedRepeat, ...
  'VariableNames', {'repeatIdx', 'taskSeed', 'maxAngleToAnchorDeg', ...
  'maxFdShiftVsW0Hz', 'minAngleToW0Deg', 'isFrozenAnchor'});

fprintf('  replay repeats                 : %s\n', localFormatNumericRow(repeatIdxList));
fprintf('  path sensitivity summary       :\n');
disp(summaryTable);
fprintf('  positive-weight solve-path rows:\n');
disp(pathDiagTable);

if any(flaggedRepeat)
  badRepeat = repeatIdxList(flaggedRepeat);
  error('regressionMsSfStaticSolverPathWeightReplay:FrozenDoaAnchor', ...
    ['Positive sat2 weights stayed frozen on the MS-SF-DoA anchor on ', ...
     'replay repeat(s) %s while fdRef still moved away from W0.'], ...
    localFormatNumericRow(badRepeat));
end

fprintf('PASS: regressionMsSfStaticSolverPathWeightReplay\n');


function result = localAnalyzeRepeat(fixture, anchorTolDeg, fdShiftTolHz, angleAwayFromW0TolDeg)
%LOCALANALYZEREPEAT Analyze one replay repeat for frozen-anchor behavior.

alphaList = reshape(fixture.weightSweepAlpha, [], 1);
weightCase = reshape(fixture.caseBundle.weightCase, [], 1);
positiveIdx = find(alphaList > 0);
if isempty(positiveIdx)
  error('regressionMsSfStaticSolverPathWeightReplay:MissingPositiveWeight', ...
    'The replay repeat does not contain positive sat2 weights.');
end

msDoaAnchor = fixture.caseBundle.caseMsDoa.estResult.doaParamEst(:);
w0Case = weightCase(alphaList == 0);
if isempty(w0Case)
  error('regressionMsSfStaticSolverPathWeightReplay:MissingZeroWeightCase', ...
    'The replay repeat must contain the alpha=0 static MS case.');
end
w0Doa = w0Case(1).estResult.doaParamEst(:);
w0FdRef = w0Case(1).estResult.fdRefEst;

numPositive = numel(positiveIdx);
rowCell = cell(numPositive, 1);
angleToAnchorDeg = nan(numPositive, 1);
angleToW0Deg = nan(numPositive, 1);
fdShiftVsW0Hz = nan(numPositive, 1);

for iPos = 1:numPositive
  caseInfo = weightCase(positiveIdx(iPos));
  localAssertResolved(caseInfo, sprintf('repeat %d alpha=%.2f', ...
    fixture.repeatIdx, alphaList(positiveIdx(iPos))));
  estResult = caseInfo.estResult;
  doaNow = estResult.doaParamEst(:);
  angleToAnchorDeg(iPos) = calcLatlonAngleError(doaNow, msDoaAnchor);
  angleToW0Deg(iPos) = calcLatlonAngleError(doaNow, w0Doa);
  fdShiftVsW0Hz(iPos) = estResult.fdRefEst - w0FdRef;

  rowCell{iPos} = table( ...
    fixture.repeatIdx, fixture.taskSeed, alphaList(positiveIdx(iPos)), ...
    calcLatlonAngleError(doaNow, fixture.truth.latlonTrueDeg(:)), ...
    estResult.fdRefEst, fdShiftVsW0Hz(iPos), angleToAnchorDeg(iPos), ...
    angleToW0Deg(iPos), localGetOptimField(estResult, 'funcCount', NaN), ...
    localGetOptimField(estResult, 'iterations', NaN), ...
    'VariableNames', {'repeatIdx', 'taskSeed', 'alphaSat2', ...
    'angleErrDeg', 'fdRefEstHz', 'fdShiftVsW0Hz', ...
    'angleToMsDoaAnchorDeg', 'angleToW0Deg', 'funcCount', 'iterations'});
end

result = struct();
result.repeatIdx = fixture.repeatIdx;
result.taskSeed = fixture.taskSeed;
result.pathDiagTable = vertcat(rowCell{:});
result.maxAngleToAnchorDeg = max(angleToAnchorDeg);
result.maxFdShiftVsW0Hz = max(abs(fdShiftVsW0Hz));
result.minAngleToW0Deg = min(angleToW0Deg);
result.isFrozenAnchor = ...
  (result.maxAngleToAnchorDeg <= anchorTolDeg) && ...
  (result.maxFdShiftVsW0Hz > fdShiftTolHz) && ...
  (result.minAngleToW0Deg > angleAwayFromW0TolDeg);
end

function value = localGetOptimField(estResult, fieldName, defaultValue)
%LOCALGETOPTIMFIELD Read one field from estResult.optimInfo.

value = defaultValue;
optimInfo = struct();
if isstruct(estResult) && isfield(estResult, 'optimInfo')
  optimInfo = estResult.optimInfo;
end
if isstruct(optimInfo) && isfield(optimInfo, fieldName)
  value = optimInfo.(fieldName);
end
end

function localAssertResolved(caseInfo, displayName)
%LOCALASSERTRESOLVED Ensure one regression case resolved successfully.

if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  error('regressionMsSfStaticSolverPathWeightReplay:MissingCase', ...
    'Missing regression case: %s.', displayName);
end
estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  if ~logical(estResult.isResolved)
    error('regressionMsSfStaticSolverPathWeightReplay:UnresolvedCase', ...
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
