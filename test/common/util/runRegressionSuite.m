function result = runRegressionSuite(suiteName, caseRelPathList, varargin)
%RUNREGRESSIONSUITE Run one grouped regression suite.
%
% Syntax:
%   result = runRegressionSuite(suiteName, caseRelPathList)
%   result = runRegressionSuite(suiteName, caseRelPathList, runOpt)
%   result = runRegressionSuite(suiteName, caseRelPathList, Name, Value)
%
% Inputs:
%   suiteName       - suite tag shown in summary messages.
%   caseRelPathList - regression case paths relative to project root.
%
% Optional inputs:
%   runOpt.projectRoot   - repository root. Default: auto resolve.
%   runOpt.stopOnFailure - stop on first failure. Default: false.
%   runOpt.verbose       - enable case-level diagnostics. Default: false.
%
% Output:
%   result          - suite result structure.
%
% Notes:
%   - Regression cases are function based and are executed with feval.
%   - Suite progress and pass/fail summary are always printed.
%   - Each case receives the shared verbose option through
%     feval(caseName, 'verbose', runOpt.verbose).
%   - This helper only orchestrates grouped execution. Regression-specific
%     assertions stay inside each individual regression case.

if ~(ischar(suiteName) || (isstring(suiteName) && isscalar(suiteName)))
  error('runRegressionSuite:InvalidSuiteName', ...
    'suiteName must be a character vector or string scalar.');
end

caseRelPathList = string(caseRelPathList);
caseRelPathList = caseRelPathList(:);
if isempty(caseRelPathList)
  error('runRegressionSuite:EmptyCaseList', ...
    'caseRelPathList must contain at least one regression case.');
end

runOpt = localParseRunOpt(varargin{:});
projectRoot = localResolveProjectRoot(runOpt.projectRoot);
addpath(genpath(projectRoot));

numCase = numel(caseRelPathList);
item = repmat(localBuildEmptyItem(), numCase, 1);
totalRunTimeSec = 0;

verboseState = localSnapshotVerboseState();
cleanupVerbose = onCleanup(@() localRestoreVerboseState(verboseState));
localAssignVerboseFlag(runOpt.verbose);

fprintf('=== Running %s regression suite (%d cases) ===\n', ...
  char(string(suiteName)), numCase);

for iCase = 1:numCase
  caseRelPath = char(caseRelPathList(iCase));
  casePath = fullfile(projectRoot, caseRelPath);
  caseName = localResolveCaseName(caseRelPath);

  itemUse = localBuildEmptyItem();
  itemUse.caseName = string(caseName);
  itemUse.caseRelPath = string(caseRelPath);
  itemUse.casePath = string(casePath);
  itemUse.scriptName = itemUse.caseName;
  itemUse.scriptRelPath = itemUse.caseRelPath;
  itemUse.scriptPath = itemUse.casePath;

  fprintf('[%d/%d] %s\n', iCase, numCase, char(itemUse.caseName));

  if exist(casePath, 'file') ~= 2
    itemUse.status = "FAIL";
    itemUse.errorId = "runRegressionSuite:CaseNotFound";
    itemUse.errorMessage = sprintf('Regression case not found: %s', caseRelPath);
    item(iCase) = itemUse;

    fprintf('  FAIL (missing file)\n');

    if runOpt.stopOnFailure
      localPrintSuiteSummary(suiteName, item(1:iCase));
      error('runRegressionSuite:CaseNotFound', '%s', itemUse.errorMessage);
    end
    continue;
  end

  runTimer = tic;
  try
    localAssignVerboseFlag(runOpt.verbose);
    localRunRegressionCase(caseName, runOpt.verbose);
    itemUse.status = "PASS";
  catch ME
    itemUse.status = "FAIL";
    itemUse.errorId = string(localNormalizeErrorId(ME.identifier));
    itemUse.errorMessage = string(ME.message);
    itemUse.errorLocation = string(localBuildErrorLocation(ME));
  end

  itemUse.runTimeSec = toc(runTimer);
  totalRunTimeSec = totalRunTimeSec + itemUse.runTimeSec;
  item(iCase) = itemUse;

  fprintf('  %s (%.3f s)\n', char(itemUse.status), itemUse.runTimeSec);
  if itemUse.status == "FAIL"
    if strlength(itemUse.errorId) > 0
      fprintf('    %s\n', char(itemUse.errorId));
    end
    fprintf('    %s\n', char(itemUse.errorMessage));
  end

  if itemUse.status == "FAIL" && runOpt.stopOnFailure
    localPrintSuiteSummary(suiteName, item(1:iCase));
    error('runRegressionSuite:RegressionFailed', ...
      'Regression failed in %s: %s\n%s', ...
      char(string(suiteName)), char(itemUse.caseName), char(itemUse.errorMessage));
  end
end

result = struct();
result.suiteName = string(suiteName);
result.projectRoot = string(projectRoot);
result.caseRelPathList = caseRelPathList;
result.scriptRelPathList = caseRelPathList; % Compatibility alias for older callers.
result.item = item;
result.summaryTable = localBuildSummaryTable(item);
result.numCase = numCase;
result.numScript = numCase; % Compatibility alias for older callers.
result.numPassed = nnz(string({item.status}) == "PASS");
result.numFailed = nnz(string({item.status}) == "FAIL");
result.allPassed = (result.numFailed == 0);
result.totalRunTimeSec = totalRunTimeSec;

localPrintSuiteSummary(suiteName, item);
if ~result.allPassed
  localPrintFailureDetail(item);
end

if ~result.allPassed && ~runOpt.stopOnFailure
  firstFailIdx = find(string({item.status}) == "FAIL", 1, 'first');
  firstFail = item(firstFailIdx);
  error('runRegressionSuite:RegressionFailed', ...
    'Regression suite %s finished with %d failure(s). First failure: %s\n%s', ...
    char(string(suiteName)), result.numFailed, char(firstFail.caseName), char(firstFail.errorMessage));
end
end


function runOpt = localParseRunOpt(varargin)
%LOCALPARSERUNOPT Parse runner options from struct or Name/Value pairs.

runOpt = struct();
runOpt.projectRoot = "";
runOpt.stopOnFailure = false;
runOpt.verbose = false;

if nargin == 0
  return;
end

if nargin == 1 && isstruct(varargin{1})
  optIn = varargin{1};
  if isfield(optIn, 'projectRoot') && ~isempty(optIn.projectRoot)
    runOpt.projectRoot = string(optIn.projectRoot);
  end
  if isfield(optIn, 'stopOnFailure') && ~isempty(optIn.stopOnFailure)
    runOpt.stopOnFailure = logical(optIn.stopOnFailure);
  end
  if isfield(optIn, 'verbose') && ~isempty(optIn.verbose)
    runOpt.verbose = logical(optIn.verbose);
  end
  if isfield(optIn, 'stopOnError') && ~isempty(optIn.stopOnError)
    runOpt.stopOnFailure = logical(optIn.stopOnError);
  end
  return;
end

if mod(nargin, 2) ~= 0
  error('runRegressionSuite:InvalidOptionInput', ...
    'Optional inputs must be one struct or Name/Value pairs.');
end

for iArg = 1:2:nargin
  name = lower(char(string(varargin{iArg})));
  value = varargin{iArg + 1};
  switch name
    case 'projectroot'
      runOpt.projectRoot = string(value);
    case 'stoponfailure'
      runOpt.stopOnFailure = logical(value);
    case 'verbose'
      runOpt.verbose = logical(value);
    case {'stoponerror', 'stoponfail'}
      runOpt.stopOnFailure = logical(value);
    otherwise
      error('runRegressionSuite:UnknownOption', ...
        'Unknown option name: %s', char(string(varargin{iArg})));
  end
end
end


function item = localBuildEmptyItem()
%LOCALBUILDEMPTYITEM Build one empty regression result record.

item = struct();
item.caseName = "";
item.caseRelPath = "";
item.casePath = "";
item.scriptName = "";
item.scriptRelPath = "";
item.scriptPath = "";
item.status = "NOT_RUN";
item.runTimeSec = NaN;
item.errorId = "";
item.errorMessage = "";
item.errorLocation = "";
end


function projectRoot = localResolveProjectRoot(projectRootIn)
%LOCALRESOLVEPROJECTROOT Resolve project root from input or file path.

if strlength(string(projectRootIn)) > 0
  projectRoot = char(string(projectRootIn));
  return;
end

thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(thisDir)));
end


function caseName = localResolveCaseName(caseRelPath)
%LOCALRESOLVECASEName Extract case base name without extension.

[~, caseName, ~] = fileparts(caseRelPath);
end


function localRunRegressionCase(caseName, verbose)
%LOCALRUNREGRESSIONCASE Run one function-style regression case.

feval(caseName, 'verbose', verbose);
end


function summaryTable = localBuildSummaryTable(item)
%LOCALBUILDSUMMARYTABLE Convert item array to a compact table.

numItem = numel(item);
caseName = strings(numItem, 1);
scriptName = strings(numItem, 1);
status = strings(numItem, 1);
runTimeSec = NaN(numItem, 1);
errorId = strings(numItem, 1);
errorLocation = strings(numItem, 1);

for iItem = 1:numItem
  caseName(iItem) = item(iItem).caseName;
  scriptName(iItem) = item(iItem).scriptName;
  status(iItem) = item(iItem).status;
  runTimeSec(iItem) = item(iItem).runTimeSec;
  errorId(iItem) = item(iItem).errorId;
  errorLocation(iItem) = item(iItem).errorLocation;
end

summaryTable = table(caseName, scriptName, status, runTimeSec, errorId, errorLocation);
end


function localPrintSuiteSummary(suiteName, item)
%LOCALPRINTSUITESUMMARY Print one compact suite summary.

numPassed = nnz(string({item.status}) == "PASS");
numFailed = nnz(string({item.status}) == "FAIL");
numTotal = numel(item);
totalRunTimeSec = sum([item.runTimeSec], 'omitnan');

fprintf('=== %s regression summary: %d/%d passed, %d failed, total %.3f s ===\n', ...
  char(string(suiteName)), numPassed, numTotal, numFailed, totalRunTimeSec);
end


function localPrintFailureDetail(item)
%LOCALPRINTFAILUREDETAIL Print one compact failure list after the suite finishes.

failIdx = find(string({item.status}) == "FAIL");
if isempty(failIdx)
  return;
end
fprintf('--- failed regressions ---\n');
for iFail = reshape(failIdx, 1, [])
  fprintf('  %s', char(item(iFail).caseName));
  if strlength(item(iFail).errorId) > 0
    fprintf(' [%s]', char(item(iFail).errorId));
  end
  fprintf('\n');
  if strlength(item(iFail).errorLocation) > 0
    fprintf('    at %s\n', char(item(iFail).errorLocation));
  end
  if strlength(item(iFail).errorMessage) > 0
    fprintf('    %s\n', char(item(iFail).errorMessage));
  end
end
end


function errorLocation = localBuildErrorLocation(ME)
%LOCALBUILDERRORLOCATION Build one compact first-stack error location string.

errorLocation = "";
if isempty(ME.stack)
  return;
end
frame = ME.stack(1);
errorLocation = string(frame.name) + ":" + string(frame.line);
end


function errorId = localNormalizeErrorId(errorIdIn)
%LOCALNORMALIZEERRORID Normalize empty error ids.

errorId = char(string(errorIdIn));
if isempty(errorId)
  errorId = 'runRegressionSuite:UnknownError';
end
end

function verboseState = localSnapshotVerboseState()
%LOCALSNAPSHOTVERBOSESTATE Snapshot inherited verbose carriers.

verboseState = struct();
verboseState.baseExists = false;
verboseState.baseValue = false;
verboseState.appDataExists = false;
verboseState.appDataValue = false;

try
  hasVerbose = evalin('base', "exist(''doaToolsVerbose'', ''var'')");
  verboseState.baseExists = isequal(hasVerbose, 1);
  if verboseState.baseExists
    verboseState.baseValue = logical(evalin('base', 'doaToolsVerbose'));
  end
catch
  verboseState.baseExists = false;
  verboseState.baseValue = false;
end

try
  verboseState.appDataExists = isappdata(0, 'doaToolsVerbose');
  if verboseState.appDataExists
    verboseState.appDataValue = logical(getappdata(0, 'doaToolsVerbose'));
  end
catch
  verboseState.appDataExists = false;
  verboseState.appDataValue = false;
end
end


function localAssignVerboseFlag(verbose)
%LOCALASSIGNVERBOSEFLAG Publish one inherited verbose flag to base/appdata.

try
  assignin('base', 'doaToolsVerbose', logical(verbose));
catch
end

try
  setappdata(0, 'doaToolsVerbose', logical(verbose));
catch
end
end


function localRestoreVerboseState(verboseState)
%LOCALRESTOREVERBOSESTATE Restore the prior inherited verbose carriers.

try
  if isstruct(verboseState) && isfield(verboseState, 'baseExists') && verboseState.baseExists
    assignin('base', 'doaToolsVerbose', logical(verboseState.baseValue));
  else
    evalin('base', 'if exist(''doaToolsVerbose'', ''var''), clear(''doaToolsVerbose''); end');
  end
catch
end

try
  if isstruct(verboseState) && isfield(verboseState, 'appDataExists') && verboseState.appDataExists
    setappdata(0, 'doaToolsVerbose', logical(verboseState.appDataValue));
  else
    if isappdata(0, 'doaToolsVerbose')
      rmappdata(0, 'doaToolsVerbose');
    end
  end
catch
end
end
