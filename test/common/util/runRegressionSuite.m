function result = runRegressionSuite(suiteName, scriptRelPathList, varargin)
%RUNREGRESSIONSUITE Run one grouped regression suite.
%
% Syntax:
%   result = runRegressionSuite(suiteName, scriptRelPathList)
%   result = runRegressionSuite(suiteName, scriptRelPathList, runOpt)
%   result = runRegressionSuite(suiteName, scriptRelPathList, Name, Value)
%
% Inputs:
%   suiteName         - suite tag shown in summary messages.
%   scriptRelPathList - regression script paths relative to project root.
%
% Optional inputs:
%   runOpt.projectRoot   - repository root. Default: auto resolve.
%   runOpt.stopOnFailure - stop on first failure. Default: false.
%   runOpt.verbose       - print per-script progress. Default: true.
%
% Output:
%   result            - suite result structure.
%
% Notes:
%   - Each regression script is executed in the base workspace so that the
%     common "clear; close all; clc;" pattern inside scripts does not wipe
%     this runner's bookkeeping variables.
%   - This helper only orchestrates grouped execution. Regression-specific
%     assertions stay inside each individual regression script.

if ~(ischar(suiteName) || (isstring(suiteName) && isscalar(suiteName)))
  error('runRegressionSuite:InvalidSuiteName', ...
    'suiteName must be a character vector or string scalar.');
end

scriptRelPathList = string(scriptRelPathList);
scriptRelPathList = scriptRelPathList(:);
if isempty(scriptRelPathList)
  error('runRegressionSuite:EmptyScriptList', ...
    'scriptRelPathList must contain at least one regression script.');
end

runOpt = localParseRunOpt(varargin{:});
projectRoot = localResolveProjectRoot(runOpt.projectRoot);
addpath(genpath(projectRoot));

numScript = numel(scriptRelPathList);
item = repmat(localBuildEmptyItem(), numScript, 1);
totalRunTimeSec = 0;

if runOpt.verbose
  fprintf('=== Running %s regression suite (%d scripts) ===\n', ...
    char(string(suiteName)), numScript);
end

for iScript = 1:numScript
  scriptRelPath = char(scriptRelPathList(iScript));
  scriptPath = fullfile(projectRoot, scriptRelPath);

  itemUse = localBuildEmptyItem();
  itemUse.scriptName = string(localResolveScriptName(scriptRelPath));
  itemUse.scriptRelPath = string(scriptRelPath);
  itemUse.scriptPath = string(scriptPath);

  if runOpt.verbose
    fprintf('[%d/%d] %s\n', iScript, numScript, char(itemUse.scriptName));
  end

  if exist(scriptPath, 'file') ~= 2
    itemUse.status = "FAIL";
    itemUse.errorId = "runRegressionSuite:ScriptNotFound";
    itemUse.errorMessage = sprintf('Regression script not found: %s', scriptRelPath);
    item(iScript) = itemUse;

    if runOpt.verbose
      fprintf('  FAIL (missing file)\n');
    end

    if runOpt.stopOnFailure
      localPrintSuiteSummary(suiteName, item(1:iScript));
      error('runRegressionSuite:ScriptNotFound', '%s', itemUse.errorMessage);
    end
    continue;
  end

  runTimer = tic;
  try
    evalin('base', localBuildRunCommand(scriptPath));
    itemUse.status = "PASS";
  catch ME
    itemUse.status = "FAIL";
    itemUse.errorId = string(localNormalizeErrorId(ME.identifier));
    itemUse.errorMessage = string(ME.message);
    itemUse.errorLocation = string(localBuildErrorLocation(ME));
  end

  itemUse.runTimeSec = toc(runTimer);
  totalRunTimeSec = totalRunTimeSec + itemUse.runTimeSec;
  item(iScript) = itemUse;

  if runOpt.verbose
    fprintf('  %s (%.3f s)\n', char(itemUse.status), itemUse.runTimeSec);
    if itemUse.status == "FAIL"
      if strlength(itemUse.errorId) > 0
        fprintf('    %s\n', char(itemUse.errorId));
      end
      fprintf('    %s\n', char(itemUse.errorMessage));
    end
  end

  if itemUse.status == "FAIL" && runOpt.stopOnFailure
    localPrintSuiteSummary(suiteName, item(1:iScript));
    error('runRegressionSuite:RegressionFailed', ...
      'Regression failed in %s: %s\n%s', ...
      char(string(suiteName)), char(itemUse.scriptName), char(itemUse.errorMessage));
  end
end

result = struct();
result.suiteName = string(suiteName);
result.projectRoot = string(projectRoot);
result.scriptRelPathList = scriptRelPathList;
result.item = item;
result.summaryTable = localBuildSummaryTable(item);
result.numScript = numScript;
result.numPassed = nnz(string({item.status}) == "PASS");
result.numFailed = nnz(string({item.status}) == "FAIL");
result.allPassed = (result.numFailed == 0);
result.totalRunTimeSec = totalRunTimeSec;

if runOpt.verbose
  localPrintSuiteSummary(suiteName, item);
  if ~result.allPassed
    localPrintFailureDetail(item);
  end
end

if ~result.allPassed && ~runOpt.stopOnFailure
  firstFailIdx = find(string({item.status}) == "FAIL", 1, 'first');
  firstFail = item(firstFailIdx);
  error('runRegressionSuite:RegressionFailed', ...
    'Regression suite %s finished with %d failure(s). First failure: %s\n%s', ...
    char(string(suiteName)), result.numFailed, char(firstFail.scriptName), char(firstFail.errorMessage));
end
end


function runOpt = localParseRunOpt(varargin)
%LOCALPARSERUNOPT Parse runner options from struct or Name/Value pairs.

runOpt = struct();
runOpt.projectRoot = "";
runOpt.stopOnFailure = false;
runOpt.verbose = true;

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


function scriptName = localResolveScriptName(scriptRelPath)
%LOCALRESOLVESCRIPTNAME Extract script base name without extension.

[~, scriptName, ~] = fileparts(scriptRelPath);
end


function runCmd = localBuildRunCommand(scriptPath)
%LOCALBUILDRUNCOMMAND Build one base-workspace run command.

scriptPath = strrep(scriptPath, '''', '''''');
runCmd = sprintf("run('%s');", scriptPath);
end


function summaryTable = localBuildSummaryTable(item)
%LOCALBUILDSUMMARYTABLE Convert item array to a compact table.

numItem = numel(item);
scriptName = strings(numItem, 1);
status = strings(numItem, 1);
runTimeSec = NaN(numItem, 1);
errorId = strings(numItem, 1);
errorLocation = strings(numItem, 1);

for iItem = 1:numItem
  scriptName(iItem) = item(iItem).scriptName;
  status(iItem) = item(iItem).status;
  runTimeSec(iItem) = item(iItem).runTimeSec;
  errorId(iItem) = item(iItem).errorId;
  errorLocation(iItem) = item(iItem).errorLocation;
end

summaryTable = table(scriptName, status, runTimeSec, errorId, errorLocation);
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
  fprintf('  %s', char(item(iFail).scriptName));
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
