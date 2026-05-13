function result = runRegressionQuick(varargin)
%RUNREGRESSIONQUICK Run the fast invariant and branch regression suite.
caseRelPathList = [ ...
  "test/regression/invariant/regressionRefStateInvariant.m"; ...
  "test/regression/invariant/regressionUserStateFromLatlon.m"; ...
  "test/regression/invariant/regressionLocalDoaMapping.m"; ...
  "test/regression/invariant/regressionSfModelBuild.m"; ...
  "test/regression/invariant/regressionSfObjectiveShape.m"; ...
  "test/regression/invariant/regressionSfStaticRefOnlyEquivalence.m"; ...
  "test/regression/invariant/regressionSfStaticReferenceDopplerInvariant.m"; ...
  "test/regression/invariant/regressionSfStaticTruthReplayMultiSat.m"; ...
  "test/regression/invariant/regressionSfStaticFdInitZeroWeightInvariant.m"; ...
  "test/regression/invariant/regressionMsSfStaticObjectiveDecomposition.m"; ...
  "test/regression/invariant/regressionMsSfStaticW0EqualsSsStatic.m"; ...
  "test/regression/invariant/regressionMsSfStaticWeightInjection.m"; ...
  "test/regression/invariant/regressionMfModelBuild.m"; ...
  "test/regression/invariant/regressionMfInitFdLine.m"; ...
  "test/regression/invariant/regressionSnapshotTruthReplayMf.m"; ...
  "test/regression/invariant/regressionOtherSatOnlyMf.m"; ...
  "test/regression/invariant/regressionMfCpSupportCollapsePenalty.m"; ...
  "test/regression/invariant/regressionMfCpIpTimeAxisInvariant.m"; ...
  "test/regression/invariant/regressionCrbKnownUnknownConsistency.m"; ...
  "test/regression/branch/regressionSfStaticFixedDoaFdRefSolve.m"; ...
  "test/regression/branch/regressionSfStaticDoaAnchorFallback.m"; ...
  "test/regression/branch/regressionSfStaticJointCouplingCaseLocked.m"; ...
  "test/regression/branch/regressionSfStaticWeightSweepCouplingCaseLocked.m"; ...
  "test/regression/branch/regressionMfSubsetSelectNoTruthLeak.m"; ...
  "test/regression/branch/regressionMfFastSubsetEscalation.m"; ...
  "test/regression/branch/regressionMfDoaBasinEntrySplit.m"; ...
  "test/regression/branch/regressionMfUnknownFinalSelectionRules.m"; ...
  "test/regression/branch/regressionMfUnknownFixedDoaWarmAnchor.m"; ...
  "test/regression/branch/regressionMfUnknownWarmStartSet.m"; ...
  "test/regression/branch/regressionMfUnknownReleaseFromCpK.m"; ...
  "test/regression/branch/regressionMfUnknownBestStartSelection.m"; ...
  "test/regression/branch/regressionMfBranchKnownUnknown.m"; ...
  "test/regression/branch/regressionMfWarmAnchorParforGate.m"; ...
  "test/regression/branch/regressionMfUnknownNoisyWrongToothGuard.m" ...
  ];
runOpt = localParseQuickRunOpt(varargin{:});
if localShouldUseParfor(runOpt)
  result = localRunQuickSuiteParfor(caseRelPathList, runOpt);
else
  result = runRegressionSuite("quick", caseRelPathList, localBuildSerialRunOpt(runOpt));
end
end

function runOpt = localParseQuickRunOpt(varargin)
%LOCALPARSEQUICKRUNOPT Parse quick-runner options.

runOpt = struct();
runOpt.projectRoot = "";
runOpt.stopOnFailure = false;
runOpt.verbose = false;
runOpt.useParfor = [];

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
  if isfield(optIn, 'stopOnError') && ~isempty(optIn.stopOnError)
    runOpt.stopOnFailure = logical(optIn.stopOnError);
  end
  if isfield(optIn, 'verbose') && ~isempty(optIn.verbose)
    runOpt.verbose = logical(optIn.verbose);
  end
  if isfield(optIn, 'useParfor') && ~isempty(optIn.useParfor)
    runOpt.useParfor = localParseUseParforValue(optIn.useParfor);
  end
  return;
end

if mod(nargin, 2) ~= 0
  error('runRegressionQuick:InvalidOptionInput', ...
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
    case {'stoponerror', 'stoponfail'}
      runOpt.stopOnFailure = logical(value);
    case 'verbose'
      runOpt.verbose = logical(value);
    case {'useparfor', 'parfor', 'parallel'}
      runOpt.useParfor = localParseUseParforValue(value);
    otherwise
      error('runRegressionQuick:UnknownOption', ...
        'Unknown option name: %s', char(string(varargin{iArg})));
  end
end
end

function useParfor = localParseUseParforValue(value)
%LOCALPARSEUSEPARFORVALUE Parse a logical-or-auto parfor option.

if ischar(value) || (isstring(value) && isscalar(value))
  valueText = lower(strtrim(string(value)));
  if valueText == "auto"
    useParfor = [];
    return;
  end
  if ismember(valueText, ["true", "on", "yes", "1"])
    useParfor = true;
    return;
  end
  if ismember(valueText, ["false", "off", "no", "0"])
    useParfor = false;
    return;
  end
  error('runRegressionQuick:InvalidUseParfor', ...
    'useParfor must be logical or "auto".');
end
useParfor = logical(value);
end

function serialOpt = localBuildSerialRunOpt(runOpt)
%LOCALBUILDSERIALRUNOPT Forward only runRegressionSuite-supported options.

serialOpt = struct();
serialOpt.projectRoot = runOpt.projectRoot;
serialOpt.stopOnFailure = runOpt.stopOnFailure;
serialOpt.verbose = runOpt.verbose;
end

function useParfor = localShouldUseParfor(runOpt)
%LOCALSHOULDUSEPARFOR Resolve quick-runner parallel execution.

autoMode = isempty(runOpt.useParfor);
if autoMode
  useParfor = localParallelToolboxAvailable();
else
  useParfor = logical(runOpt.useParfor);
end
if useParfor && autoMode && runOpt.verbose
  fprintf(['runRegressionQuick: verbose=true is clearer in serial mode; ', ...
    'parfor is disabled for this auto run.\n']);
  useParfor = false;
end
if useParfor && runOpt.stopOnFailure
  fprintf(['runRegressionQuick: stopOnFailure=true requires serial execution; ', ...
    'parfor is disabled for this run.\n']);
  useParfor = false;
end
if useParfor && ~localParallelToolboxAvailable()
  fprintf('runRegressionQuick: parallel toolbox unavailable; falling back to serial execution.\n');
  useParfor = false;
end
end

function tf = localParallelToolboxAvailable()
%LOCALPARALLELTOOLBOXAVAILABLE Check whether parfor can be used.

tf = false;
try
  tf = exist('parpool', 'file') == 2 && exist('gcp', 'file') == 2 && ...
    license('test', 'Distrib_Computing_Toolbox');
catch
  tf = false;
end
end

function result = localRunQuickSuiteParfor(caseRelPathList, runOpt)
%LOCALRUNQUICKSUITEPARFOR Run quick regression cases in parallel.

suiteName = "quick";
caseRelPathList = string(caseRelPathList(:));
projectRoot = localResolveProjectRoot(runOpt.projectRoot);
addpath(genpath(projectRoot));

numCase = numel(caseRelPathList);
item = repmat(localBuildEmptyItem(), numCase, 1);

fprintf('=== Running quick regression suite (%d cases, parfor) ===\n', numCase);
wallTimer = tic;
parfor iCase = 1:numCase
  item(iCase) = localRunOneQuickCase(projectRoot, caseRelPathList(iCase), runOpt.verbose);
end
wallTimeSec = toc(wallTimer);

result = struct();
result.suiteName = suiteName;
result.projectRoot = string(projectRoot);
result.caseRelPathList = caseRelPathList;
result.scriptRelPathList = caseRelPathList;
result.item = item;
result.summaryTable = localBuildSummaryTable(item);
result.numCase = numCase;
result.numScript = numCase;
result.numPassed = nnz(string({item.status}) == "PASS");
result.numFailed = nnz(string({item.status}) == "FAIL");
result.allPassed = (result.numFailed == 0);
result.totalRunTimeSec = sum([item.runTimeSec], 'omitnan');
result.wallTimeSec = wallTimeSec;
result.useParfor = true;

localPrintSuiteSummary(suiteName, item, wallTimeSec);
if ~result.allPassed
  localPrintFailureDetail(item);
  firstFailIdx = find(string({item.status}) == "FAIL", 1, 'first');
  firstFail = item(firstFailIdx);
  error('runRegressionQuick:RegressionFailed', ...
    'Regression suite quick finished with %d failure(s). First failure: %s\n%s', ...
    result.numFailed, char(firstFail.caseName), char(firstFail.errorMessage));
end
end

function itemUse = localRunOneQuickCase(projectRoot, caseRelPath, verbose)
%LOCALRUNONEQUICKCASE Run one function-style regression case on a worker.

caseRelPath = char(string(caseRelPath));
casePath = fullfile(projectRoot, caseRelPath);
caseName = localResolveCaseName(caseRelPath);

itemUse = localBuildEmptyItem();
itemUse.caseName = string(caseName);
itemUse.caseRelPath = string(caseRelPath);
itemUse.casePath = string(casePath);
itemUse.scriptName = itemUse.caseName;
itemUse.scriptRelPath = itemUse.caseRelPath;
itemUse.scriptPath = itemUse.casePath;

try
  addpath(genpath(projectRoot));
catch
end

if exist(casePath, 'file') ~= 2
  itemUse.status = "FAIL";
  itemUse.errorId = "runRegressionQuick:CaseNotFound";
  itemUse.errorMessage = sprintf('Regression case not found: %s', caseRelPath);
  itemUse.runTimeSec = 0;
  localPrintQuickCaseStatus(itemUse);
  return;
end

runTimer = tic;
consoleText = "";
try
  localAssignVerboseFlag(verbose);
  consoleText = string(evalc('feval(caseName, ''verbose'', verbose);'));
  itemUse.status = "PASS";
catch ME
  itemUse.status = "FAIL";
  itemUse.errorId = string(localNormalizeErrorId(ME.identifier));
  itemUse.errorMessage = string(ME.message);
  itemUse.errorLocation = string(localBuildErrorLocation(ME));
end
itemUse.runTimeSec = toc(runTimer);
itemUse.consoleText = consoleText;
localPrintQuickCaseStatus(itemUse);
end

function localPrintQuickCaseStatus(itemUse)
%LOCALPRINTQUICKCASESTATUS Print one completed parfor case line.

caseName = char(itemUse.caseName);
statusText = char(itemUse.status);
if itemUse.status == "PASS"
  fprintf('[quick/parfor] %s %s (%.3f s)\n', statusText, caseName, itemUse.runTimeSec);
  return;
end

if strlength(itemUse.errorId) > 0
  failText = sprintf('%s: %s', char(itemUse.errorId), char(itemUse.errorMessage));
elseif strlength(itemUse.errorMessage) > 0
  failText = char(itemUse.errorMessage);
else
  failText = 'Unknown failure';
end
fprintf('[quick/parfor] %s %s (%.3f s) - %s\n', ...
  statusText, caseName, itemUse.runTimeSec, failText);
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
item.consoleText = "";
end

function projectRoot = localResolveProjectRoot(projectRootIn)
%LOCALRESOLVEPROJECTROOT Resolve project root from input or file path.

if strlength(string(projectRootIn)) > 0
  projectRoot = char(string(projectRootIn));
  return;
end
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(thisDir));
end

function caseName = localResolveCaseName(caseRelPath)
%LOCALRESOLVECASEName Extract case base name without extension.

[~, caseName, ~] = fileparts(caseRelPath);
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

function localPrintSuiteSummary(suiteName, item, wallTimeSec)
%LOCALPRINTSUITESUMMARY Print one compact suite summary.

numPassed = nnz(string({item.status}) == "PASS");
numFailed = nnz(string({item.status}) == "FAIL");
numTotal = numel(item);
totalRunTimeSec = sum([item.runTimeSec], 'omitnan');
fprintf('=== %s regression summary: %d/%d passed, %d failed, total %.3f s, wall %.3f s ===\n', ...
  char(string(suiteName)), numPassed, numTotal, numFailed, totalRunTimeSec, wallTimeSec);
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
  errorId = 'runRegressionQuick:UnknownError';
end
end

function localAssignVerboseFlag(verbose)
%LOCALASSIGNVERBOSEFLAG Publish one inherited verbose flag to worker state.

try
  assignin('base', 'doaToolsVerbose', logical(verbose));
catch
end
try
  setappdata(0, 'doaToolsVerbose', logical(verbose));
catch
end
end
