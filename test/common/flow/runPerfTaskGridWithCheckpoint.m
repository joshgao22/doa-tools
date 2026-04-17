function runState = runPerfTaskGridWithCheckpoint(taskGrid, sharedData, taskRunner, opt)
%RUNPERFTASKGRIDWITHCHECKPOINT Run one outer task grid with resumable checkpoint files.
% This helper keeps checkpoint logic at the perf orchestration layer. Each
% outer task is evaluated independently and saved to one small temporary MAT
% file so interrupted runs can resume without changing the estimator path.
%
% One task file is written per outer task index. A lightweight status MAT
% file is also updated on the client side so the run directory always shows
% visible progress while the task grid is still running.

arguments
  taskGrid (:, 1) struct
  sharedData
  taskRunner (1, 1) function_handle
  opt (1, 1) struct = struct()
end

opt = localResolveOpt(opt);

runDir = localBuildRunDir(opt.outputRoot, opt.runName, opt.runKey);
taskDir = fullfile(runDir, "task");
manifestFile = fullfile(runDir, "manifest.mat");
statusFile = fullfile(runDir, "status.mat");

if ~opt.resume && isfolder(runDir)
  rmdir(runDir, 's');
end
if ~isfolder(taskDir)
  mkdir(taskDir);
end

manifestInfo = struct();
manifestInfo.runName = opt.runName;
manifestInfo.runKey = opt.runKey;
manifestInfo.meta = opt.meta;
manifestInfo.taskGrid = taskGrid;
manifestInfo.numTask = numel(taskGrid);
manifestInfo.createdAt = datetime('now');
localWriteManifest(manifestFile, manifestInfo);

numTask = numel(taskGrid);
doneMask = localLoadDoneMask(taskDir, numTask);
numDone = nnz(doneMask);
localPrintStartStatus(runDir, numTask, numDone);
localWriteStatus(statusFile, opt.runName, opt.runKey, numTask, doneMask, numDone);
if ~isempty(opt.progressFcn) && numDone > 0
  opt.progressFcn(numDone);
end

todoIdx = find(~doneMask);
if opt.useParfor && ~isempty(todoIdx)
  progressQueue = parallel.pool.DataQueue;
  afterEach(progressQueue, @localHandleTaskCompleted);

  parfor iTodo = 1:numel(todoIdx)
    iTask = todoIdx(iTodo);
    taskInfo = taskGrid(iTask);
    taskResult = taskRunner(taskInfo, sharedData);
    localSaveTaskResult(taskDir, iTask, taskResult);
    send(progressQueue, iTask);
  end
else
  for iTodo = 1:numel(todoIdx)
    iTask = todoIdx(iTodo);
    taskInfo = taskGrid(iTask);
    taskResult = taskRunner(taskInfo, sharedData);
    localSaveTaskResult(taskDir, iTask, taskResult);
    localHandleTaskCompleted(iTask);
  end
end

resultCell = localCollectTaskResults(taskDir, numTask);
doneMask = localLoadDoneMask(taskDir, numTask);
numDone = nnz(doneMask);
localWriteStatus(statusFile, opt.runName, opt.runKey, numTask, doneMask, numDone);

runState = struct();
runState.runDir = runDir;
runState.taskDir = taskDir;
runState.manifestFile = manifestFile;
runState.statusFile = statusFile;
runState.resultCell = resultCell;
runState.doneMask = doneMask;
runState.numTask = numTask;
runState.numDone = numDone;
runState.isComplete = all(doneMask);

  function localHandleTaskCompleted(taskIndex)
  %LOCALHANDLETASKCOMPLETED Update progress and persistent status on client side.

    if taskIndex >= 1 && taskIndex <= numTask
      doneMask(taskIndex) = true;
    end
    numDone = nnz(doneMask);
    localWriteStatus(statusFile, opt.runName, opt.runKey, numTask, doneMask, numDone);
    if ~isempty(opt.progressFcn)
      opt.progressFcn(1);
    end
  end
end


function opt = localResolveOpt(opt)
%LOCALRESOLVEOPT Fill defaults and validate one option struct.

optDefault = struct();
optDefault.runName = "runPerfTaskGridWithCheckpoint";
optDefault.runKey = "default";
optDefault.outputRoot = fullfile("out", "checkpoint_tmp");
optDefault.useParfor = true;
optDefault.resume = true;
optDefault.meta = struct();
optDefault.progressFcn = [];

optField = fieldnames(opt);
validField = fieldnames(optDefault);
extraField = setdiff(optField, validField);
if ~isempty(extraField)
  error('runPerfTaskGridWithCheckpoint:UnknownOption', ...
    'Unknown option field(s): %s', strjoin(extraField, ', '));
end

for iField = 1:numel(validField)
  fieldName = validField{iField};
  if ~isfield(opt, fieldName) || isempty(opt.(fieldName))
    opt.(fieldName) = optDefault.(fieldName);
  end
end

if ~(ischar(opt.runName) || isStringScalar(opt.runName))
  error('runPerfTaskGridWithCheckpoint:InvalidRunName', ...
    'opt.runName must be a character vector or string scalar.');
end
if ~(ischar(opt.runKey) || isStringScalar(opt.runKey))
  error('runPerfTaskGridWithCheckpoint:InvalidRunKey', ...
    'opt.runKey must be a character vector or string scalar.');
end
if ~(ischar(opt.outputRoot) || isStringScalar(opt.outputRoot))
  error('runPerfTaskGridWithCheckpoint:InvalidOutputRoot', ...
    'opt.outputRoot must be a character vector or string scalar.');
end
if ~islogical(opt.useParfor) || ~isscalar(opt.useParfor)
  error('runPerfTaskGridWithCheckpoint:InvalidUseParfor', ...
    'opt.useParfor must be a logical scalar.');
end
if ~islogical(opt.resume) || ~isscalar(opt.resume)
  error('runPerfTaskGridWithCheckpoint:InvalidResume', ...
    'opt.resume must be a logical scalar.');
end
if ~isstruct(opt.meta) || ~isscalar(opt.meta)
  error('runPerfTaskGridWithCheckpoint:InvalidMeta', ...
    'opt.meta must be a scalar struct.');
end
if ~(isempty(opt.progressFcn) || isa(opt.progressFcn, 'function_handle'))
  error('runPerfTaskGridWithCheckpoint:InvalidProgressFcn', ...
    'opt.progressFcn must be empty or a function handle.');
end

opt.runName = string(opt.runName);
opt.runKey = string(opt.runKey);
opt.outputRoot = string(opt.outputRoot);
end


function runDir = localBuildRunDir(outputRoot, runName, runKey)
%LOCALBUILDRUNDIR Build one stable temporary checkpoint directory.

runDir = fullfile(outputRoot, runName, runKey);
end


function localWriteManifest(manifestFile, manifestInfo)
%LOCALWRITEMANIFEST Create or validate one run manifest.

if ~isfile(manifestFile)
  save(manifestFile, 'manifestInfo');
  return;
end

manifestData = load(manifestFile, 'manifestInfo');
if ~isfield(manifestData, 'manifestInfo')
  error('runPerfTaskGridWithCheckpoint:InvalidManifest', ...
    'Manifest file is missing the expected manifestInfo field: %s', manifestFile);
end

manifestSaved = manifestData.manifestInfo;
if ~isequaln(manifestSaved.runName, manifestInfo.runName) || ...
    ~isequaln(manifestSaved.runKey, manifestInfo.runKey) || ...
    ~isequaln(manifestSaved.meta, manifestInfo.meta) || ...
    ~isequaln(manifestSaved.taskGrid, manifestInfo.taskGrid)
  error('runPerfTaskGridWithCheckpoint:ManifestMismatch', ...
    ['Checkpoint manifest mismatch detected for runName=%s, runKey=%s. ', ...
     'Delete the temporary checkpoint directory or change runKey before rerunning.'], ...
    manifestInfo.runName, manifestInfo.runKey);
end
end


function localWriteStatus(statusFile, runName, runKey, numTask, doneMask, numDone)
%LOCALWRITESTATUS Persist one lightweight progress snapshot.

statusInfo = struct();
statusInfo.runName = runName;
statusInfo.runKey = runKey;
statusInfo.numTask = numTask;
statusInfo.numDone = numDone;
statusInfo.doneMask = doneMask;
statusInfo.updatedAt = datetime('now');
save(statusFile, 'statusInfo');
end


function localPrintStartStatus(runDir, numTask, numDone)
%LOCALPRINTSTARTSTATUS Print one concise start banner for new or resumed runs.

if numDone > 0
  fprintf('[checkpoint] Resume run from %s (%d/%d done).\n', runDir, numDone, numTask);
else
  fprintf('[checkpoint] Start new run at %s (0/%d done).\n', runDir, numTask);
end
end


function doneMask = localLoadDoneMask(taskDir, numTask)
%LOCALLOADDONEMASK Return true for completed outer tasks.

doneMask = false(numTask, 1);
for iTask = 1:numTask
  doneMask(iTask) = isfile(localTaskFile(taskDir, iTask));
end
end


function localSaveTaskResult(taskDir, taskIndex, taskResult)
%LOCALSAVETASKRESULT Save one completed task through a temporary file.

taskFile = localTaskFile(taskDir, taskIndex);
taskTmpFile = taskFile + ".tmp";
save(taskTmpFile, 'taskResult', '-v7.3');
movefile(taskTmpFile, taskFile, 'f');
end


function resultCell = localCollectTaskResults(taskDir, numTask)
%LOCALCOLLECTTASKRESULTS Load all available task results back to the client.

resultCell = cell(numTask, 1);
for iTask = 1:numTask
  taskFile = localTaskFile(taskDir, iTask);
  if ~isfile(taskFile)
    continue;
  end
  taskData = load(taskFile, 'taskResult');
  if isfield(taskData, 'taskResult')
    resultCell{iTask} = taskData.taskResult;
  end
end
end


function taskFile = localTaskFile(taskDir, taskIndex)
%LOCALTASKFILE Build the per-task checkpoint filename.

taskFile = fullfile(taskDir, sprintf('task_%06d.mat', taskIndex));
end
