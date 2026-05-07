function runState = runPerfTaskGridWithCheckpoint(taskGrid, sharedData, taskRunner, opt)
%RUNPERFTASKGRIDWITHCHECKPOINT Run one outer task grid with resumable checkpoint files.
% This helper keeps checkpoint logic at the perf orchestration layer. Each
% outer task is evaluated independently and saved to one small temporary MAT
% file so interrupted runs can resume without changing the estimator path.
%
% One task file is written per outer task index. On success, callers can set
% cleanupOnSuccess=true for immediate cleanup or call cleanupPerfTaskGridCheckpoint
% after assembling and saving their final lightweight result struct.

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
todoIdx = find(~doneMask);
localPrintStartStatus(runDir, numDone, numTask);

if opt.useParfor && ~isempty(todoIdx)
  localRunParallelCheckpointTasks(taskGrid, sharedData, taskRunner, taskDir, todoIdx, opt.progressFcn);
else
  for iTodo = 1:numel(todoIdx)
    iTask = todoIdx(iTodo);
    taskInfo = taskGrid(iTask);
    taskResult = taskRunner(taskInfo, sharedData);
    localSaveTaskResult(taskDir, iTask, taskResult);
    if ~isempty(opt.progressFcn)
      opt.progressFcn(1);
    end
  end
end

resultCell = localCollectTaskResults(taskDir, numTask);
doneMask = localLoadDoneMask(taskDir, numTask);

runState = struct();
runState.runName = opt.runName;
runState.runKey = opt.runKey;
runState.outputRoot = opt.outputRoot;
runState.runDir = runDir;
runState.taskDir = taskDir;
runState.manifestFile = manifestFile;
runState.resultCell = resultCell;
runState.doneMask = doneMask;
runState.numTask = numTask;
runState.numDone = nnz(doneMask);
runState.isComplete = all(doneMask);
runState.cleanupOnSuccess = opt.cleanupOnSuccess;
runState.cleanedOnSuccess = false;
if opt.cleanupOnSuccess && runState.isComplete
  runState.cleanupReport = cleanupPerfTaskGridCheckpoint(runState, opt.cleanupOpt);
  runState.cleanedOnSuccess = true;
else
  runState.cleanupReport = struct([]);
end
end


function opt = localResolveOpt(opt)
%LOCALRESOLVEOPT Resolve defaults and validate checkpoint runner options.

allowedField = { ...
  'runName', ...
  'runKey', ...
  'outputRoot', ...
  'useParfor', ...
  'resume', ...
  'meta', ...
  'progressFcn', ...
  'cleanupOnSuccess', ...
  'cleanupOpt'};

optField = fieldnames(opt);
extraField = setdiff(optField, allowedField);
if ~isempty(extraField)
  error('runPerfTaskGridWithCheckpoint:UnknownOption', ...
    'Unknown checkpoint option field: %s', strjoin(extraField, ', '));
end

if ~isfield(opt, 'runName') || isempty(opt.runName)
  error('runPerfTaskGridWithCheckpoint:MissingRunName', ...
    'checkpointOpt.runName is required.');
end
if ~isfield(opt, 'runKey') || isempty(opt.runKey)
  opt.runKey = "default";
end
if ~isfield(opt, 'outputRoot') || isempty(opt.outputRoot)
  opt.outputRoot = fullfile("out", "checkpoint_tmp");
end
if ~isfield(opt, 'useParfor') || isempty(opt.useParfor)
  opt.useParfor = true;
end
if ~isfield(opt, 'resume') || isempty(opt.resume)
  opt.resume = true;
end
if ~isfield(opt, 'meta') || isempty(opt.meta)
  opt.meta = struct();
end
if ~isfield(opt, 'progressFcn')
  opt.progressFcn = [];
end
if ~isfield(opt, 'cleanupOnSuccess') || isempty(opt.cleanupOnSuccess)
  opt.cleanupOnSuccess = false;
end
if ~isfield(opt, 'cleanupOpt') || isempty(opt.cleanupOpt)
  opt.cleanupOpt = struct();
end

opt.runName = string(opt.runName);
opt.runKey = string(opt.runKey);
opt.outputRoot = string(opt.outputRoot);
opt.useParfor = logical(opt.useParfor);
opt.resume = logical(opt.resume);
opt.cleanupOnSuccess = logical(opt.cleanupOnSuccess);
if ~isa(opt.meta, 'struct')
  error('runPerfTaskGridWithCheckpoint:InvalidMeta', ...
    'checkpointOpt.meta must be a struct.');
end
if ~isempty(opt.progressFcn) && ~isa(opt.progressFcn, 'function_handle')
  error('runPerfTaskGridWithCheckpoint:InvalidProgressFcn', ...
    'checkpointOpt.progressFcn must be empty or a function handle.');
end
if ~isstruct(opt.cleanupOpt)
  error('runPerfTaskGridWithCheckpoint:InvalidCleanupOpt', ...
    'checkpointOpt.cleanupOpt must be a struct.');
end
end


function localPrintStartStatus(runDir, numDone, numTask)
%LOCALPRINTSTARTSTATUS Print one concise checkpoint start banner.

if numDone > 0
  fprintf('[checkpoint] Resume run from %s (%d/%d done).\n', runDir, numDone, numTask);
else
  fprintf('[checkpoint] Start new run at %s (0/%d done).\n', runDir, numTask);
end
end


function localRunParallelCheckpointTasks(taskGrid, sharedData, taskRunner, taskDir, todoIdx, progressFcn)
%LOCALRUNPARALLELCHECKPOINTTASKS Run checkpoint tasks without returning large payloads.
%
% Each worker writes the full task result to its checkpoint MAT file and only
% returns a tiny status struct to the client. This avoids the parfor gather
% path deserializing large replay task payloads after all workers finish.

poolObj = gcp('nocreate');
if isempty(poolObj)
  poolObj = parpool;
end

sharedDataConst = parallel.pool.Constant(sharedData);
numTodo = numel(todoIdx);
futureList = parallel.FevalFuture.empty(numTodo, 0);
for iTodo = 1:numTodo
  iTask = todoIdx(iTodo);
  taskInfo = taskGrid(iTask);
  futureList(iTodo, 1) = parfeval(poolObj, @localRunAndSaveCheckpointTask, 1, ...
    taskDir, iTask, taskInfo, sharedDataConst, taskRunner);
end
futureCleanup = onCleanup(@() localCancelUnfinishedFutures(futureList)); %#ok<NASGU>

for iDone = 1:numTodo
  [~, taskStatus] = fetchNext(futureList);
  localValidateCheckpointTaskStatus(taskStatus);
  if ~isempty(progressFcn)
    progressFcn(1);
  end
end
end


function taskStatus = localRunAndSaveCheckpointTask(taskDir, taskIndex, taskInfo, sharedDataConst, taskRunner)
%LOCALRUNANDSAVECHECKPOINTTASK Run one task on a worker and save its result.

sharedData = sharedDataConst.Value;
taskResult = taskRunner(taskInfo, sharedData);
localSaveTaskResult(taskDir, taskIndex, taskResult);
taskStatus = struct();
taskStatus.taskIndex = taskIndex;
taskStatus.taskFile = string(localTaskFile(taskDir, taskIndex));
end


function localValidateCheckpointTaskStatus(taskStatus)
%LOCALVALIDATECHECKPOINTTASKSTATUS Validate the tiny status returned by workers.

if ~(isstruct(taskStatus) && isfield(taskStatus, 'taskIndex') && isfield(taskStatus, 'taskFile'))
  error('runPerfTaskGridWithCheckpoint:InvalidTaskStatus', ...
    'Parallel checkpoint task returned an invalid status payload.');
end
if ~(isfinite(double(taskStatus.taskIndex)) && isfile(taskStatus.taskFile))
  error('runPerfTaskGridWithCheckpoint:MissingTaskFile', ...
    'Parallel checkpoint task %d did not create the expected task file: %s', ...
    double(taskStatus.taskIndex), char(taskStatus.taskFile));
end
end


function localCancelUnfinishedFutures(futureList)
%LOCALCANCELUNFINISHEDFUTURES Cancel pending futures after an error.

for iFuture = 1:numel(futureList)
  if ~strcmp(futureList(iFuture).State, 'finished')
    cancel(futureList(iFuture));
  end
end
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
