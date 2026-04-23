% Development script for validating saveExpSnapshot / loadExpSnapshot.
% The script exercises default full-save behavior, include / exclude modes,
% maxVarBytes filtering, metadata-only load, selective restore, caller/base
% workspace restore, and legacy snapshot compatibility.
clear(); close all; clc;

repoRoot = localRepoRoot();
sessionDir = fullfile(repoRoot, 'test', 'data', 'cache', ...
  ['devSaveLoadExpSnapshot_' datestr(now, 'yyyymmdd-HHMMSS')]);
mkdir(sessionDir);

results = struct('name', {}, 'passed', {}, 'detail', {});

results(end + 1) = localRunTest('saveAll_customDir', @() localTestSaveAll(sessionDir));
results(end + 1) = localRunTest('defaultOutputDir_defaultPrefix', @() localTestDefaultOutputDir());
results(end + 1) = localRunTest('includeVars_only', @() localTestIncludeVars(sessionDir));
results(end + 1) = localRunTest('includeVars_overrideMaxVarBytes', @() localTestIncludeOverridesMaxBytes(sessionDir));
results(end + 1) = localRunTest('includeVars_missingVar', @() localTestIncludeMissingVar(sessionDir));
results(end + 1) = localRunTest('excludeVars_plusMaxVarBytes', @() localTestExcludePlusMaxBytes(sessionDir));
results(end + 1) = localRunTest('includeExclude_conflict', @() localTestIncludeExcludeConflict(sessionDir));
results(end + 1) = localRunTest('load_varNames_strictModes', @() localTestLoadVarNames(sessionDir));
results(end + 1) = localRunTest('load_metaOnly', @() localTestLoadMetaOnly(sessionDir));
results(end + 1) = localRunTest('load_toCaller', @() localTestLoadToCaller(sessionDir));
results(end + 1) = localRunTest('load_toBase', @() localTestLoadToBase(sessionDir));
results(end + 1) = localRunTest('legacySnapshot_compatibility', @() localTestLegacySnapshot(sessionDir));

resultTable = struct2table(results);
disp(resultTable);

numPass = nnz([results.passed]);
numFail = numel(results) - numPass;
fprintf('\nSummary: %d / %d tests passed.\n', numPass, numel(results));
if numFail > 0
  fprintf('Artifacts preserved under: %s\n', sessionDir);
  error('devSaveLoadExpSnapshot:TestFailed', ...
    '%d tests failed. Inspect resultTable for details.', numFail);
end

cleanupRunArtifacts(sessionDir, struct( ...
  'requiredPrefix', 'devSaveLoadExpSnapshot_', ...
  'verbose', false));
fprintf('Temporary cache artifacts cleaned: %s\n', sessionDir);


function result = localRunTest(name, testFcn)
result = struct('name', string(name), 'passed', false, 'detail', "");
try
  detail = testFcn();
  result.passed = true;
  result.detail = string(detail);
catch err
  result.passed = false;
  result.detail = string(sprintf('%s: %s', err.identifier, err.message));
end
end


function detail = localTestSaveAll(sessionDir)
outDir = fullfile(sessionDir, 'saveAll');
mkdir(outDir);

alpha = 42;
beta = magic(3);
gammaStruct = struct('tag', 'all', 'vec', [1 3 5]); %#ok<NASGU>
noteText = "saveAllCase"; %#ok<NASGU>
resultFile = saveExpSnapshot('devSaveLoad_all', struct( ...
  'outputDir', outDir, ...
  'verbose', false));

[data, meta, loadedVar, inventory] = loadExpSnapshot(resultFile, 'none', struct());

localAssert(isfield(data, 'alpha') && isequal(data.alpha, alpha), ...
  'devSaveLoadExpSnapshot:SaveAllAlpha', 'alpha was not saved correctly.');
localAssert(isfield(data, 'beta') && isequal(data.beta, beta), ...
  'devSaveLoadExpSnapshot:SaveAllBeta', 'beta was not saved correctly.');
localAssert(isfield(data, 'gammaStruct'), ...
  'devSaveLoadExpSnapshot:SaveAllGamma', 'gammaStruct was not saved.');
localAssert(strcmp(meta.selectionMode, 'all'), ...
  'devSaveLoadExpSnapshot:SaveAllMode', 'selectionMode should be all.');
localAssert(any(strcmp(loadedVar, 'alpha')) && any(strcmp(loadedVar, 'beta')), ...
  'devSaveLoadExpSnapshot:SaveAllLoadedVar', 'loadedVar is missing expected names.');
localAssert(any(strcmp({inventory.name}, 'alpha')) && any([inventory.isSaved]), ...
  'devSaveLoadExpSnapshot:SaveAllInventory', 'inventory did not record saved variables.');

detail = sprintf('saved %d vars, selectionMode=%s', numel(fieldnames(data)), meta.selectionMode);
end


function detail = localTestDefaultOutputDir()
probeValue = 12345; %#ok<NASGU>
resultFile = saveExpSnapshot('', struct( ...
  'includeVars', {{'probeValue'}}, ...
  'verbose', false));
cleanupObj = onCleanup(@() cleanupRunArtifacts(resultFile, struct( ...
  'requiredPrefix', 'saveExpSnapshot_', ...
  'verbose', false))); %#ok<NASGU>
[data, meta, ~, ~] = loadExpSnapshot(resultFile, 'none', struct());

expectedDir = localCanonicalPath(fullfile(localRepoRoot(), 'test', 'data', 'cache'));
[saveDir, saveName, saveExt] = fileparts(resultFile);
actualDir = localCanonicalPath(saveDir);
metaDir = localCanonicalPath(meta.outputDir);
localAssert(strcmpi(actualDir, expectedDir), ...
  'devSaveLoadExpSnapshot:DefaultDir', 'Default outputDir is incorrect.');
localAssert(strcmpi(metaDir, expectedDir), ...
  'devSaveLoadExpSnapshot:MetaOutputDir', 'meta.outputDir is inconsistent with default outputDir.');
localAssert(isfield(data, 'probeValue') && data.probeValue == probeValue, ...
  'devSaveLoadExpSnapshot:DefaultSaveData', 'Default save path did not preserve probeValue.');
localAssert(~isempty(regexp(saveName, '^[A-Za-z0-9_]+_[0-9]{8}-[0-9]{6}$', 'once')), ...
  'devSaveLoadExpSnapshot:DefaultNamePattern', 'Default snapshot name should match prefix_timestamp.');
localAssert(strcmp(saveExt, '.mat'), ...
  'devSaveLoadExpSnapshot:DefaultExt', 'Snapshot extension should be .mat.');
localAssert(strcmp(meta.selectionMode, 'include'), ...
  'devSaveLoadExpSnapshot:DefaultMode', 'selectionMode should be include.');

nameParts = split(string(saveName), '_');
resolvedPrefix = strjoin(nameParts(1:end-1), '_');
detail = sprintf('default dir=%s, prefix=%s', actualDir, resolvedPrefix);
end


function detail = localTestIncludeVars(sessionDir)
outDir = fullfile(sessionDir, 'includeOnly');
mkdir(outDir);

keepA = 11; %#ok<NASGU>
keepB = struct('v', [2 4 6]); %#ok<NASGU>
dropC = pi; %#ok<NASGU>
resultFile = saveExpSnapshot('devSaveLoad_include', struct( ...
  'includeVars', {{'keepA', 'keepB'}}, ...
  'outputDir', outDir, ...
  'verbose', false));

[data, meta, ~, inventory] = loadExpSnapshot(resultFile, 'none', struct());
nameList = sort(fieldnames(data));
localAssert(isequal(nameList, sort({'keepA'; 'keepB'})), ...
  'devSaveLoadExpSnapshot:IncludeFields', 'includeVars did not limit saved fields correctly.');
localAssert(strcmp(meta.selectionMode, 'include'), ...
  'devSaveLoadExpSnapshot:IncludeMode', 'selectionMode should be include.');
idxDrop = find(strcmp({inventory.name}, 'dropC'), 1, 'first');
localAssert(~isempty(idxDrop) && strcmp(inventory(idxDrop).skipReason, 'notIncluded'), ...
  'devSaveLoadExpSnapshot:IncludeInventory', 'dropC should be marked notIncluded.');

detail = sprintf('saved fields=%s', strjoin(nameList, ','));
end


function detail = localTestIncludeOverridesMaxBytes(sessionDir)
outDir = fullfile(sessionDir, 'includeMaxOverride');
mkdir(outDir);

largeKeep = rand(1, 4096); %#ok<NASGU>
smallDrop = 7; %#ok<NASGU>
resultFile = saveExpSnapshot('devSaveLoad_includeLarge', struct( ...
  'includeVars', {{'largeKeep'}}, ...
  'maxVarBytes', 16, ...
  'outputDir', outDir, ...
  'verbose', false));

[data, ~, ~, inventory] = loadExpSnapshot(resultFile, 'none', struct());
localAssert(isfield(data, 'largeKeep'), ...
  'devSaveLoadExpSnapshot:IncludeLargeField', 'largeKeep should be preserved by includeVars.');
idxKeep = find(strcmp({inventory.name}, 'largeKeep'), 1, 'first');
localAssert(~isempty(idxKeep) && inventory(idxKeep).isSaved, ...
  'devSaveLoadExpSnapshot:IncludeLargeInventory', 'largeKeep should be marked saved.');

detail = sprintf('largeKeep bytes=%d preserved with maxVarBytes=%d', ...
  inventory(idxKeep).bytes, 16);
end


function detail = localTestIncludeMissingVar(sessionDir)
outDir = fullfile(sessionDir, 'includeMissing');
mkdir(outDir);

a = 1; %#ok<NASGU>
err = localExpectError(@() saveExpSnapshot('devSaveLoad_missing', struct( ...
  'includeVars', {{'a', 'missingVar'}}, ...
  'outputDir', outDir, ...
  'verbose', false)));
localAssert(strcmp(err.identifier, 'saveExpSnapshot:MissingIncludeVar'), ...
  'devSaveLoadExpSnapshot:MissingIncludeId', 'Unexpected error for missing includeVars variable.');

detail = err.identifier;
end


function detail = localTestExcludePlusMaxBytes(sessionDir)
outDir = fullfile(sessionDir, 'excludeMax');
mkdir(outDir);

keepSmall = 21; %#ok<NASGU>
excludeMe = magic(4); %#ok<NASGU>
largeAutoSkip = rand(1, 4096); %#ok<NASGU>
resultFile = saveExpSnapshot('devSaveLoad_exclude', struct( ...
  'excludeVars', {{'excludeMe'}}, ...
  'maxVarBytes', 128, ...
  'outputDir', outDir, ...
  'verbose', false));

[data, meta, ~, inventory] = loadExpSnapshot(resultFile, 'none', struct());
localAssert(isfield(data, 'keepSmall'), ...
  'devSaveLoadExpSnapshot:ExcludeKeep', 'keepSmall should be saved.');
localAssert(~isfield(data, 'excludeMe'), ...
  'devSaveLoadExpSnapshot:ExcludeDrop', 'excludeMe should be excluded.');
localAssert(~isfield(data, 'largeAutoSkip'), ...
  'devSaveLoadExpSnapshot:ExcludeLarge', 'largeAutoSkip should be skipped by maxVarBytes.');
idxExclude = find(strcmp({inventory.name}, 'excludeMe'), 1, 'first');
idxLarge = find(strcmp({inventory.name}, 'largeAutoSkip'), 1, 'first');
localAssert(strcmp(meta.selectionMode, 'exclude'), ...
  'devSaveLoadExpSnapshot:ExcludeMode', 'selectionMode should be exclude.');
localAssert(strcmp(inventory(idxExclude).skipReason, 'excluded'), ...
  'devSaveLoadExpSnapshot:ExcludeReason', 'excludeMe should be marked excluded.');
localAssert(strcmp(inventory(idxLarge).skipReason, 'tooLarge'), ...
  'devSaveLoadExpSnapshot:TooLargeReason', 'largeAutoSkip should be marked tooLarge.');

detail = sprintf('excluded=%s, tooLarge=%s', inventory(idxExclude).name, inventory(idxLarge).name);
end


function detail = localTestIncludeExcludeConflict(sessionDir)
outDir = fullfile(sessionDir, 'conflict');
mkdir(outDir);

a = 1; %#ok<NASGU>
err = localExpectError(@() saveExpSnapshot('devSaveLoad_conflict', struct( ...
  'includeVars', {{'a'}}, ...
  'excludeVars', {{'a'}}, ...
  'outputDir', outDir, ...
  'verbose', false)));
localAssert(strcmp(err.identifier, 'saveExpSnapshot:ConflictingSelection'), ...
  'devSaveLoadExpSnapshot:ConflictId', 'Unexpected error for include/exclude conflict.');

detail = err.identifier;
end


function detail = localTestLoadVarNames(sessionDir)
outDir = fullfile(sessionDir, 'loadVarNames');
mkdir(outDir);
resultFile = localCreateRestoreSnapshot(outDir, 'devSaveLoad_varNames');

[dataSubset, ~, loadedVarSubset, ~] = loadExpSnapshot(resultFile, 'none', struct( ...
  'varNames', {{'restoredA'}}, ...
  'strict', false));
localAssert(numel(loadedVarSubset) == 1 && strcmp(loadedVarSubset{1}, 'restoredA'), ...
  'devSaveLoadExpSnapshot:VarNamesSubset', 'strict=false subset load returned wrong variables.');
localAssert(isfield(dataSubset, 'restoredA') && isfield(dataSubset, 'restoredB'), ...
  'devSaveLoadExpSnapshot:VarNamesData', 'data output should still contain full saved data.');

err = localExpectError(@() loadExpSnapshot(resultFile, 'none', struct( ...
  'varNames', {{'restoredA', 'missingVar'}}, ...
  'strict', true)));
localAssert(strcmp(err.identifier, 'loadExpSnapshot:MissingRequestedVar'), ...
  'devSaveLoadExpSnapshot:VarNamesStrictId', 'Unexpected error for strict missing varNames.');

detail = sprintf('subset loaded=%s; strict missing id=%s', ...
  strjoin(loadedVarSubset, ','), err.identifier);
end


function detail = localTestLoadMetaOnly(sessionDir)
outDir = fullfile(sessionDir, 'metaOnly');
mkdir(outDir);
resultFile = localCreateRestoreSnapshot(outDir, 'devSaveLoad_metaOnly');

[data, meta, loadedVar, inventory] = loadExpSnapshot(resultFile, 'none', struct('metaOnly', true));
localAssert(isstruct(data) && isstruct(meta), ...
  'devSaveLoadExpSnapshot:MetaOnlyStruct', 'metaOnly should still return data/meta structs.');
localAssert(isempty(loadedVar), ...
  'devSaveLoadExpSnapshot:MetaOnlyLoadedVar', 'loadedVar should be empty in metaOnly mode.');
localAssert(~isempty(inventory), ...
  'devSaveLoadExpSnapshot:MetaOnlyInventory', 'inventory should still be available in metaOnly mode.');

detail = sprintf('metaOnly selectionMode=%s', meta.selectionMode);
end


function detail = localTestLoadToCaller(sessionDir)
outDir = fullfile(sessionDir, 'loadCaller');
mkdir(outDir);
resultFile = localCreateRestoreSnapshot(outDir, 'devSaveLoad_caller');
clear restoredA restoredB;
loadExpSnapshot(resultFile, 'caller', struct('varNames', {{'restoredA', 'restoredB'}}));
localAssert(exist('restoredA', 'var') == 1 && restoredA == 101, ...
  'devSaveLoadExpSnapshot:CallerA', 'restoredA was not restored into caller workspace.');
localAssert(exist('restoredB', 'var') == 1 && isequal(restoredB, [3 5 7]), ...
  'devSaveLoadExpSnapshot:CallerB', 'restoredB was not restored into caller workspace.');

detail = 'caller restore verified';
end


function detail = localTestLoadToBase(sessionDir)
outDir = fullfile(sessionDir, 'loadBase');
mkdir(outDir);
resultFile = localCreateRestoreSnapshot(outDir, 'devSaveLoad_base');

if evalin('base', 'exist(''restoredA'', ''var'')')
  evalin('base', 'clear restoredA');
end
if evalin('base', 'exist(''restoredB'', ''var'')')
  evalin('base', 'clear restoredB');
end

loadExpSnapshot(resultFile, 'base', struct('varNames', {{'restoredA', 'restoredB'}}));
baseA = evalin('base', 'restoredA');
baseB = evalin('base', 'restoredB');
localAssert(baseA == 101, ...
  'devSaveLoadExpSnapshot:BaseA', 'restoredA was not restored into base workspace.');
localAssert(isequal(baseB, [3 5 7]), ...
  'devSaveLoadExpSnapshot:BaseB', 'restoredB was not restored into base workspace.');

evalin('base', 'clear restoredA restoredB');
detail = 'base restore verified';
end


function detail = localTestLegacySnapshot(sessionDir)
outDir = fullfile(sessionDir, 'legacy');
mkdir(outDir);
legacyFile = fullfile(outDir, 'legacySnapshot.mat');

data = struct();
data.legacyA = 8;
data.legacyB = [1 2 4 8];
meta = struct('tag', 'legacyOnly');
save(legacyFile, 'data', 'meta');

[~, metaOut, loadedVar, inventory] = loadExpSnapshot(legacyFile, 'none', struct());
localAssert(numel(loadedVar) == 2, ...
  'devSaveLoadExpSnapshot:LegacyLoadedVar', 'Legacy snapshot should expose both saved variables.');
localAssert(numel(inventory) == 2 && all([inventory.isSaved]), ...
  'devSaveLoadExpSnapshot:LegacyInventory', 'Legacy inventory fallback failed.');
localAssert(isfield(metaOut, 'tag') && strcmp(metaOut.tag, 'legacyOnly'), ...
  'devSaveLoadExpSnapshot:LegacyMeta', 'Legacy meta should remain available.');

detail = 'legacy compatibility verified';
end


function resultFile = localCreateRestoreSnapshot(outDir, prefix)
restoredA = 101; %#ok<NASGU>
restoredB = [3 5 7]; %#ok<NASGU>
helperNoise = 'ignoreMe'; %#ok<NASGU>
resultFile = saveExpSnapshot(prefix, struct( ...
  'includeVars', {{'restoredA', 'restoredB'}}, ...
  'outputDir', outDir, ...
  'verbose', false));
end


function err = localExpectError(fcn)
try
  fcn();
  error('devSaveLoadExpSnapshot:ExpectedError', 'Expected an error, but none was thrown.');
catch err
  if strcmp(err.identifier, 'devSaveLoadExpSnapshot:ExpectedError')
    rethrow(err);
  end
end
end


function localAssert(condition, errId, errMsg)
if ~condition
  error(errId, errMsg);
end
end


function repoRoot = localRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end


function pathStr = localCanonicalPath(pathStr)
pathStr = char(string(pathStr));
if isempty(pathStr)
  return;
end
try
  pathStr = char(java.io.File(pathStr).getCanonicalPath());
catch
  pathStr = char(string(pathStr));
end
pathStr = strrep(pathStr, '/', filesep);
while numel(pathStr) > 1 && (pathStr(end) == '/' || pathStr(end) == '\')
  pathStr(end) = [];
end
end
