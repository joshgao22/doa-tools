function runDynamicRepresentativeReplayReport(repTable, rerunBundleFn, flowOpt, enableWeightSweep, runtimeOpt)
%RUNDYNAMICREPRESENTATIVEREPLAYREPORT Replay representative repeats and print diagnostics.

arguments
  repTable table
  rerunBundleFn function_handle
  flowOpt struct
  enableWeightSweep (1, 1) logical = false
  runtimeOpt struct = struct()
end

maxRepresentativeRepeats = localGetRuntimeOpt(runtimeOpt, 'maxRepresentativeRepeats', height(repTable));
enableTimingSummary = localGetRuntimeOpt(runtimeOpt, 'enableTimingSummary', true);
if isempty(repTable) || maxRepresentativeRepeats <= 0
  return;
end
repTable = repTable(1:min(height(repTable), maxRepresentativeRepeats), :);

for iRep = 1:height(repTable)
  repLabel = repTable.label(iRep);
  repRepeatIdx = repTable.repeatIdx(iRep);
  repTaskSeed = repTable.taskSeed(iRep);
  repBundle = rerunBundleFn(repTaskSeed);

  estTable = buildDoaDopplerSummaryTable(repBundle.caseResult, repBundle.truth, struct('mode', 'dynamic'));
  ablationTruth = [ ...
    buildDoaDopplerCaseTruthFromScene(repBundle.truth, repBundle.periodicFixture.sceneRefOnly), ...
    buildDoaDopplerCaseTruthFromScene(repBundle.truth, repBundle.periodicFixture.sceneOtherOnly)];
  ablationTable = buildDoaDopplerCaseSummaryTable( ...
    [repBundle.caseStaticRefAbl, repBundle.caseStaticOtherOnly], ablationTruth);

  weightTable = table();
  weightSolvePathTable = table();
  if enableWeightSweep
    weightTable = buildDoaDopplerWeightSweepTable(flowOpt.weightSweepAlpha, repBundle.weightCase, repBundle.truth);
    weightSolvePathTable = buildDynamicWeightSolvePathDiagTable( ...
      repBundle.caseMsDoa, repBundle.caseStaticRefAbl, repBundle.weightCase, ...
      flowOpt.weightSweepAlpha, repBundle.truth);
  end

  fprintf('\n========== Representative repeat: %s (repeat %d, taskSeed %d) ==========%s', ...
    repLabel, repRepeatIdx, repTaskSeed, newline);
  disp(estTable);

  fprintf('---------- Static ablation ----------%s', newline);
  disp(ablationTable);

  if enableWeightSweep
    fprintf('---------- Static weight-sweep summary ----------%s', newline);
    disp(weightTable);
  end

  fprintf('---------- Dynamic tooth-selection summary ----------%s', newline);
  disp(repBundle.toothFlowTable);
  fprintf('  selected subset label      : %s\n', repBundle.selectedSubsetLabel);
  fprintf('  selected subset offsets    : %s\n', localFormatIntegerRow(repBundle.selectedSubsetOffsetIdx));
  fprintf('  selected final tag         : %s\n', repBundle.selectedFinalTag);
  if isfield(repBundle, 'finalSelectReason') && strlength(repBundle.finalSelectReason) > 0
    fprintf('  final-select reason        : %s\n', repBundle.finalSelectReason);
  end
  fprintf('  in-tooth fdRange           : %s\n', localFormatNumericRow(repBundle.fdRangeInTooth));
  fprintf('  in-tooth fdRateRange       : %s\n', localFormatNumericRow(repBundle.fdRateRangeInTooth));
  inToothRunMeta = localGetFieldOrDefault(repBundle, 'inToothRunMeta', struct());
  inToothDoaHalfWidthDeg = localGetFieldOrDefault(inToothRunMeta, 'doaHalfWidthDeg', flowOpt.inToothDoaHalfWidthDeg(:).');
  fprintf('  in-tooth DoA half width    : %s\n', localFormatNumericRow(inToothDoaHalfWidthDeg));
  inToothDoaSeedSource = string(localGetFieldOrDefault(inToothRunMeta, 'doaSeedSource', "selected-subset"));
  fprintf('  in-tooth DoA seed source   : %s\n', inToothDoaSeedSource);

  if ~isempty(repBundle.subsetCandidateTable)
    fprintf('---------- Subset candidate summary ----------%s', newline);
    disp(repBundle.subsetCandidateTable);
  end

  if isfield(repBundle, 'finalSelectTable') && ~isempty(repBundle.finalSelectTable)
    fprintf('---------- Dynamic final-select summary ----------%s', newline);
    disp(repBundle.finalSelectTable);
  end

  if enableWeightSweep
    fprintf('---------- Static weight solve-path diagnostics ----------%s', newline);
    disp(weightSolvePathTable);
  end

  fprintf('---------- Dynamic estimator diagnostics ----------%s', newline);
  disp(repBundle.dynSummary.diagTable);

  fprintf('---------- Dynamic objective summary ----------%s', newline);
  disp(repBundle.dynSummary.objectiveTable);

  fprintf('---------- Dynamic per-satellite summary ----------%s', newline);
  disp(repBundle.dynSummary.perSatTable);

  fprintf('---------- Dynamic block summary ----------%s', newline);
  disp(repBundle.dynSummary.blockTable);

  fprintf('---------- Dynamic local-state compare ----------%s', newline);
  disp(repBundle.dynSummary.localStateTable);

  fprintf('---------- Dynamic truth-state compare ----------%s', newline);
  disp(repBundle.dynSummary.truthStateTable);

  fprintf('---------- Dynamic satellite-order summary ----------%s', newline);
  disp(repBundle.dynSummary.satOrderTable);

  if enableTimingSummary && isfield(repBundle, 'timingTable') && ~isempty(repBundle.timingTable)
    fprintf('---------- Dynamic stage timing (sec) ----------%s', newline);
    disp(repBundle.timingTable);
    if isfield(repBundle, 'subsetTimingTable') && ~isempty(repBundle.subsetTimingTable)
      fprintf('---------- Subset candidate timing (sec) ----------%s', newline);
      disp(repBundle.subsetTimingTable);
    end
  end

  if ~isempty(repBundle.dynMultiStartTable)
    fprintf('---------- Dynamic multi-start summary ----------%s', newline);
    disp(repBundle.dynMultiStartTable);
  end
end
end

function textOut = localFormatIntegerRow(valueRow)
%LOCALFORMATINTEGERROW Format one integer row for compact diagnostics.

if isempty(valueRow)
  textOut = '[]';
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), valueRow, 'UniformOutput', false), ', '), ']'];
end

function textOut = localFormatNumericRow(valueRow)
%LOCALFORMATNUMERICROW Format one numeric row for compact diagnostics.

if isempty(valueRow)
  textOut = '[]';
  return;
end
valueRow = reshape(valueRow, 1, []);
textOut = ['[', strjoin(arrayfun(@(x) sprintf('%.6g', x), valueRow, 'UniformOutput', false), ', '), ']'];
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field with a default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function value = localGetRuntimeOpt(runtimeOpt, fieldName, defaultValue)
%LOCALGETRUNTIMEOPT Read one runtime option with default fallback.

value = defaultValue;
if isempty(runtimeOpt) || ~isstruct(runtimeOpt)
  return;
end
if isfield(runtimeOpt, fieldName) && ~isempty(runtimeOpt.(fieldName))
  value = runtimeOpt.(fieldName);
end
end
