function resultFile = saveExpSnapshot(prefix)
%SAVE_EXPERIMENT_SNAPSHOT Save current function workspace and entry code into one MAT file.
%
%   resultFile = save_experiment_snapshot()
%   resultFile = save_experiment_snapshot(prefix)
%
% This function packs:
%   1) All variables in the *caller* function workspace
%   2) The entry m-file source code (as text)
% into a single .mat file for experiment reproducibility.
%
% Input:
%   prefix (optional) - filename prefix (default: 'exp_snapshot')
%
% Output:
%   resultFile - saved .mat filename

  if nargin < 1 || isempty(prefix)
    prefix = 'exp_snapshot';
  end

  % ---------- filename ----------
  timestamp  = datestr(now,'yyyymmdd-HHMMSS');
  resultFile = sprintf('%s_%s.mat', prefix, timestamp);

  % ---------- collect caller workspace ----------
  vars = evalin('caller', 'who');
  data = struct();
  for k = 1:numel(vars)
    data.(vars{k}) = evalin('caller', vars{k});
  end

  % ---------- meta information ----------
  meta = struct();
  meta.time  = datetime('now');
  meta.entry = evalin('caller', 'mfilename(''fullpath'')');

  try
    meta.code = fileread([meta.entry '.m']);
  catch
    meta.code = '';
  end

  % ---------- save ----------
  save(resultFile, 'data', 'meta');

  fprintf('Saved experiment snapshot to %s\n', resultFile);
end
