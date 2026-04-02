function [tleSel, outPath] = filterTleByIndex(tleInput, satIdx, outPath)
%FILTERTLEBYINDEX Select satellites from a TLE catalog by index.
% Extracts a subset of TLE entries according to the specified indices.
% When an output file path is provided together with an input TLE file path,
% the selected TLE entries are also written to a new filtered TLE file.
%
%Syntax:
%  tleSel = filterTleByIndex(tleInput, satIdx)
%
%  [tleSel, outPath] = filterTleByIndex(tleInput, satIdx, outPath)
%
%Inputs:
%  tleInput          - TLE source, supported forms:
%                      - char / string scalar : input TLE file path
%                      - struct vector        : TLE structure from tleread
%
%  satIdx            - selected satellite indices
%                      - vector : indices into the original TLE order
%                      - logical: logical mask with length equal to numSat
%
%  outPath           - (optional) output filtered TLE file path
%                      This is only supported when tleInput is a TLE file
%                      path, because exact TLE text lines are preserved from
%                      the original file.
%
%Outputs:
%  tleSel            - filtered TLE structure vector
%
%  outPath           - output filtered TLE file path, or [] if not written
%
%Notes:
%  - If tleInput is a file path, the returned tleSel is obtained through
%    tleread on the original file and then indexed by satIdx.
%  - If outPath is specified, the function copies the original 3-line TLE
%    blocks (name line + line 1 + line 2) for the selected satellites.
%  - For downstream propagation only, keeping tleSel in memory is enough;
%    you do not need to write a new .tle file.
%
%Example:
%  tleAll = tleread("./tle/starlink_all.tle");
%  availIdx = satIdx.available;
%  tleAvail = filterTleByIndex(tleAll, availIdx);
%
%  [tleAvail, outFile] = filterTleByIndex("./tle/starlink_all.tle", ...
%    availIdx, "./tle/starlink_available.tle");
%
%See also:
%  tleread, propagateOrbit, findVisibleSatFromTle

arguments
  tleInput
  satIdx
  outPath = []
end

% -------------------------------------------------------------------------
% Parse TLE input
% -------------------------------------------------------------------------
[parsedTle, srcPath] = localParseTleInput(tleInput);
numSat = numel(parsedTle);

% -------------------------------------------------------------------------
% Parse and validate selected indices
% -------------------------------------------------------------------------
satIdx = localParseIndex(satIdx, numSat);
tleSel = parsedTle(satIdx);

% -------------------------------------------------------------------------
% Write filtered TLE file if requested
% -------------------------------------------------------------------------
if isempty(outPath)
  outPath = [];
  return;
end

if isempty(srcPath)
  error('filterTleByIndex:WriteRequiresFileInput', ...
    ['Writing a filtered TLE file requires tleInput to be the original ', ...
     'TLE file path, so the exact TLE text lines can be preserved.']);
end

outPath = char(string(outPath));
localWriteFilteredTle(srcPath, satIdx, outPath);
end

function [parsedTle, srcPath] = localParseTleInput(tleInput)
%LOCALPARSETLEINPUT Parse TLE file path or struct input.

srcPath = '';

if ischar(tleInput) || (isstring(tleInput) && isscalar(tleInput))
  srcPath = char(string(tleInput));
  if exist(srcPath, 'file') ~= 2
    error('filterTleByIndex:TleFileNotFound', ...
      'The TLE file does not exist: %s', srcPath);
  end
  parsedTle = tleread(srcPath);
  parsedTle = parsedTle(:);
  return;
end

if isstruct(tleInput)
  parsedTle = tleInput(:);
  return;
end

error('filterTleByIndex:UnsupportedTleInput', ...
  'tleInput must be a TLE file path or a struct vector from tleread.');
end

function satIdx = localParseIndex(satIdx, numSat)
%LOCALPARSEINDEX Parse numeric indices or logical mask.

if islogical(satIdx)
  if numel(satIdx) ~= numSat
    error('filterTleByIndex:LogicalIndexSizeMismatch', ...
      'Logical satIdx must have length equal to the number of TLE entries.');
  end
  satIdx = find(satIdx(:)).';
  return;
end

if ~isnumeric(satIdx) || ~isreal(satIdx)
  error('filterTleByIndex:InvalidIndexType', ...
    'satIdx must be a numeric index vector or a logical mask.');
end

satIdx = satIdx(:).';
if isempty(satIdx)
  return;
end

if any(~isfinite(satIdx)) || any(mod(satIdx, 1) ~= 0)
  error('filterTleByIndex:InvalidIndexValue', ...
    'satIdx must contain finite integer indices.');
end

if any(satIdx < 1) || any(satIdx > numSat)
  error('filterTleByIndex:IndexOutOfRange', ...
    'satIdx must lie within [1, %d].', numSat);
end

satIdx = unique(satIdx, 'stable');
end

function localWriteFilteredTle(srcPath, satIdx, outPath)
%LOCALWRITEFILTEREDTLE Write selected 3-line TLE blocks to a new file.

allLines = localReadTextLines(srcPath);
lineBlock = localParseTleLineBlocks(allLines);
numBlock = numel(lineBlock);

if numBlock < max([satIdx, 0])
  error('filterTleByIndex:TleBlockCountMismatch', ...
    ['The number of parsed TLE blocks in the source file is smaller than ', ...
     'the requested maximum index.']);
end

outFolder = fileparts(outPath);
if ~isempty(outFolder) && exist(outFolder, 'dir') ~= 7
  mkdirOk = mkdir(outFolder);
  if ~mkdirOk
    error('filterTleByIndex:CreateOutputDirFailed', ...
      'Failed to create output folder: %s', outFolder);
  end
end

fid = fopen(outPath, 'w');
if fid < 0
  error('filterTleByIndex:OpenOutputFailed', ...
    'Failed to open output TLE file: %s', outPath);
end
cleanupObj = onCleanup(@() fclose(fid));

for iSat = 1:numel(satIdx)
  curBlock = lineBlock{satIdx(iSat)};
  for iLine = 1:numel(curBlock)
    fprintf(fid, '%s\n', curBlock{iLine});
  end
end

clear cleanupObj;
end

function allLines = localReadTextLines(filePath)
%LOCALREADTEXTLINES Read a text file into a cell array of lines.

fid = fopen(filePath, 'r');
if fid < 0
  error('filterTleByIndex:OpenInputFailed', ...
    'Failed to open input TLE file: %s', filePath);
end
cleanupObj = onCleanup(@() fclose(fid));

allLines = {};
while true
  curLine = fgetl(fid);
  if isequal(curLine, -1)
    break;
  end
  allLines{end + 1, 1} = curLine; %#ok<AGROW>
end

clear cleanupObj;
end

function lineBlock = localParseTleLineBlocks(allLines)
%LOCALPARSETLELINEBLOCKS Parse a text TLE file into 3-line blocks.

numLine = numel(allLines);
lineBlock = {};
lineIdx = 1;

while lineIdx <= numLine
  nameLine = strtrim(allLines{lineIdx});

  if isempty(nameLine)
    lineIdx = lineIdx + 1;
    continue;
  end

  if startsWith(nameLine, '1 ') || startsWith(nameLine, '2 ')
    error('filterTleByIndex:MissingNameLine', ...
      ['Invalid TLE text format. Expected a satellite name line before ', ...
       'line 1 / line 2.']);
  end

  if lineIdx + 2 > numLine
    error('filterTleByIndex:IncompleteTleBlock', ...
      'The TLE file ends with an incomplete 3-line block.');
  end

  line1 = strtrim(allLines{lineIdx + 1});
  line2 = strtrim(allLines{lineIdx + 2});

  if ~startsWith(line1, '1 ') || ~startsWith(line2, '2 ')
    error('filterTleByIndex:InvalidTleBlock', ...
      ['Invalid TLE text format. Each satellite must contain one name ', ...
       'line followed by line 1 and line 2.']);
  end

  lineBlock{end + 1, 1} = {nameLine; line1; line2}; %#ok<AGROW>
  lineIdx = lineIdx + 3;
end
end
