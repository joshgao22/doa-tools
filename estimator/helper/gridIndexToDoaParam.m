function doaParam = gridIndexToDoaParam(gridMeta, gridIdx)
%GRIDINDEXTODOAPARAM Convert DoA grid indices to continuous parameters.
% gridMeta can be a model or any struct carrying
%   - doaType / type
%   - eciAngleGrid for angle mode
%   - latlonGrid   for lat-lon mode

arguments
  gridMeta (1,1) struct
  gridIdx
end

gridIdx = reshape(gridIdx, 1, []);

if isfield(gridMeta, 'doaType')
  doaType = string(gridMeta.doaType);
elseif isfield(gridMeta, 'type')
  doaType = string(gridMeta.type);
else
  error('gridIndexToDoaParam:MissingDoaType', ...
    'gridMeta must contain field doaType or type.');
end

switch lower(char(doaType))
  case 'angle'
    if ~isfield(gridMeta, 'eciAngleGrid')
      error('gridIndexToDoaParam:MissingAngleGrid', ...
        'Angle mode requires field eciAngleGrid.');
    end
    doaParam = gridMeta.eciAngleGrid(:, gridIdx);

  case 'latlon'
    if ~isfield(gridMeta, 'latlonGrid')
      error('gridIndexToDoaParam:MissingLatlonGrid', ...
        'Lat-lon mode requires field latlonGrid.');
    end
    doaParam = gridMeta.latlonGrid(:, gridIdx);

  otherwise
    error('gridIndexToDoaParam:InvalidDoaType', ...
      'Unsupported doa type: %s.', doaType);
end
end
