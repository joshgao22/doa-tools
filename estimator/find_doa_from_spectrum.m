function [est, est_idx, resolved] = find_doa_from_spectrum(dg, y, n, varargin)
%FIND_DOA_FROM_SPECTRUM Identifies DOA estimates by finding largest
%peaks in the given spectrum, supporting both 1D and 2D grids.
%Inputs:
%   dg - DOA grid.
%        For 1D: a 1-dimensional row vector.
%        For 2D: a 2-row matrix, where the first row is azimuth and the 
%               second row is elevation.
%   y - Computed spectrum. Should have the same length as dg (for 1D) or
%       correspond to the flattened 2D grid (for 2D).
%   n - Expected number of DOAs.
%   varargin - Optional arguments:
%       grid_size - Explicit grid size.
%                   1D: scalar (number of points)
%                   2D: [azimuth_points, elevation_points]
%Outputs:
%   est - Estimated DOAs.
%         For 1D: a 1 × n vector.
%         For 2D: a n × 2 matrix, each row is [azimuth, elevation].
%   est_idx - Estimated DOA indices.
%         For 1D: a 1 × n vector of indices.
%         For 2D: a n × 2 matrix, each row [azimuth_idx, elevation_idx].
%   resolved - A boolean indicating whether the direction finding is
%              successful.

% Parse optional arguments
if ~isempty(varargin)
    grid_size = varargin{1};
else
    grid_size = [];
end

% Basic check
if size(dg, 2) ~= length(y)
    error('The size of dg and y must match.');
end

[peaks, idx_est] = findpeaks(y);

if length(idx_est) == n
    % Exactly n peaks found
    [est, est_idx] = format_estimates(dg, idx_est, grid_size);
    resolved = true;
    
elseif length(idx_est) > n
    % More than expected, pick the n largest
    [~, sort_idx] = sort(peaks, 'descend');
    selected_idx = idx_est(sort_idx(1:n));
    [est, est_idx] = format_estimates(dg, selected_idx, grid_size);
    resolved = true;
    
else
    % fewer than expected
    est = [];
    est_idx = [];
    resolved = false;
end
end

function [est, est_idx] = format_estimates(dg, idx, grid_size)
%FORMAT_ESTIMATES Internal function to format estimates based on grid dimension.
switch size(dg,1)
    case 1
        % 1D case
        est = dg(idx);
        est_idx = idx;
    case 2
        % 2D case
        est = dg(:, idx).'; % n × 2
        if isempty(grid_size)
            % If not provided, infer grid size
            [nx, ny] = grid_size_from_dg(dg);
        else
            % User provided grid size
            if isscalar(grid_size)
                nx = grid_size;
                ny = 1;
            else
                nx = grid_size(1);
                ny = grid_size(2);
            end
        end
        [az_idx, el_idx] = ind2sub([nx, ny], idx);
        est_idx = [az_idx, el_idx]; % n × 2
    otherwise
        error('Unsupported grid dimension. dg must have 1 or 2 rows.');
end
end

function [nx, ny] = grid_size_from_dg(dg)
%GRID_SIZE_FROM_DG Estimates the original grid size [nx, ny] from flattened dg.
%Assumptions:
%   - The 2D grid is constructed by grid2().
%   - In grid2(), azimuth varies slowly, elevation varies quickly.
%     (i.e., elevation is the faster-changing dimension.)
%Method:
%   - The number of unique elevation values determines ny.
%   - The total number of points divided by ny gives nx.

elev = dg(2,:); 
unique_el = numel(unique(round(elev, 8))); % safe rounding
ny = unique_el;
nx = size(dg,2) / ny;
end
