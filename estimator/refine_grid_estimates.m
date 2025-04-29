function est = refine_grid_estimates(f_obj, grid, est_idx)
%REFINE_GRID_ESTIMATES Refines DOA estimates obtained from a grid.
%Inputs:
%   f_obj - Objective function. Its local minimums identify the DOAs.
%   grid - Grid used for the original estimation. For the 1D case, this is 
%          an 1-row vector. For the 2D case, this is a 2-row matrix, where 
%          the first row represents the range of the azimuth angle, and the
%          second row represents that of the elevation angle.
%   est_idx - Indices (corresponding to the grid) of the original
%             estimates.
%Output:
%   est - Refined estimates.
est = zeros(size(est_idx));
if size(grid, 1) == 1
    % 1d
    n_iter = 10;
    subgrid_size = 10;
    for kk = 1:length(est_idx)
        % k-th DOA
        % init bounds
        if est_idx(kk) > 1
            lb = grid(est_idx(kk) - 1);
        else
            lb = grid(1);
        end
        if est_idx(kk) < length(grid)
            ub = grid(est_idx(kk) + 1);
        else
            ub = grid(end);
        end
        % refine
        for ii = 1:n_iter
            % find minimum over the subgrid
            subgrid = linspace(lb, ub, subgrid_size);
            obj_vals = f_obj(subgrid);
            [~, min_idx] = min(obj_vals);
            % update bounds
            if min_idx > 1
                lb = subgrid(min_idx - 1);
            else
                lb = subgrid(1);
            end
            if min_idx < subgrid_size
                ub = subgrid(min_idx + 1);
            else
                ub = subgrid(end);
            end
        end
        est(kk) = subgrid(min_idx);
    end
else
    % 2d
    n_iter = 10;
    subgrid_size = [4 4]; % 2D subgrid size
    
    for kk = 1:size(est_idx,1)
        % k-th DOA
        % init bounds for azimuth and elevation angles
        azimuth_idx = est_idx(kk, 1);
        elevation_idx = est_idx(kk, 2);
        
        if azimuth_idx > 1
            lb_azimuth = grid(1, azimuth_idx - 1);
        else
            lb_azimuth = grid(1, 1);
        end
        if azimuth_idx < size(grid, 2)
            ub_azimuth = grid(1, azimuth_idx + 1);
        else
            ub_azimuth = grid(1, end);
        end
        
        if elevation_idx > 1
            lb_elevation = grid(2, elevation_idx - 1);
        else
            lb_elevation = grid(2, 1);
        end
        if elevation_idx < size(grid, 2)
            ub_elevation = grid(2, elevation_idx + 1);
        else
            ub_elevation = grid(2, end);
        end
        
        % refine
        for ii = 1:n_iter
            % generate 2D subgrid
            azimuth_subgrid = linspace(lb_azimuth, ub_azimuth, subgrid_size(1));
            elevation_subgrid = linspace(lb_elevation, ub_elevation, subgrid_size(2));
            
            % evaluate the objective function on the 2D subgrid
            [azimuth_grid, elevation_grid] = meshgrid(azimuth_subgrid, elevation_subgrid);
            obj_vals = f_obj([azimuth_grid(:), elevation_grid(:)]);
            
            % find the minimum
            [~, min_idx] = min(obj_vals);
            [min_azimuth_idx, min_elevation_idx] = ind2sub(subgrid_size, min_idx);
            
            % update bounds for the next iteration
            lb_azimuth = azimuth_subgrid(min_azimuth_idx);
            ub_azimuth = azimuth_subgrid(min_azimuth_idx + 1);
            lb_elevation = elevation_subgrid(min_elevation_idx);
            ub_elevation = elevation_subgrid(min_elevation_idx + 1);
        end
        
        % store refined estimate
        est(kk, :) = [azimuth_subgrid(min_azimuth_idx), elevation_subgrid(min_elevation_idx)];
    end
end
end

