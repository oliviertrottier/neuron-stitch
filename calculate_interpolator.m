function grid_interpolator = calculate_interpolator(mask,varargin)
% Function that solves the heat equation using the pde toolbox to find a 2D interpolator
% that interpolates boundary values defined in a mask.
% Only linear boundaries or rectangular boundaries are allowed in mask.

% Input
% mask = double array whose entries define the boundary types (1 or 2).

% Determine the boundaries in cartesian coordinates from the mask.
mask_size = size(mask);
mask_interior = mask(2:end-1,2:end-1);
%% Parse optional parameters
p = inputParser;
addParameter(p, 'Plots', nargout==0); % Plot info plots detailing the interpolator calculation.
addParameter(p, 'Method', 'Heat_numeric' ,@(x) ismember(x,{'Heat_numeric','Laplace_theoretical'})); % Method to compute the interpolator. 'Laplace_solution' is faster but works only for rectangular boundaries.
parse(p, varargin{:});
options = p.Results;
%% Define geometric rectangles for each region type.
% Plot the original mask.
if options.Plots
    figure;
    imagesc(mask);colorbar;
    title('Original Mask')
end

Region_types = unique(mask(mask>0));
Region_rectangles = [];
for i=1:numel(Region_types)
    type = Region_types(i);
    region_mask = mask_interior==type;
    if nnz(region_mask) > 0
        region_properties = regionprops(region_mask);
        
        % For each region of type==Region_types(i), find its bounding box.
        for j=1:numel(region_properties)
            BoundingBox = repmat(ceil(region_properties(j).BoundingBox(1:2)),1,2) + [0 0 region_properties(j).BoundingBox(3:4)-1];
            
            % Correct the coordinates for the truncation.
            BoundingBox = BoundingBox + 1;
            BoundingBox(BoundingBox==2)=1;
            BoundingBox(3) = BoundingBox(3) + double((BoundingBox(3)==mask_size(2)-1));
            BoundingBox(4) = BoundingBox(4) + double((BoundingBox(4)==mask_size(1)-1));
            
            % Define the rectangle in geometric form.
            rect = [3 4 BoundingBox(1) BoundingBox(3) BoundingBox(3) BoundingBox(1) BoundingBox(2) BoundingBox(2) BoundingBox(4) BoundingBox(4)];
            Region_rectangles = [rect' Region_rectangles];
        end
    end
end

%% Form the interior region by removing rectangles belonging to each region type.
Initial_rectangle = [3 4 1 mask_size(2) mask_size(2) 1 1 1 mask_size(1) mask_size(1)]';
N_rectangles_subtracted = size(Region_rectangles,2);

rectangle_labels = cellfun(@(x) ['R',num2str(x)],num2cell(1:N_rectangles_subtracted+1),'Uni',0);
set_formula = join(rectangle_labels,'-');
boundaries = decsg([Initial_rectangle Region_rectangles],set_formula{1},cell2mat(rectangle_labels')');
N_boundaries = size(boundaries,2);

% Determine the type of each boundary.
Boundary_type = nan(1,N_boundaries);
for i=1:N_boundaries
    edge_row_ind = round(min(boundaries(4:5,i)):max(boundaries(4:5,i)));
    edge_col_ind = round(min(boundaries(2:3,i)):max(boundaries(2:3,i)));
    edge_values = mask(edge_row_ind,edge_col_ind);
    type_counts = histc(edge_values,Region_types);
    [~,max_count_ind] = max(type_counts);
    Boundary_type(i) = Region_types(max_count_ind);
end

% Plot the interior region.
if options.Plots
    figure;
    Model = createpde;
    geometryFromEdges(Model,boundaries);
    pdegplot(Model,'EdgeLabels','on');
    title('Edge Labels');
    a=gca;a.YDir='reverse';
end

%% Solve for the interpolator
switch options.Method
    case 'Heat_numeric'
        thermalmodelS = createpde('thermal','steadystate');
        geometryFromEdges(thermalmodelS,boundaries);
        
        % Specify conductivity.
        thermalProperties(thermalmodelS,'ThermalConductivity',1);
        
        % Set the temperature for each boundary.
        for i=1:N_boundaries
            thermalBC(thermalmodelS,'Edge',i,'Temperature',Boundary_type(i));
        end
        
        % Create mesh.
        hmin = max(mask_size)/100;
        warning('off','MATLAB:subscripting:noSubscriptsSpecified');
        triangular_mesh = generateMesh(thermalmodelS,'Hmin',hmin,'Hmax',1.5*hmin);
        warning('on','MATLAB:subscripting:noSubscriptsSpecified');
        
        % Solve for the steady state solution.
        R = solve(thermalmodelS);
        T = R.Temperature;
        assert(all(~isnan(T)),'The temperature calculation is erroneous.');
        
        % Interpolate the coarse interpolator to obtain a finer interpolator on a uniform grid.
        [mask_grid_y,mask_grid_x] = ndgrid(mask_size(1):-1:1,1:mask_size(2));
        Interpolant = scatteredInterpolant(triangular_mesh.Nodes(1,:)',triangular_mesh.Nodes(2,:)',T,'linear','linear');
        grid_interpolator = Interpolant(mask_grid_x,mask_grid_y);
        grid_interpolator(mask>0) = mask(mask>0);
        
        if any(isnan(grid_interpolator(:)))
            error('The interpolator has nan values.');
        end
        
        if options.Plots
            figure;
            pdeplot(thermalmodelS,'XYData',T,'Contour','on','ColorMap','hot');
            title('Temperature field')
            axis equal tight
        end
        
    case 'Laplace_theoretical'
        % Check that the boundaries in the mask are on the edge of the mask
        % and 1 pixel wide.
        assert(N_boundaries==4,'The Laplace method can only be used for rectangular boundaries.');
        are_boundaries_on_edge = all(ismember(unique(boundaries(2:3,:)),[1 mask_size(2)])) & all(ismember(unique(boundaries(4:5,:)),[1 mask_size(1)]));
        assert(are_boundaries_on_edge,'In the Laplace method, boundaries must be 1 pixel wide and located at the edge of the mask array.')
        
        % Create an interpolator between the 2D boundaries using solutions of
        % the Laplace equation.
        
        % V(0,y) = 1
        % V(1,y) = 0
        [mask_grid_y,mask_grid_x] = ndgrid(mask_size(1):-1:1,1:mask_size(2));
        
        % Normalize the grid coordinates to (0,1).
        grid_y = (mask_grid_y - 1)/(mask_size(1) - 1);
        grid_x = (mask_grid_x - 1)/(mask_size(2) - 1);
        grid_interpolator = ones(mask_size);
        
        % Add the solution that fixes the horizontal boundary.
        bottom_boundary_type = Boundary_type(all(boundaries(4:5,:)==1,1));
        top_boundary_type = Boundary_type(all(boundaries(4:5,:)==mask_size(1),1));
        if bottom_boundary_type==2
            % V(x,0) = 2
            % V(x,1) = 1
            grid_interpolator = 2/pi*atan(sin(pi*grid_x)./sinh(pi*(grid_y))) + grid_interpolator;
        else
            % V(x,0) = 1
            % V(x,1) = 2
            grid_interpolator = 2/pi*atan(sin(pi*grid_x)./sinh(pi*(1-grid_y))) + grid_interpolator;
        end
        
         % Add the solution that fixes the vertical boundary.
        left_boundary_type = Boundary_type(all(boundaries(2:3,:)==1,1));
        right_boundary_type = Boundary_type(all(boundaries(2:3,:)==mask_size(2),1));
        if left_boundary_type==2
            % V(0,y) = 2
            % V(1,y) = 1
            grid_interpolator = 2/pi*atan(sin(pi*grid_y)./sinh(pi*(grid_x))) + grid_interpolator;
        else
            % V(0,y) = 1
            % V(1,y) = 2
            grid_interpolator = 2/pi*atan(sin(pi*grid_y)./sinh(pi*(1-grid_x))) + grid_interpolator;
        end
        
        % Fix the interpolator at the corners.
        grid_interpolator(1,1) = (left_boundary_type + top_boundary_type)/2;
        grid_interpolator(mask_size(1),1) = (left_boundary_type + bottom_boundary_type)/2;
        grid_interpolator(1,mask_size(2)) = (right_boundary_type + bottom_boundary_type)/2;
        grid_interpolator(mask_size(1),mask_size(2)) = (right_boundary_type + bottom_boundary_type)/2;
end

% Plot the final interpolator.
if options.Plots
    figure;
    imagesc([mask_grid_x(1) mask_grid_x(end)],[mask_grid_y(1) mask_grid_y(end)],grid_interpolator); a=gca; a.YDir='normal';
    colorbar;
    axis equal tight
    title('Interpolator field')
end
end