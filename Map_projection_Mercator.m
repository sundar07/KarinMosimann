function Map_projection_Mercator()
% This code peforms map projection (Mercator) of fused spherical data
% Pre-requisites - Image files in ics/ids format and sphere parameters (radius and center) that best fits the sample
% Created by Sundar Naganathan - 2018

% Go to the folder, where all images are stored in ics/ids format
% You should have already obtained a point cloud representing your sample and performed a sphere fit

Map_parameters.mfile = 'Map_projection_Mercator'; % Record the name of the mfile
Bin = 4; % By how much the data was binned

Folder_path = cd; % get the Folder_Folder_path of the current folder

% Make a new folder called 'Maps_eqdcylin' under 'Sphere_fit' if it does not exist already
if exist(strcat(Folder_path,'/binned/Sphere_fit/Maps_eqdcylin'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/binned/Sphere_fit/Maps_eqdcylin'))
end

load(strcat(Folder_path,'/binned/Sphere_fit/Sphere_params.mat')) % Load sphere parameters mat file

% Centre and radius of sphere obtained from 'Get_embryo_surface' code
Map_parameters.Centre_sphere = round(Sphere_params.Center(1,:) .* Bin);
Map_parameters.Radius_sphere = round(Sphere_params.Radius(1) .* Bin);

Map_parameters.Recenter_map = 160; % Depending on region of interest, recenter the map
% By how much the map has to be scaled
Map_parameters.Scale_factor = 5;

% Latitude longitude limits. These values can be changed according to your specific region of interest in the sample
Map_parameters.Lat_start = 90;
Map_parameters.Lat_end = -90;

Map_parameters.Long_start = 0;
Map_parameters.Long_end = 360;

% Projections will be performed for spheres of different radii. Each sphere will be 'radius_step' units apart from each other
Map_parameters.Radius_step = 2;

% Usually I run a dummy map projection by putting arbitrary values in radius_iterate surrounding the sphere radius to know the extent
% to which I want to project and then enter the desired values
Map_parameters.Radius_iterate = 264:Map_parameters.Radius_step:370;

% Read ics files and sort them in alphabetical order
clear Num_chars
Ch0 = dir(strcat(Folder_path,'/*ch_0.ids'));
Ch0_names = {Ch0.name}';
for i = 1:numel(Ch0)
    Num_chars(i) = numel(Ch0_names{i});
end
[~,idx] = sort(Num_chars);
Ch0_names = Ch0_names(idx);

%% Map projection

% x_map,y_map are the x,y-coordinates in the projected map respectively
Start_x = Map_parameters.Long_start*Map_parameters.Scale_factor; % Multiply start/end of latitudes and longitudes by scale factor to get the size
End_x = (Map_parameters.Long_end*Map_parameters.Scale_factor)-1;
Start_y = Map_parameters.Lat_start*Map_parameters.Scale_factor; 
End_y = (Map_parameters.Lat_end*Map_parameters.Scale_factor)-1;

x_map = (Start_x:1:End_x);
y_map = (Start_y:-1:End_y);

% Get lambda (longitude) and phi (latitude) that corresponds to each position in the projected map. The formula will vary depending on the desired projection. 
%% Mercator projection
% Inverse formulas for this projection were obtained here: https://mathworld.wolfram.com/MercatorProjection.html

x_map1 = repmat(x_map,length(y_map),1);
y_map1 = repmat(y_map',1,length(x_map));

% Convert to radian
y_map1 = deg2rad(y_map1);
Lambda_rad = deg2rad(Map_parameters.Recenter_map+x_map1./Map_parameters.Scale_factor);
Phi_rad = (atan(sinh(y_map1./Map_parameters.Scale_factor)));

%% Loop through each time point

for frame = 1:numel(Ch0_names)
    
    % I_map will host the projected images. Clear previous instances of I_map and initialize with NaN.
    clear I_map_ch0
    I_map_ch0 = nan(size(Lambda_rad,1),size(Lambda_rad,2),numel(Map_parameters.Radius_iterate));
    
    % Download bfopen from Mathworks to open ics/ids files
    data_ch0 = bfopen(strcat(Folder_path,'/',Ch0_names{frame}));

    data_dummy = cell2mat(data_ch0{1,1}(1));
    
    % Initialize I_ch0, which will contain images to be map projected
    I_ch0 = nan(size(data_dummy,1),size(data_dummy,2),size(data_ch0{1,1},1));
    clear data_dummy
    % Convert cell to mat and to double format and store all images
    for slice = 1:size(data_ch0{1,1},1)
        I_ch0(:,:,slice) = im2double(cell2mat(data_ch0{1,1}(slice)));
    end
    clear data_ch0

    % Loop through each radius iteration and obtain projections
    for Np = 1:numel(Map_parameters.Radius_iterate)

        % Get x,y,z coordinates that correspond to a certain lambda and phi
        % These are standard spherical to cartesian coordinate conversions
        x = round(Map_parameters.Radius_iterate(Np) .* sin(Lambda_rad) .* cos(Phi_rad) + Map_parameters.Centre_sphere(1));
        y = round(Map_parameters.Centre_sphere(2) - Map_parameters.Radius_iterate(Np) .* sin(Phi_rad));
        z = round(Map_parameters.Centre_sphere(3) - Map_parameters.Radius_iterate(Np) .* cos(Lambda_rad) .* cos(Phi_rad));
        
        % Ensure that there are no indices less than zero or greater than the size of the desired image
        xidx = find(x <= 0 | x > size(I_ch0,2));
        x(xidx) = 1;
        y(xidx) = 1;
        z(xidx) = 1;

        yidx = find(y <= 0 | y > size(I_ch0,1));
        y(yidx) = 1;
        x(yidx) = 1;
        z(yidx) = 1;

        zidx = find(z <= 0 | z > size(I_ch0,3));
        z(zidx) = 1;
        x(zidx) = 1;
        y(zidx) = 1;

        % Get pixel values of the x,y,z coordinates and store them in map projection variables
        I_map_ch0(:,:,Np) = I_ch0(sub2ind(size(I_ch0),y,x,z));
        
        % In case there were indices less than zero or greater than the size of the desired image, convert them back to zero
        I_map_ch0_dummy = I_map_ch0(:,:,Np);
        I_map_ch0_dummy([xidx;yidx;zidx]) = 0;
        I_map_ch0(:,:,Np) = I_map_ch0_dummy;
        
        % You need the 'Mapping' toolbox to run the following 5 lines. If you do not have this toolbox, you can comment it out
        % Visualize the map in classic map style
        figure,axesm('mercator','frame','on','grid','on','glinestyle','-','gcolor',[0 1 1],'ParallelLabel','on','MeridianLabel','on','FLatLimit',[-45 45],'FLonLimit',[20 180])
        geoshow(Phi,Lambda,I_map_ch0(:,:,26),'DisplayType','texturemap')
        set(gcf,'Position',[10 10 2000 1000])
        set(gcf,'Color','w')
        colormap(hot(100))
        
        % Save maps in tif format
        A = I_map_ch0(:,:,Np);
        imwrite(A,strcat(Folder_path,'/binned/Sphere_fit/Maps_eqdcylin/I_map_ch0_fr_',num2str(frame),'.tif'),'WriteMode','append','Compression','none')

    end
    % Save maps in mat format
    save(strcat(Folder_path,'/binned/Sphere_fit/Maps_eqdcylin/I_map_ch0_fr_',num2str(frame),'.mat'),'I_map_ch0','-v7.3')
    
    % Download cprintf function from Mathworks
    cprintf('comment',['Map projection for frame ' num2str(frame) ' completed\n'])
    clear I_ch0 I_ch1
end
% Save parameters that were used for performing map projection
save(strcat(Folder_path,'/binned/Sphere_fit/Maps_eqdcylin/Map_parameters.mat'),'Map_parameters')

