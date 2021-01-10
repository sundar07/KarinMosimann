function Get_isosurface()
% This code detects isosurface in fused images and performs a sphere fit of the point cloud
% Created by Sundar Naganathan - 2018

% Go to the folder, where all images are stored
Folder_path = cd; % get the Folder_path of the current folder

% Get names of all tif files to be analyzed
DirOutput = dir(strcat(Folder_path,'/c*.tif'));
FileNames = {DirOutput.name};

% If Sphere_fit folder does not exist, make it
if exist(strcat(Folder_path,'/Sphere_fit'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/Sphere_fit'))
end

% Loop through each image
for i = 1:numel(FileNames)
    
    Pt_Cloud.mfile = 'Get_isosurface'; % Record the name of the mfile
    FileInfo = imfinfo(strcat(Folder_path,'/',FileNames{i}));
    
    Name = strrep(FileNames{i},'.tif','');
    
    % Pre-initalize the variable Img_thresh_open, which will contain all tif files
    Img_thresh_open = nan(FileInfo(1).Height,FileInfo(1).Width,numel(FileInfo));
    
    % Loop through each image, read the images followed by thresholding
    for j = 1:numel(FileInfo)        
        Img = im2double(mat2gray(imread(strcat(Folder_path,'/',FileNames{i}),j)));
        Img_thresh = imbinarize(Img,'adaptive','Sensitivity',0.4); % Perform adaptive thresholding. The sensitivity can be changed according to your image
        Img_thresh_open(:,:,j) = bwareaopen(Img_thresh,500); % Remove small particles
    end
    
    Img_clear_border = imclearborder(Img_thresh_open); % Clear the border of the images, which is usually noisy
%     [L,~] = bwlabeln(Img_thresh_open);
%     figure,imshow(L(:,:,100))
    
    map = bwdist(Img_clear_border); % Compute Euclidean distance transform of the binary image
    surf = isosurface(map,10); % Obtain the isosurface. The isovalue (here, 10) can be changed and played around until you obtain the isosurface that best represents your data

    % Make point cloud from vertices
    p1 = surf.vertices;

    %p2c = p2 - c2 + c1;

    P1 = pointCloud(p1); % Make a point cloud object
    Pt_Cloud.Positions = P1.Location;
%     pcshow(P1)    
    
    % Fit a sphere to the point cloud
    % ellipsoid_fit was downloaded from Mathworks developed by Yury Petrov (Sep 2015)
    [center,radii,evecs,v,chi2] = ellipsoid_fit(P1.Location,'xyz'); 
    Pt_Cloud.Center = center;
    Pt_Cloud.Radii = radii;
    Pt_Cloud.evecs = evecs;
    Pt_Cloud.v = v;
    Pt_Cloud.chi2 = chi2;

    %%%%%%%%% Plot the point cloud and the fit sphere %%%%%%%%%

    x = P1.Location(1:50:end,1);
    y = P1.Location(1:50:end,2);
    z = P1.Location(1:50:end,3);
    figure(2);clf;
    plot3(x,y,z,'Marker','o','MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
    hold on
    %draw fit
    mind = min([x y z])-200; 
    maxd = max([x y z])+200; 
    nsteps = 50;
    step = (maxd - mind) / nsteps;
    [x,y,z] = meshgrid(linspace(mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps));
    % Obtain the ellipsoid values from the fit
    Ellipsoid = v(1)*x.*x + v(2)*y.*y + v(3)*z.*z + ...
    2*v(4)*x.*y + 2*v(5)*x.*z + 2*v(6)*y.*z + ...
    2*v(7)*x + 2*v(8)*y + 2*v(9)*z;
    % Use patch to draw the ellipsoid
    p = patch(isosurface(x,y,z,Ellipsoid,-v(10)));
    hold off
    set(p,'FaceColor','g','EdgeColor','none','FaceAlpha',0.3);
    view(34.4,14.4); % these values can be changed according to how you want to rotate your data 
    axis vis3d equal;
    camlight;
    lighting phong;
    
    % Save image of sphere fit
    saveas(gcf,strcat(Folder_path,'/Sphere_fit/',Name,'_SphereFit'),'tif')
    save(strcat(Folder_path,'/Sphere_fit/',Name,'_ptCloud.mat'),'Pt_Cloud') % Save the point cloud
end