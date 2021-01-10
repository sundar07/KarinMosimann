function Compare_embryo_fit_sphere_over_time()
% This code compares embryo fit sphere over time
% Pre-requisites - You should have already run the sphere fit code and obtained all point cloud mat files
% Created by Sundar Naganathan - 2018

Folder_path = cd;
Sphere_params.Frame_rate = 5; % min

% Get names of all point cloud mat files and sort them
Matfiles = dir(strcat(Folder_path,'/binned/Sphere_fit/*ptCloud.mat'));
Matfiles_names = {Matfiles.name}';
for i = 1:numel(Matfiles)
    Num_chars(i) = numel(Matfiles_names{i});
end
[~,idx] = sort(Num_chars);
Matfiles_names = Matfiles_names(idx);

% Loop through each point cloud mat file and obtain the radius and centre of the sphere
for i = 1:numel(Matfiles)
    load(strcat(Folder_path,'/binned/Sphere_fit/',Matfiles_names{i})) % Load mat file 
    
    Sphere_params.Radius(i) = Pt_Cloud.Radii(1);
    Sphere_params.Center(i,:) = Pt_Cloud.Center;
    
end

% Plot sphere radius over time
Sphere_params.Time = 0:Sphere_params.Frame_rate:(numel(Matfiles)-1)*Sphere_params.Frame_rate;
figure(1);clf;plot(Sphere_params.Time,Sphere_params.Radius,'Marker','o','LineStyle','--','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0.75 0.75])
hold on
plot(Sphere_params.Time,repmat(nanmean(Sphere_params.Radius),1,length(Sphere_params.Time)),'--r')
set(gcf,'Color','w')
% set(gca,'XLim',[-Sphere_params.Frame_rate Sphere_params.Time(end)+Sphere_params.Frame_rate],'YLim',[280 320],'YTick',[280 290 300 round(nanmean(Sphere_params.Radius)) 310 320])
saveas(gcf,strcat(Folder_path,'/binned/Sphere_fit/Fit_sphere_radii'),'tif')

save(strcat(Folder_path,'/binned/Sphere_fit/Sphere_params.mat'),'Sphere_params')