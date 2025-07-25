% This script demonstrates 3D weighted reconstruction using the cone-beam DF-CT data acquired by us.
%% add paths
addpath("../functions/utils/");
addpath("../functions/recon/");
addpath("../functions/forward_projection/");
addpath("../functions/unweighted_recon/");addpath("../functions/weighted_iterative/");
%% load data
addpath("..\data\real_data\CT_sino3_0423\")
fid = fopen("df_sino3_gap_bh.raw",'r');df_sino3=fread(fid,'double');df_sino3 = reshape(df_sino3,[35 1200 1200]);fclose(fid);
%% permute 3D sinogram
df_sino3 = permute(df_sino3,[2 1 3]);
%% callibrate detector center
offset = 2; 
df_sino3 = df_sino3([offset+1:end 1:offset],:,:);
%% rebinning
% rebin_factor = [2,2];
% df_sino3 = imresize3(df_sino3,size(df_sino3)./[rebin_factor  1]);

%% define system geometry
clear sys_geo
sys_geo.voxel_num = 512; % voxel number of FOV descretization in x-y plane
sys_geo.phantom_length = 250; % the size of FOV in x-y plane, unit: mm
sys_geo.sid0 = 937;  % source-isocenter distance, unit: mm
sys_geo.sdd0 = 1217; % source-detector distance, unit: mm
sys_geo.pixel_num = size(df_sino3,1); % pixel number of the detector along u direction
sys_geo.detector_length = size(df_sino3,1) * 300e-3; % total length of the detector along u direction, mm. pixel size in u: 300 um
sys_geo.pixel_num3 = size(df_sino3,2); % pixel number of the detector along v direction
sys_geo.detector_length3 = size(df_sino3,2) * 300e-3; % total length of the detector along v direction, mm. pixel size in v: 300 um
sys_geo.voxel_num3 = 20; % voxel number of FOV descretization in z direction
sys_geo.phantom_length3 = sys_geo.voxel_num3 * sys_geo.phantom_length / sys_geo.voxel_num; % the size of FOV in z direction, unit: mm
sys_geo = update_sys_geo(sys_geo);

angle_num = 1200;
theta_list = linspace(0,360,angle_num+1); theta_list(end) = [];

%% reconstruction
% direct FDK
% recon_icone = myIFDKPlanar(df_sino_gap_offset,sys_geo,'theta_list',[0:1:359]);
recon_astra = astraIFDK(df_sino3,sys_geo,theta_list);
% weighted analytical reconstruction
weight_function_propagate = @(idd)(@(x,y,theta)(((y*cos(theta)-x*sin(theta))/idd)+1));

recon_result = myWeightedIConebeamPlanar(df_sino3,sys_geo,'weight_function',weight_function_propagate(sys_geo.idd),'theta_list',theta_list);
