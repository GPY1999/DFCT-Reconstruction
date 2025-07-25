% This script demonstrates 3D weighted reconstruction using the helical DF-CT data acquired by us.
% Unfortunately the dark-field noise of helical data is rather strong.
%% add paths
addpath("../functions/utils/");
addpath("../functions/recon/");
addpath("../functions/unweighted_recon/");
%% load data
load('../data/real_data/CT_sino3_1124.mat','df_sino3');
%% permute 3D sinogram
df_sino3 = df_sino3(21:end,:,:);
df_sino3_p = permute(df_sino3,[2 1 3]);
%% callibrate
% bad gap
gap_list = [595, 606, 626];
df_sino3 = df_sino3_p;
for igap = 1 : length(gap_list)
    gap = gap_list(igap);
    df_sino3(gap,:,:) = 2/3*df_sino3_p(gap-1,:,:) + 1/3*df_sino3_p(gap+2,:,:);
    df_sino3(gap+1,:,:) = 2/3*df_sino3_p(gap+2,:,:) + 1/3*df_sino3_p(gap-1,:,:);
end

% hangle the nan
df_sino3 = handle_nan_inf3(df_sino3);

% callibrate detector center
offset = 0;
% tran_sino = tran_sino([offset+1:end 1:offset],:);
df_sino3 = df_sino3([offset+1:end 1:offset],:,:);


%% define system geometry
rotation_num = 3; % number of rotations
clear sys_geo

sys_geo.voxel_num = 512; % descretize the FOV into 512 voxels in x-y plane
sys_geo.phantom_length = 250; % the size of FOV in x-y plane, unit: mm
sys_geo.sid0 = 939;  % source-isocenter distance, unit: mm
sys_geo.sdd0 = 1219; % source-detector distance, unit: mm
sys_geo.pixel_num = size(df_sino3,1); % pixel number of the detector along u direction
sys_geo.detector_length = size(df_sino3,1) * 300e-3; % total length of the detector along u direction, mm. pixel size in u: 300 um
sys_geo.pixel_num3 = size(df_sino3,2); % pixel number of the detector along v direction
sys_geo.detector_length3 = size(df_sino3,2) * 300e-3; % total length of the detector along v direction, mm. pixel size in v: 300 um
sys_geo.voxel_num3 = 80; % descretize the FOV into 80 voxels in z direction
sys_geo.phantom_length3 = sys_geo.voxel_num3 * sys_geo.phantom_length / sys_geo.voxel_num; % the size of FOV in z direction, unit: mm

sys_geo.pitch0 = 10.08; % pitch of the helical scan trajectory, unit: mm
sys_geo.z_start0 = sys_geo.pitch0 /2 * rotation_num; % the z position of scanning starting point, unit: mm

sys_geo = update_sys_geo(sys_geo);

theta_list = [0:1440-1]*0.75;

%% recon
% direct FDK
recon_my = myIFDKHelicalPlanar(df_sino3,sys_geo,theta_list)
;% weighted analytical reconstruction
weight_function_propagate = @(idd)(@(x,y,theta)(exp(((y*cos(theta)-x*sin(theta))/idd)+1)));
recon_myweight = myWeightedIHelicalPlanar(df_sino3,sys_geo,theta_list,'weight_function',weight_function_propagate(sys_geo.idd));
%% some functions
function interp_data = handle_nan_inf3(data)
% replace the nan in 2D sinogram by interpolation
[rows, cols, depth] = size(data);
[X, Y] = meshgrid(1:cols, 1:rows);
interp_data = data;
% go through every exposure
for k = 1:depth
    disp(['currently working on ' num2str(k)]);
    slice = data(:, :, k);
    validMask = (~isnan(slice)) & (~isinf(slice));
    if any(validMask(:))
        interp_slice = griddata(X(validMask), Y(validMask), slice(validMask), X, Y, 'linear');
        interp_data(:, :, k) = interp_slice;
    end
end
end
