% This script demonstrates 2D weighted reconstruction using the DF-CT data acquired by us.
%% add paths
addpath("../functions/utils/");
addpath("../functions/recon/");
addpath("../functions/unweighted_recon/");

%% load data
load("../data/real_data/CT_sino_0423.mat",'df_sino_gap_bh','tran_sino_gap');
%% callibrate
df_sino = df_sino_gap_bh;
tran_sino = tran_sino_gap;
% the loaded data have already been calibrated for ring artefact and beam-hardening

%% define system geometry
clear sys_geo
sys_geo.voxel_num = 512; % voxel number of FOV descretization
sys_geo.phantom_length = 250; % the size of FOV, unit: mm

sys_geo.sid0 = 937;  % source-isocenter distance, unit: mm
sys_geo.sdd0 = 1217; % source-detector distance, unit: mm
sys_geo.pixel_num = size(df_sino,1); % pixel number of the detector
sys_geo.detector_length = size(df_sino,1) * 300e-3; % total length of the detector, mm. pixel size: 300 um
sys_geo = update_sys_geo(sys_geo);

angle_num = 1200;
theta_list = linspace(0,360,angle_num+1); theta_list(end) = [];

% callibrate the detector center
offset = 2;
% tran_sino = tran_sino([offset+1:end 1:offset],:);
df_sino = df_sino([offset+1:end 1:offset],:);

%% reconstruction
% direct FBP
% note: In this `ifanbeam` function in MATLAB, the sinogram is rebinned into parallel geometry and only 0~pi is used. Therefore we recommend other recon function such as ASTRA toolbox.
recon_ifan = ifanbeam(flipud(df_sino(:,[1+round(size(df_sino,2)/2):end 1:round(size(df_sino,2)/2)])),...
    sys_geo.sid,'fanrotationincrement',360/size(df_sino,2),'fansensorgeometry','line','fansensorspacing',sys_geo.p2v*sys_geo.sid/sys_geo.sdd,'outputsize',sys_geo.voxel_num);

recon_astra = flipud(astraIFanbeam(df_sino,sys_geo,theta_list));


% weighted analytical reconstruction
weight_function_propagate = @(idd)(@(x,y,theta)(((y*cos(theta)-x*sin(theta))/idd)+1));
theta_list = linspace(0,360,angle_num+1); theta_list(end) = [];

recon_result = myWeightedIFanbeamPlanar(df_sino,sys_geo,'weight_function',weight_function_propagate(sys_geo.idd),'theta_list',theta_list);


% conjugate averaging preprocessing
recon_compensate = conjuCompenPlanar(df_sino,sys_geo,'theta_list',theta_list);


