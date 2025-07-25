% This script demonstrates numerical simulation with weighted reconstruction.
%% add paths
addpath("../functions/utils/");
addpath("../functions/recon/");
addpath("../functions/forward_projection/");
addpath("../functions/unweighted_recon/");addpath("../functions/weighted_iterative/");

%% define system geometry
angle_num = 1200;
theta_list = linspace(0,360,angle_num+1)+0.0128; theta_list(end) = []; 

clear sys_geo
sys_geo.sid0 = 400; % source-isocenter distance, mm
sys_geo.sdd0 = 600; % source-detector distance, mm
sys_geo.pixel_num = 2048 ; % pixel number of the detector
sys_geo.detector_length = 2048 * 0.2 ; % total length of the detector, mm
sys_geo.voxel_num = 1024; % descretized voxel number of the FOV
sys_geo.phantom_length = 1024 * 0.25; % size of the FOV, mm

sys_geo = update_sys_geo(sys_geo);

%% load digital phantom
addpath("../data/digital_phantom/")
load("advance_320.mat","phant");
phant = interp2(linspace(1,sys_geo.voxel_num,size(phant,1)),linspace(1,sys_geo.voxel_num,size(phant,2)),phant,1:sys_geo.voxel_num,(1:sys_geo.voxel_num)',"nearest");

% if you want another phantom:
% load("Forbild_phantom_mod2.mat",'phant');
% phant = interp2(linspace(1,sys_geo.voxel_num,size(phant,1)),linspace(1,sys_geo.voxel_num,size(phant,2)),phant,1:sys_geo.voxel_num,(1:sys_geo.voxel_num)',"nearest");

%% fanbeam forward projection and reconstruction
idd = sys_geo.sdd-sys_geo.sid; % distance between iso-center and detector

% forward model weight
% exponentia weight
% weight_function_propagate = @(idd)(@(x,y,theta)(exp(((y*cos(theta)-x*sin(theta))/idd).^1)));
% proportional weight
weight_function_propagate = @(idd)(@(x,y,theta)((+((y*cos(theta)-x*sin(theta))/idd)+1)));
% no weight
% weight_function_propagate = @(idd)(@(x,y,theta)(1));

% weight for reconstruction
% use proportional weight
weight_function_recon = @(idd)(@(x,y,theta)((+((y*cos(theta)-x*sin(theta))/idd)+1)));
% use the same weight as forward model
% weight_function_recon = weight_function_propagate;

sino_my = myWeightedFanbeamPlanar(phant,sys_geo,'theta_list',theta_list,'weight_function',weight_function_propagate(idd));

% weighted analytical reconstruction
recon = myWeightedIFanbeamPlanar(sino_my,sys_geo,'weight_function',weight_function_propagate(idd),'theta_list',theta_list);
% weighted analytical reconstruction - using proportional weight
recon_zerointercept = myWeightedIFanbeamPlanar(sino_my,sys_geo,'weight_function',weight_function_recon(idd),'theta_list',theta_list);
% direct FBP
% note: In this `ifanbeam` function in MATLAB, the sinogram is rebinned into parallel geometry and only 0~pi is used. Therefore we recommend other recon function such as ASTRA toolbox.
recon_ifan = ifanbeam(flipud(sino_my(:,[(round(length(theta_list)/2)+1):length(theta_list) 1:round(length(theta_list)/2)])),sys_geo.sid,'fanrotationincrement',theta_list(2)-theta_list(1),'fansensorgeometry','line','fansensorspacing',sys_geo.p2v*sys_geo.sid/sys_geo.sdd,'outputsize',sys_geo.voxel_num);
recon_astra = flipud(astraIFanbeam(sino_my,sys_geo,theta_list));
% conjugate averaging preprocessing
[recon_compen,sino_compen] = conjuCompenPlanar(sino_my,sys_geo,'theta_list',theta_list);
% weighted iterative reconstruction
% recon_artzero = myWeightedARTFanbeamPlanar(sino_my,sys_geo,'theta_list',theta_list,'weight_function',weight_function_recon(idd),'nART',10,'recon_init',recon_astra);
recon_art = myWeightedARTFanbeamPlanar(sino_my,sys_geo,'theta_list',theta_list,'weight_function',weight_function_propagate(idd),'nART',10,'recon_init',recon_astra);

