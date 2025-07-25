function [recon_compensate,sino_compensate] = conjuCompenPlanar(sino,sys_geo,varargin)
%CONJUCOMPENPLANAR Reconstruction with conjugate averaging pre-processing. For fan-beam and planar detector.
%   For every ray in the sinogram, find its conjugate ray and average with it, then normalize to iso-center (this normalization is for planar detector)
%   INPUT:
%     sino: 2-D sinogram, shape: pixel * angle
%     sys_geo: a struct of the system's geometrical parameters, including:
%         pixel_num: number of the pixels of the detector
%         phantom_length: size of the phantom, scalar, unit: mm
%         voxel_num: number of the voxels of the FOV descretization
%         detector_length: length of the detector, unit: mm
%         d2v, v2d, p2d, p2v: the transformation coefficient from unit Distance (mm) to Voxel size, Voxel size to unit Distance, Pixel size to unit Distance, and Pixel size to Voxel size
%         sid0, sdd0, sid, sdd: the distance between Source and Iso-center, and between Source and Detector, unit: mm or voxel
%         
%     varargin: 
%         theta_list: projection angles, 1D vector, unit: degree, default value: [0:359]
% UPDATES:
%   2025/07/25 (GUO Peiyuan): complete comments
%   2024/11/09 (GUO Peiyuan): second version, revised from conjuCompen.m
%   2024/07/05 (GUO Peiyuan): first version


% set varargin
defaults = {'theta_list',[0:359]};
params = parseKeyValuePairs(defaults,varargin{:});
theta_list = params.theta_list;

field_angle = -(atan(((1:sys_geo.pixel_num)-sys_geo.pixel_num/2-1/2)*sys_geo.p2v./sys_geo.sdd));
angle_num = length(theta_list);
angle_increment_rad = 2 * pi / angle_num;
[angle_rad_mesh,field_angle_mesh] = meshgrid(theta_list*pi/180,field_angle);
% search conjugate ray pair
field_angle_conjugate = - field_angle_mesh;
angle_rad_conjugate = mod(angle_rad_mesh + (pi + 2 * abs(field_angle_mesh)) .* ((field_angle_mesh<=0)-1/2)*2 , 2*pi);
angle_rad_extended = [-angle_increment_rad*ones(sys_geo.pixel_num,1), angle_rad_mesh, (angle_num)*angle_increment_rad*ones(sys_geo.pixel_num,1)];
field_angle_extended = [field_angle_mesh(:,[1 1 1:end])];
sino_extended = [sino(:,end) sino sino(:,1)];
sino_conjugate = interp2(angle_rad_extended,field_angle_extended,sino_extended,angle_rad_conjugate,field_angle_conjugate,"linear");
% calculate fan-angle dependent weight
weight_field_angle = (sys_geo.sdd./cos(field_angle_mesh) - (sys_geo.sid)*cos(field_angle_mesh)) / (sys_geo.sdd - sys_geo.sid) ./ cos(field_angle_mesh);
% average and recon
sino_compensate = (sino_conjugate + sino)/2 ./ weight_field_angle;
recon_compensate=ifanbeam(flipud(sino_compensate(:,[angle_num/2+1:end 1:angle_num/2])),sys_geo.sid,'fanrotationincrement',angle_increment_rad/pi*180,'fansensorgeometry','line','fansensorspacing',sys_geo.p2v*sys_geo.sid/sys_geo.sdd,'outputsize',sys_geo.voxel_num);

end

