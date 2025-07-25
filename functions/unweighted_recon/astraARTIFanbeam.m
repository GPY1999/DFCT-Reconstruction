function [recon] = astraARTIFanbeam(sino,sys_geo,angle_list,varargin)
%ASTRAARTIFANBEAM Iterative Fanbeam reconstruction (ART, SART, SIRT...) using ASTRA toolbox
%   INPUT:
%     sino: 2-D sinogram, shape: pixel * angle. Notice that the sinogram acquired by MATLAB need to be transposed, and the real data needs flip left-right and then transposed
%     sys_geo: a struct of the system's geometrical parameters, including:
%         pixel_num: number of the pixels of the detector
%         phantom_length: size of the phantom, scalar, unit: mm
%         voxel_num: number of the voxels of the FOV descretization
%         detector_length: length of the detector, unit: mm
%         d2v, v2d, p2d, p2v: the transformation coefficient from unit Distance (mm) to Voxel size, Voxel size to unit Distance, Pixel size to unit Distance, and Pixel size to Voxel size
%         sid0, sdd0, sid, sdd: the distance between Source and Iso-center, and between Source and Detector, unit: mm or voxel
%     angle_list: projection angles, 1D vector, unit: degree
%         
%     varargin:
%         algorithm: type of iterative algorithm, 'ART' (default), 'SART', 'SIRT'
%         recon_init: initial value of iteration, default value: zeros_like(recon)
%         iter_num: number of iteration for every ray, default value: 10

%{
If you have trouble with ASTRA in linux platform reading "GLIBCXX_3.4.xx",
run the command below to start MATLAB with a lower version libstdc
```
LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
```
%}

% set varargin
defaults = {'algorithm',"ART",'recon_init',[],'iter_num',10};
params = parseKeyValuePairs(defaults,varargin{:});
algorithm = upper(params.algorithm);
iter_num = double(int8(params.iter_num));
recon_init = params.recon_init;
if isempty(recon_init)
    recon_init = 0;
end

proj_geom = astra_create_proj_geom('fanflat', sys_geo.p2v,sys_geo.pixel_num,angle_list/180*pi,sys_geo.sid,sys_geo.sdd-sys_geo.sid);
vol_geom = astra_create_vol_geom(sys_geo.voxel_num,sys_geo.voxel_num);
proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);


sino_id = astra_mex_data2d('create','-sino',proj_geom,(fliplr(sino))');
recon_id = astra_mex_data2d('create','-vol',vol_geom,recon_init);
% sinomask_id = astra_mex_data2d('create','-sino',proj_geom,(fliplr(sino_mask))');

switch algorithm
case "ART"
    cfg = astra_struct('ART');
case "SART"
    cfg = astra_struct('SART');
case "SIRT"
    cfg = astra_struct('SIRT');
case "SART_CUDA"
    cfg = astra_struct('SART_CUDA');
case "SIRT_CUDA"
    cfg = astra_struct('SIRT_CUDA');
otherwise
    error('Unknown reconstruction algorithm. Only ART, SART, and SIRT supported.')
end

cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = recon_id;

cfg.option.MinConstraint = 0;
% cfg.option.SinogramMaskId = sinomask_id;

[p,q] = meshgrid(0:length(angle_list)-1, 0:sys_geo.pixel_num-1);
rayOrderList = [p(:) q(:)];
rayOrderList = rayOrderList(randperm(numel(p)),:);
cfg.option.RayOrder = 'custom';
cfg.option.RayOrderList = rayOrderList;

art_id = astra_mex_algorithm('create',cfg);
astra_mex_algorithm('iterate', art_id, iter_num * sys_geo.pixel_num * length(angle_list));

recon = astra_mex_data2d('get', recon_id);

end