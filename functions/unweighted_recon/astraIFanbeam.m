function [recon] = astraIFanbeam(sino,sys_geo,angle_list,varargin)
%ASTRAIFANBEAM Fanbeam reconstruction (FBP) using ASTRA toolbox
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
%         filter_type: type of the filter, 'ram-lak' (default), 'shepp-logan', 'cosine', 'hamming', 'hann', 'none'

%{
If you have trouble with ASTRA in linux platform reading "GLIBCXX_3.4.xx",
run the command below to start MATLAB with a lower version libstdc
```
LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
```
%}

% set varargin
defaults = {'filter_type','ram-lak'};
params = parseKeyValuePairs(defaults,varargin{:});
filter_type = params.filter_type;

proj_geom = astra_create_proj_geom('fanflat', sys_geo.p2v,sys_geo.pixel_num,angle_list/180*pi,sys_geo.sid,sys_geo.sdd-sys_geo.sid);
vol_geom = astra_create_vol_geom(sys_geo.voxel_num,sys_geo.voxel_num);

sino_id = astra_mex_data2d('create','-sino',proj_geom,(fliplr(sino))');
recon_id = astra_mex_data2d('create','-vol',vol_geom,0);

cfg = astra_struct('FBP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = recon_id;
cfg.option.FilterType = filter_type;

alg_id = astra_mex_algorithm('create',cfg);
astra_mex_algorithm('run', alg_id);

recon = astra_mex_data2d('get', recon_id);
end

