function [recon_astra] = astraIFDK(sino,sys_geo,angle_list)
%ASTRAIFDK Conebeam reconstruction (FDK) using ASTRA toolbox
%   INPUT:
%       sino: 3-D sinogram, shape: pixel (along axis u) * pixel (along axis v) * angle
%     sys_geo: a struct of the system's geometrical parameters, including:
%         pixel_num: number of the pixels of the detector
%         phantom_length: size of the phantom, scalar, unit: mm
%         voxel_num: number of the voxels of the FOV descretization
%         detector_length: length of the detector, unit: mm
%         d2v, v2d, p2d, p2v: the transformation coefficient from unit Distance (mm) to Voxel size, Voxel size to unit Distance, Pixel size to unit Distance, and Pixel size to Voxel size
%         sid0, sdd0, sid, sdd: the distance between Source and Iso-center, and between Source and Detector, unit: mm or voxel
%         pixel_num3, phantom_length3, voxel_num3, detector_length3, d2v3, v2d3, p2d3, p2v3: geometry parameters along axis z or v
%       angle_list: projection angles, 1D vector, unit: degree
%  Notice that the filter type cannot be selected due to limitation of ASTRA in 3D recon

%{
If you have trouble with ASTRA in linux platform reading "GLIBCXX_3.4.xx",
run the command below to start MATLAB with a lower version libstdc
```
LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
```
%}

proj_geom = astra_create_proj_geom('cone', sys_geo.p2v,sys_geo.p2v3,sys_geo.pixel_num3,sys_geo.pixel_num,angle_list/180*pi,sys_geo.sid,sys_geo.sdd-sys_geo.sid);
vol_geom = astra_create_vol_geom(sys_geo.voxel_num,sys_geo.voxel_num,sys_geo.voxel_num3); % parameter order: rows, columns, slices (y,x,z)

sino_id = astra_mex_data3d('create','-sino',proj_geom,permute(flipud(sino),[1 3 2])); % adapt to ASTRA geoemtry (the rotation angle is at the second axis)
recon_id = astra_mex_data3d('create','-vol',vol_geom,0);

cfg = astra_struct('FDK_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = recon_id;

alg_id = astra_mex_algorithm('create',cfg);
astra_mex_algorithm('run', alg_id);

recon_astra = astra_mex_data3d('get', recon_id);
recon_astra = rot90(recon_astra,-1);
end

