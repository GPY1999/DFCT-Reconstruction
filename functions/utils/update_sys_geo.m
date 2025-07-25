function [sys_geo_new] = update_sys_geo(sys_geo_old)
%UPDATE_SYS_GEO Update system geometry. 3D supported
%   Input a struct of system geometry, and complete the absent attributes

sys_geo_new = sys_geo_old;
sys_geo_new.v2d = sys_geo_new.phantom_length / sys_geo_new.voxel_num; % transfer factor from voxel size to unit distance
test_consistency(sys_geo_old,'v2d',sys_geo_new.v2d);
sys_geo_new.d2v = 1 / sys_geo_new.v2d; % transfer factor from unit distance to voxel size
test_consistency(sys_geo_old,'d2v',sys_geo_new.d2v);
sys_geo_new.p2d = sys_geo_new.detector_length / sys_geo_new.pixel_num; % transfer factor from pixel size to unit distance
test_consistency(sys_geo_old,'p2d',sys_geo_new.p2d);
sys_geo_new.p2v = sys_geo_new.p2d * sys_geo_new.d2v; % transfer factor from pixel size to voxel size
test_consistency(sys_geo_old,'p2v',sys_geo_new.p2v);
sys_geo_new.sdd = sys_geo_new.sdd0 * sys_geo_new.d2v; % distance between source and detector
test_consistency(sys_geo_old,'sdd',sys_geo_new.sdd);
sys_geo_new.sid = sys_geo_new.sid0 * sys_geo_new.d2v; % distance between source and iso-center
test_consistency(sys_geo_old,'sid',sys_geo_new.sid);

sys_geo_new.idd = sys_geo_new.sdd - sys_geo_new.sid; % distance between iso-center and detector

% parameters for the third dimension
if isfield(sys_geo_new,'detector_length3') && isfield(sys_geo_new,'pixel_num3') && isfield(sys_geo_new,'phantom_length3') && isfield(sys_geo_new,'voxel_num3')
    sys_geo_new.p2d3 = sys_geo_new.detector_length3 / sys_geo_new.pixel_num3; % pixel size along v direction
    sys_geo_new.v2d3 = sys_geo_new.phantom_length3 / sys_geo_new.voxel_num3; % voxel size along z direction
    sys_geo_new.d2v3 = 1 / sys_geo_new.v2d3;
    sys_geo_new.p2v3 = sys_geo_new.p2d3 * sys_geo_new.d2v3;
end

% pitch for helical scan
if isfield(sys_geo_new,'pitch0')
    sys_geo_new.pitch = sys_geo_new.pitch0 * sys_geo_new.d2v; % pitch
    sys_geo_new.z_start = sys_geo_new.z_start0 * sys_geo_new.d2v; % the starting point of source in z direction
end
end

function test_consistency(sys_geo,var_name,var_value)
% check if the value of an attribute consists with a given value给定的值
% var_name: string
if isfield(sys_geo,var_name)
    temp = eval(strcat("sys_geo.",var_name,";"));
    if temp ~= var_value
        error(strcat("data inconsistent for ",var_name));
    end
end
end
