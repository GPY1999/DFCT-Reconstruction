function [recon] = myWeightedIRadon(sino,sys_geo,varargin)
%MYWEIGHTEDIRADON Weighted inverted Radon transform
%       The inversion formula is based on Boman & Strömberg's article
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
%     varargin: 其他
%         theta_list: projection angles, 1D vector, unit: degree, default value: [0:179]
%         weight_function: the weight function for every voxel and angle, with an input of (x,y,theta), where x and y are in the unit of voxel, and theta in the unit of radian
% UPDATES:
%   2025/07/25 (GUO Peiyuan): complete comments
%   2024/06/11 (GUO Peiyuan): first version


% set varargin
defaults = {'theta_list',[0:179], 'weight_function',@(x,y,theta)(1)};
params = parseKeyValuePairs(defaults,varargin{:});
theta_list = params.theta_list;
weight_function = params.weight_function;

% Hilbert filtering
s_extended = -sys_geo.pixel_num+1 : sys_geo.pixel_num;
h_of_s = zeros(length(s_extended),1); % kernel of convolution
h_of_s = 1/pi ./ (s_extended*sys_geo.p2v);
h_of_s(sys_geo.pixel_num) = 0;
h_of_s = reshape(h_of_s,[],1); % adjust shape to convolute with the sinogram

sino_conv = convn(sino,h_of_s,"same") * sys_geo.p2v; % convolution

% back-projection (adjacent projection)
I1 = zeros(sys_geo.voxel_num);I2 = I1;
x = -sys_geo.voxel_num/2+0.5 : sys_geo.voxel_num/2-0.5; % recon image voxel coordinate
y = flip(-sys_geo.voxel_num/2+0.5 : sys_geo.voxel_num/2-0.5); 
[x_mesh, y_mesh] = meshgrid(x,y);
% detector pixel coordinate
for iangle = 1 : length(theta_list)
    if mod(iangle-1,10)==0
        disp(['progress: ',num2str(iangle),'/',num2str(length(theta_list))]);
    end
    theta = (theta_list(iangle) )/180*pi; % transfer into radian

    % V represents each voxel's projection coordinate on detector plane
    V = (cos(theta) * x_mesh + sin(theta) * y_mesh)/sys_geo.p2v; % unit: pixel
    V = V + 0.5 + sys_geo.pixel_num/2; % turn it into [1,det_num]
    weight_map = 1 ./ weight_function(x_mesh,y_mesh,theta); 
    I1 = I1 + cos(theta) * interp1((1:sys_geo.pixel_num)',sino_conv(:,iangle),V) .* weight_map;
    I2 = I2 + sin(theta) * interp1((1:sys_geo.pixel_num)',sino_conv(:,iangle),V) .* weight_map; 
end


I1 = [I1(:,1) I1]; I2 = -[I2(1,:) ; I2];

recon = diff(I1,1,2) + diff(I2,1,1); 
recon = recon * (theta_list(2) - theta_list(1)) / 180 * pi; % consider integration interval
recon = recon / 4 /pi;

end
