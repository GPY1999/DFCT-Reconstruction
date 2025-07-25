function [recon] = myWeightedIConebeamPlanar(sino,sys_geo,varargin)
%MYWEIGHTEDICONEBEAMPLANAR Weighted inverted Radon transform, suitable for cone-beam and planar detector
%       The inversion formula is based on Boman & Strömberg's and Huang's articles
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
% 
%     varargin: 
%         theta_list: projection angles, 1D vector, unit: degree, default value: [0:359]
%         weight_function: the weight function for every voxel and angle, with an input of (x,y,theta), where x and y are in the unit of voxel, and theta in the unit of radian
% UPDATES:
%   2025/07/25 (GUO Peiyuan): complete comments
%   2024/11/02 (GUO Peiyuan): first version

% set varargin
defaults = {'theta_list',[0:359], 'weight_function',@(x,y,theta)(1)};
params = parseKeyValuePairs(defaults,varargin{:});
beta_list = params.theta_list;
weight_function = params.weight_function;

% cosine weighting
s = ceil(-sys_geo.pixel_num/2)+1 : ceil(sys_geo.pixel_num/2); % handle odd pixel number
t = ceil(-sys_geo.pixel_num3/2)+1 : ceil(sys_geo.pixel_num3/2); % pixel coordinate along axis v
[t_mesh,s_mesh] = meshgrid(t,s); 
det_spacing3 = sys_geo.p2v3 / sys_geo.sdd * sys_geo.sid; % pixel shape of equivalent detector along axis v
det_spacing = sys_geo.p2v / sys_geo.sdd * sys_geo.sid; % pixel shape of equivalent detector along axis u
sino_cosweight = sino .* sqrt(sys_geo.sid.^2+((t_mesh-0.5)*det_spacing3).^2) .* sys_geo.sid ...
    ./(sys_geo.sid^2+((s_mesh-0.5)*det_spacing).^2+((t_mesh-0.5)*det_spacing3).^2); % refer to my draft

% Hilbert filtering
s_extended = -sys_geo.pixel_num+1 : sys_geo.pixel_num; 
h_of_s = zeros(length(s_extended),1); % kernel of convolution
h_of_s = 1/pi ./ (s_extended*det_spacing);
h_of_s(sys_geo.pixel_num) = 0;
h_of_s = reshape(h_of_s,[],1); % adjust shape to convolute with the sinogram

sino_conv = convn(sino_cosweight,h_of_s,"same") * det_spacing; % convolution

% back-projection (adjacent projection)
I1 = zeros(sys_geo.voxel_num,sys_geo.voxel_num,sys_geo.voxel_num3);I2 = I1;
x = -sys_geo.voxel_num/2+0.5 : sys_geo.voxel_num/2-0.5; % recon image voxel coordinate
y = flip(-sys_geo.voxel_num/2+0.5 : sys_geo.voxel_num/2-0.5); 
z = -sys_geo.voxel_num3/2+0.5 : sys_geo.voxel_num3/2-0.5; 
[x_mesh, y_mesh, z_mesh] = meshgrid(x,y,z);

for iangle = 1 : length(beta_list)
    if mod(iangle-1,10)==0
        disp(['progress: ',num2str(iangle),'/',num2str(length(beta_list))]);
    end
    % beta: the angle between source-isocenter and Y-axis
    beta = (beta_list(iangle) )/180*pi; % transfer into radian
    % U: each voxel's RATIO between its projection length on the
    % optical axis and the source-isocenter distance
    U = (sys_geo.sid + sin(beta)*x_mesh - cos(beta)*y_mesh) ./ sys_geo.sid;
    % V: each voxel's horizontal projection coordinate on (equivalent) detector plane
    V = (cos(beta) * x_mesh + sin(beta) * y_mesh)/det_spacing ./ U; % unit: pixel
    Vp = V + 0.5 + sys_geo.pixel_num/2; % turn it into [1,det_num]
    % W: each voxel's vertical projection coordinate on equivalent detector plane
    W = z_mesh / det_spacing3 ./ U; % unit: pixel3
    Wp = W + 0.5 + sys_geo.pixel_num3 / 2; % turn into [1,det_num3]

    weight_map = 1 ./ weight_function(x_mesh,y_mesh,beta); 

    % theta: the angle between source-voxel and Y-axis
    theta = beta + atan(V*det_spacing./sys_geo.sid);
    I1 = I1 + cos(theta) .* interp2((1:sys_geo.pixel_num3)',(1:sys_geo.pixel_num)',sino_conv(:,:,iangle),Wp,Vp) .* 1./U .* weight_map;
    I2 = I2 + sin(theta) .* interp2((1:sys_geo.pixel_num3)',(1:sys_geo.pixel_num)',sino_conv(:,:,iangle),Wp,Vp) .* 1./U .* weight_map; 
end

I1 = cat(2,I1(:,1,:),I1); I2 = -cat(1,I2(1,:,:),I2);

recon = diff(I1,1,2) + diff(I2,1,1); 
recon = recon * (beta_list(2) - beta_list(1)) / 180 * pi; % consider integration interval
recon = recon / 4 /pi;
end

