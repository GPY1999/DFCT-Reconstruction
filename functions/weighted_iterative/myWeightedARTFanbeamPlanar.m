function [recon] = myWeightedARTFanbeamPlanar(sino,sys_geo,varargin)
%MYWEIGHTEDARTFANBEAMPLANAR_MOD2 Iterative recon: Algebraic reconstruction techinique (ART), for fan-beam and planar detector
%   Weighted recon also supported
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
%         weight_function: the weight function for every voxel and angle, with an input of (x,y,theta), where x and y are in the unit of voxel, and theta in the unit of radian
%         nART: number of iteration, default value: 1
%         recon_init: initial value of iteration, default value: zeros_like(recon)
% UPDATES:
%   2025/07/25 (GUO Peiyuan): complete comments
%   2025/05/29 (GUO Peiyuan): second version, random backprojection
%   2025/05/14 (GUO Peiyuan): first version

% set varargin
defaults = {'theta_list',[0:359], 'weight_function',@(x,y,theta)(1), 'nART',1, 'recon_init',zeros(sys_geo.voxel_num)};
params = parseKeyValuePairs(defaults,varargin{:});
beta_list = params.theta_list;
weight_function = params.weight_function;
nART = params.nART;

recon = reshape(params.recon_init,[],1);
sino_reshape = reshape(sino,[],1);

[p,q] = meshgrid(1:length(beta_list), 1:sys_geo.pixel_num); % to generate ray order for back projection

% iterative recon
for nart = 1 : nART
    disp(strcat('nart=',num2str(nart)));
    rayOrderList = [p(:) q(:)];
    rayOrderList = rayOrderList(randperm(numel(p)),:); % random order
    for iray = 1 : size(rayOrderList,1)
        iangle = rayOrderList(iray,1);
        beta = beta_list(iangle) / 180 * pi; % transfer into radian
        if mod(iray,size(rayOrderList,1)/20)==0 
            disp(strcat('progress in this iteration:',32,num2str(round(iray/size(rayOrderList,1)*100)),'%'))
        end
        ipixel = rayOrderList(iray,2);
        [point_index,~,point_weight] = get_sys_mat(beta,ipixel,sys_geo,weight_function);

        if ~isempty(point_weight)
            correct = (sum(recon(point_index).*point_weight) ...
                - sino_reshape((iangle-1)*sys_geo.pixel_num+ipixel)) ./ sum((point_weight).^2);
            recon(point_index) = recon(point_index) - correct.*point_weight;
        end
    end
    recon(recon<0) = 0;
    figure,imshow(reshape(recon,sys_geo.voxel_num,sys_geo.voxel_num),[]); title(num2str(nart));
end
recon = reshape(recon,sys_geo.voxel_num,sys_geo.voxel_num);
end

function [point_index, point_length, point_weight] = get_sys_mat(theta,ipixel,sys_geo,weight_function)
% GET_SYS_MAT Get one row of the system matrix, in the sparse form, which means the intersected voxel index and length of one ray
% The phantom is reshaped into 1-D vectorl
%   INPUT:
%     theta, ipixel: current source rotation angle (unit: radian) and pixel coordinate (unit: pixel) of this ray
%       Specifically, theta rotates counter clockwise from positive y axis, and ipixel increases from 1 counter clockwise
%       Notice that theta is not the angle of the ray, but the angle of source-isocenter
%   OUTPUT:
%     point_length: the ray length intersected with every voxel
%     point_weight: ray length * weight function

t = (ipixel - 0.5 - sys_geo.pixel_num /2)*sys_geo.p2d; % set the center of one pixel as the endpoint of the ray, then turn it into coordinates along u direction


point = [0 nan nan]; % initialize a matrix to save the intersection point of the ray and voxels

det = sys_geo.d2v*[(sys_geo.sdd0-sys_geo.sid0)*sin(theta)+t*cos(theta), -(sys_geo.sdd0-sys_geo.sid0)*cos(theta)+t*sin(theta)]; % the coordinate of the pixel of this rat, unit: voxel
src = sys_geo.d2v*sys_geo.sid0*[sin(-theta), cos(-theta)]; % the coordinate of the source, unit: voxel


% express the ray in the form of k*det+(1-k)*src, k\in(0,1), larger k means further from source
% calculate intersection point of vertical grids
x0 = (-sys_geo.voxel_num/2:sys_geo.voxel_num/2)';
k_x = (x0 - src(1))/(det(1)-src(1));
y_x = k_x*det(2)+(1-k_x)*src(2);
x0(y_x<=-sys_geo.voxel_num/2 | y_x>=sys_geo.voxel_num/2) = [];
k_x(y_x<=-sys_geo.voxel_num/2 | y_x>=sys_geo.voxel_num/2) = [];
y_x(y_x<=-sys_geo.voxel_num/2 | y_x>=sys_geo.voxel_num/2) = [];

y0 = (-sys_geo.voxel_num/2:sys_geo.voxel_num/2)';
k_y = (y0 - src(2))/(det(2)-src(2));
x_y = k_y*det(1)+(1-k_y)*src(1);
y0(x_y<=-sys_geo.voxel_num/2 | x_y>=sys_geo.voxel_num/2) = [];
k_y(x_y<=-sys_geo.voxel_num/2 | x_y>=sys_geo.voxel_num/2) = [];
x_y(x_y<=-sys_geo.voxel_num/2 | x_y>=sys_geo.voxel_num/2) = [];
point = [[k_x, x0, y_x];[k_y, x_y, y0]];

point = point(2:end,:);
point = sortrows(point,1); % sort the intersection points according to k, such that every two adjacent points are in one voxel
total_dist = norm(det-src); % the total length of this ray, unit: voxel

point_index = zeros(size(point,1)-1 ,1);
point_length = zeros(size(point,1)-1 ,1);
point_weight = ones(size(point,1)-1 ,1);

for i = 1:(size(point,1)-1) % calculate the forward projection
    middle = (point(i,2:3)+point(i+1,2:3))/2; % the middle point of two intersection points, to determine the voxel where they locate
    dist = (point(i+1,1)-point(i,1))*total_dist; % the distance, unit: voxel
    if dist <= 1e-7 % skip repeated intersections 
        continue;
    end
    coordinate_old = ceil(middle);
    coordinate = ceil(middle)+sys_geo.voxel_num/2; % the coordinate of voxels in x-y plane
    coordinate_new = [sys_geo.voxel_num + 1 - coordinate(2) , coordinate(1) ]; % the index of voxels in the phantom matrix

    point_index(i) = (coordinate_new(2)-1)*sys_geo.voxel_num + coordinate_new(1);
    point_length(i) = dist;
    point_weight(i) = dist * weight_function(coordinate_old(1),coordinate_old(2),theta);
end
point_length(point_index==0)=[];
point_weight(point_index==0)=[];
point_index(point_index==0)=[]; % skip repeated intersections 
end