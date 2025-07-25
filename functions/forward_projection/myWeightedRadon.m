function [sino] = myWeightedRadon(phant,sys_geo,varargin)
%MYWEIGHTEDRADON Weighted Radon Transform - parallel beam
%   INPUT:
%     phant: matrix of the phantom, 2D
%     sys_geo: a struct of the system's geometrical parameters, including:
%         pixel_num: number of the pixels of the detector
%         phantom_length: size of the phantom, scalar, unit: mm
%         voxel_num: number of the voxels of the FOV descretization
%         detector_length: length of the detector, unit: mm
%         d2v, v2d, p2d, p2v: the transformation coefficient from unit Distance (mm) to Voxel size, Voxel size to unit Distance, Pixel size to unit Distance, and Pixel size to Voxel size
%         sid0, sdd0, sid, sdd: the distance between Source and Iso-center, and between Source and Detector, unit: mm or voxel
%         
%     varargin: 
%         theta_list: projection angles, 1D vector, unit: degree, default value: [0:179]
%         weight_function: the weight function for every voxel and angle, with an input of (x,y,theta), where x and y are in the unit of voxel, and theta in the unit of radian

% set varargin
defaults = {'theta_list',[0:179], 'weight_function',@(x,y,theta)(1)};
params = parseKeyValuePairs(defaults,varargin{:});
theta_list = params.theta_list;
weight_function = params.weight_function;
% disp(theta_list);
% disp(weight_function);

% forward projection
sino = zeros(sys_geo.pixel_num,length(theta_list));

% go through every angle and pixel (ray-driven)
for itheta = 1:length(theta_list)
    if mod(itheta-1,10)==0
        disp(['progress: ',num2str(itheta),'/',num2str(length(theta_list))]);
    end
    theta = (theta_list(itheta) )/180*pi; % transfer into radian

    parfor ipixel = 1 : sys_geo.pixel_num
        t = (ipixel - 0.5 - sys_geo.pixel_num /2)*sys_geo.p2d; % set the center of one pixel as the endpoint of the ray, then turn it into coordinates along u direction

        point = [0 nan nan]; % initialize a matrix to save the intersection point of the ray and voxels

        det = sys_geo.d2v*[(sys_geo.sdd0-sys_geo.sid0)*sin(theta)+t*cos(theta), -(sys_geo.sdd0-sys_geo.sid0)*cos(theta)+t*sin(theta)]; % the coordinate of the pixel of this rat, unit: voxel
        src = sys_geo.d2v*[t*cos(theta), t*sin(theta)]; % the coordinate of the source, unit: voxel. note: for parallel beam this is the coordinate of the perpendicular foot from iso-center to the ray


        % express the ray in the form of k*det+(1-k)*src, k\in(0,1), larger k means further from source
        for x0 = (-sys_geo.voxel_num/2:sys_geo.voxel_num/2) % calculate intersection point of vertical grids
            k = (x0 - src(1))/(det(1)-src(1));
            y = k*det(2)+(1-k)*src(2);
            if y>=-sys_geo.voxel_num/2&&y<=sys_geo.voxel_num/2
                point = [point; [k, x0, y]];
            end
        end
        for y0 = (-sys_geo.voxel_num/2:sys_geo.voxel_num/2) % calculate for horizontal grids
            k = (y0 - src(2))/(det(2)-src(2));
            x = k*det(1)+(1-k)*src(1);
            if x<=sys_geo.voxel_num/2&&x>=-sys_geo.voxel_num/2
                point = [point; [k, x, y0]];
            end
        end
        point = point(2:end,:);
        point = sortrows(point,1); % sort the intersection points according to k, such that every two adjacent points are in one voxel
        total_dist = norm(det-src); % the total length of this ray, unit: voxel


        for i = 1:(size(point,1)-1) % calculate the forward projection
            middle = (point(i,2:3)+point(i+1,2:3))/2; % the middle point of two intersection points, to determine the voxel where they locate
            dist = (point(i+1,1)-point(i,1))*total_dist; % the distance between the two points, unit: voxel
            if dist <= 1e-7 % skip repeated intersections 
                continue; 
            end
            coordinate_old = ceil(middle);
            coordinate = ceil(middle)+sys_geo.voxel_num/2; coordinate_new = coordinate; % the coordinate of voxels in x-y plane
            coordinate_new = [sys_geo.voxel_num + 1 - coordinate(2) , coordinate(1) ]; % the index of voxels in the phantom matrix

            sino(ipixel,itheta) = sino(ipixel,itheta) + dist * phant(coordinate_new(1),coordinate_new(2)) * weight_function(coordinate_old(1),coordinate_old(2),theta);

        end

    end

end

end



