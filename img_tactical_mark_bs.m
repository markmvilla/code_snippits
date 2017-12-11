img_tactical%**************************************************************************
%
%   Author  Mark Villa
%   Date    Nov. 13, 2014
%
%   Pupose  Provide a class (functions & data) to generate a tactical
%           image from a given position and look angle.
%
%   Constructor img_tactical(geo_pt, pitch_rad, roll_rad, az_rad, camera)
%
%   Functions	function ndx  = get_pixel_point(obj, geopt)
%               function calculate_camera_parameters(
%                   geo_pt, pitch_rad, roll_rad, az_rad, camera)
%               function img_ndx = generate_perspective_image_ndx(
%                   geo_pt, pitch_rad, roll_rad, az_rad, camera)
%               function img_rgb = generate_perspective_image_rgb(
%                   geo_pt, pitch_rad, roll_rad, az_rad, camera)
%
%   Data    sky_index       Image color value of the sky.
%           default_index   Image default color value.
%           sky_rgb         Image RGB color value of the sky.
%           default_rgb     Image default RGB color value.
%           is_available    TRUE if image available, FALSE if not available
%           azimuth         + Right (deg)
%           elevation       + Down (deg)
%           focal_length    Focal length (m)
%           fov             Field of view (deg)
%           num_rows        Number of image rows. (pixels)
%           num_cols        Number of image cols. (pixels)
%           row_cen         Center pixel along rows. (pixels)
%           col_cen         Center pixel along columns. (pixels)
%           dy_m            Delta Y-axis spacing (m)
%           dx_m            Delta X-axis spacing (m)
%           yvec            Vector containing the Y-values image (m, along rows).
%           xvec            Vector containing the X-values image (m, along cols).
%
%**************************************************************************

classdef img_tactical
    properties
        %------------------------------------------------------------------
        % Class Data
        %------------------------------------------------------------------
        sky_index           % Image color value of the sky.
        default_index       % Image default color value.
        sky_rgb             % Image RGB color value of the sky.
        default_rgb         % Image default RGB color value.
        is_available        % TRUE if image available, FALSE if not available
        geo_0               % Geodetic point of camera projection (LLA)
        ecf_0               % ECF point of camera projection (m)
        focal_length        % Focal length (m)
        fov                 % Field of view (deg)
        C_cam_to_ecf        % Camera-to-ECF transformation [3x3, orthonormal]
        num_rows            % Number of image rows. (pixels)
        num_cols            % Number of image cols. (pixels)
        row_cen             % Center pixel along rows. (pixels)
        col_cen             % Center pixel along columns. (pixels)
        dy_m                % Delta Y-axis spacing (m)
        dx_m                % Delta X-axis spacing (m)
        yvec                % Vector containing the Y-values image (m, along rows).
        xvec                % Vector containing the X-values image (m, along cols).
    end

    methods (Access=public)

        %------------------------------------------------------------------
        % Class Constructor :  Create a camera model.
        %
        % Input geopt       Point of projection [Lat(rad),Lon(rad),hae(m)].
        %       camera      Structure of camera parameters.
        %       q           Vector of body angles [roll,pitch,az] (rad).
        %       cmap        Map of image colors.
        %------------------------------------------------------------------

        function obj = img_tactical(...
                geo_pt, pitch_rad, roll_rad, az_rad, camera)

            obj.default_index   = 1;
            obj.sky_index       = 256;
            obj.default_rgb     = [0 0 0];
            obj.sky_rgb         = [150 150 256];

            if (nargin == 0)
                % Default Constrcctor
                obj.is_available = false;
                return
            elseif (nargin ~= 5)
                % Error statement
                disp('Error: Number of arguments not equal to 0 or 5 in img_tactical()')
                obj.is_available = false;
                return
            end

            % Set the camera parameters
            obj = calculate_camera_parameters(obj, ...
                geo_pt, pitch_rad, roll_rad, az_rad, camera);

            % Set the availabilty flag
            obj.is_available = true;
        end

        %------------------------------------------------------------------
        % Pupose    Calculate the camera parameters.
        %
        % Input     geo_pt      Geodetic point of light ray intersection.
        %                       [rad,rad,hae(m)]
        %           roll_rad    (+ right-wing-down).
        %           pitch_rad   (+ nose-up).
        %           az_rad      (+ cw or azimuth-right-turn).
        %           camera      Camera parameter structure.
        %------------------------------------------------------------------

        function obj = calculate_camera_parameters(obj, ...
                geo_pt, pitch_rad, roll_rad, az_rad, camera)

            global environ

            % Set the camera parameters
            obj.num_rows        = camera.num_rows;
            obj.num_cols        = camera.num_cols;
            obj.focal_length    = camera.f;
            obj.fov             = camera.fov;

            sfov2 = sin(obj.fov/2);
            obj.dy_m = 2 * obj.focal_length * sfov2 / (obj.num_rows-1);
            obj.dx_m = 2 * obj.focal_length * sfov2 / (obj.num_cols-1);

            % Set the pixel reference values (m)
            % - obj.img(1,1) =
            % - obj.yvec(1)  =
            % - obj.xvec(1)  =

            % Number of image rows and cols is ODD and
            % the image center is at (N+1)/2.
            obj.yvec = obj.dy_m * (-(obj.num_rows-1)/2:(obj.num_rows-1)/2);
            obj.xvec = obj.dx_m * (-(obj.num_cols-1)/2:(obj.num_cols-1)/2);
            obj.row_cen = (obj.num_rows + 1) / 2;
            obj.col_cen = (obj.num_cols + 1) / 2;

            % Number of image rows and cols is EVEN and
            % the image center is at (N/2)+1.
            % obj.yvec = obj.dy_m * (-(obj.num_rows-1)/2:(obj.num_rows-1)/2);
            % obj.xvec = obj.dx_m * (-(obj.num_cols-1)/2:(obj.num_cols-1)/2);
            % obj.row_cen = obj.num_rows / 2 + 1;
            % obj.col_cen = obj.num_cols / 2 + 1;

            % Set the positional information
            obj.geo_0 = geo_pt;
            obj.ecf_0 = geo_to_ecf(obj.geo_0, environ.ellip, environ.ustr);

            % Camera-to-Body transformation
            C_cam_to_body = get_dcm_cam_to_body(camera.el, camera.roll, camera.az);

            % Body-to-NWU transformation
            C_body_to_nwu = get_dcm_body_to_nwu(pitch_rad, roll_rad, az_rad);

            % NWU-to-ENU transformation
            C_nwu_to_enu = get_dcm_nwu_to_enu();

            % ENU-to-ECF transformation
            C_enu_to_ecf = get_dcm_enu_to_ecf(obj.geo_0.lat, obj.geo_0.lon);

            % Camera-to-ECF transformation
            obj.C_cam_to_ecf ...
                = C_enu_to_ecf * C_nwu_to_enu * C_body_to_nwu * C_cam_to_body;

        end

        %------------------------------------------------------------------
        % Pupose    Calculate a geodetic position intersection by ray tracing
        %           from a projection point, through a camera pixel location
        %           to an interpolated terrain surface. The algorithm
        %           iterates until a predetermined accuracy is achieved.
        %
        % Input     geo         Geodetic point of light ray intersection.
        %                       [rad,rad,hae(m)]
        %           roll_rad    (+ right-wing-down).
        %           pitch_rad   (+ nose-up).
        %           az_rad      (+ cw or azimuth-right-turn).
        %           camera      Camera parameter structure.
        %
        % Output    img         Indexed image array containing data
        %------------------------------------------------------------------

        function img_ndx = generate_perspective_image(obj, ...
                geo_pt, pitch_rad, roll_rad, az_rad, camera, traj_ndx)

            global environ
            global ned_ter
            global img_ter
            global options
            global last_hit

            gcf;
            hold on;
            %scatter3(rad_to_deg(geo_pt.lon),rad_to_deg(geo_pt.lat),geo_pt.alt,'k*');

            % Set the camera parameters
            obj = calculate_camera_parameters(obj, ...
                geo_pt, pitch_rad, roll_rad, az_rad, camera);

            % Pre-allocate image memory and give it a default value
            img_ndx = uint8(obj.default_index * ones(obj.num_rows,obj.num_cols));
            if (traj_ndx == 1),
                last_hit.target = 5000; % Starting range
                disp('last_hit = 5000!')
            end

            for col = 1:obj.num_cols
                disp([num2str(col,'%3d') ' : ' num2str(obj.num_cols,'%3d')])
                for row = 1:obj.num_rows
                    disp(['r,c =',num2str(col,'%2d'),num2str(row,'%2d')])
                    % Initialize iteration parameters
                    alternate = -1;          %initializing variables
                    range = last_hit.target;       % sets it to last terrain hit
                    delta  = 0;           % Iteration step in (m)
                    increment = 1500;
                    fine_delta = increment;
                    thresh = 0.1;           % Altitude accuracy threshold
                    count  = 0;             % Iteration step counter
                    iter_finished = false;
                    found = false;
                    backtracking = true;

                    % Calculate the unit vector of the image point
                    cam = [obj.xvec(col); obj.yvec(row); obj.focal_length];
                    unit_cam = cam / norm(cam);

                    while (~iter_finished)

                        % Calculate ECF position of the end of the ray
                        % trace. Scale the unit vector to the full range
                        % and then transform the unit vector to ECF coord.
                        ecf_k = obj.ecf_0 + obj.C_cam_to_ecf * (range * unit_cam);

                        % Calculate Geodetic position of the end of the ray trace
                        geo_k = ecf_to_geo_olson(ecf_k, environ.ellip, environ.ustr);

                        %scatter3(rad_to_deg(geo_k.lon),rad_to_deg(geo_k.lat),geo_k.alt,'r.');

                        % Get the height of the terrain at the end of the ray trace
                        ter_ht = ned_ter.get_elevation_point(geo_k);

                        % Calculate the height of the ray trace above the terrain
                        diff = geo_k.alt - ter_ht;

                        % Check the vertical distance from the ray trace
                        % endpoint to the interpolated terrain surface. If the
                        % accuracy is not sufficient, then adjust the estimate
                        % of range and try again.
                        if (abs(diff) > thresh),
                            iter_finished = false;
                            if(((count>0) && (found == true))||(diff<0)||(found == true)),
                                found = true;
                                if(((diff<0)&&(count==0)&&(backtracking == true))||((alternate<0)&&(diff<0)&&(backtracking == true))||((backtracking == true)&&(diff<0))),
%                                     disp('hit spot-backing')
                                    backtracking = true;
                                    range = range - 500;
%                                     disp(['range =',num2str(range,'%2d')])
                                else
                                    backtracking = false;
%                                     disp('hit spot-binomial')
                                    fine_delta = fine_delta/2;
                                    if (diff > 0),
                                        range = range + fine_delta;
%                                         disp(['range =',num2str(range,'%2d')])
                                    else
                                        range = range - fine_delta;
%                                         disp(['range =',num2str(range,'%2d')])
                                    end
                                end
                            else
                                disp('missing')
                                found = false;
                                count = count + 1;
%                                 disp(['count',num2str(count,'%2d')])
                                alternate = -alternate;
                                delta = delta + increment;
                                range = range + alternate*delta;
                                disp(['range =',num2str(range,'%2d')])
                            end
                        else
                            iter_finished = true;
                            ndx = img_ter.get_pixel_nn_index(...
                            rad_to_deg(geo_k.lat), rad_to_deg(geo_k.lon));
                            img_ndx(row,col) = ndx;
                            disp(['*** Final range = ', num2str(range,'%.6f') ...
                            ', Cnt=', num2str(count,'%2d')])
                        end
                        if (range > 20000),
                            img_ndx(row,col) = obj.sky_index;
                            iter_finished = true;
                            % disp(['*** Too far!, Final range = ', num2str(range,'%.6f') ...
                            %  ', Cnt=', num2str(count,'%2d')])
                        end
                    end
                    if (range < 20000),
                        last_hit.target = range;
                    end
                    if (options.produce_camera_rays),
                        if     ((row==1            && col==1) || ...
                                (row==obj.num_rows && col==1) || ...
                                (row==1            && col==obj.num_cols) || ...
                                (row==obj.num_rows && col==obj.num_cols)),
                            plot3([rad_to_deg(geo_pt.lon) rad_to_deg(geo_k.lon)], ...
                                [rad_to_deg(geo_pt.lat) rad_to_deg(geo_k.lat)], ...
                                [geo_pt.alt geo_k.alt], 'm-')
                        end
                        if (row==obj.row_cen && col==obj.col_cen),
                            plot3([rad_to_deg(geo_pt.lon) rad_to_deg(geo_k.lon)], ...
                                [rad_to_deg(geo_pt.lat) rad_to_deg(geo_k.lat)], ...
                                [geo_pt.alt geo_k.alt], 'k-')
                        end
                    end
                end
            end
        end
    end
end

%**************************************************************************
