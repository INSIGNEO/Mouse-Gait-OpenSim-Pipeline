clc; clear; close all;

%% ==============================
% CONFIGURATION
% ==============================
% Mouse number
mouse_age = "22weeks";
mouse_name = "Cage5none";

% Anatomical scaling (scale the bones and muscle path points)
anatomicalScaling = true;

unit_conversion = 0.1; % CT(mm) to OpenSim(cm)

% Path to your template OpenSim model
osim_input_file = "../Model_mouse_right.osim";

output_folder = strcat(mouse_name,'/Model_data/');
output_file = fullfile(output_folder, ...
    sprintf('Model_mouse_right_markers.osim'));

% ---- MANUAL PRE-ROTATION (DEGREES) ----
% Adjust these if the CT bone is "upside down" or "sideways" vs OpenSim
% [Roll(X), Pitch(Y), Yaw(Z)]
rot_adjust.pelvis = [90, 90, 90];
rot_adjust.thigh_r  = [90, 0, 0];  % e.g., [90, 0, 0] if rotated 90 deg on X
rot_adjust.leg_r  = [0, 90, 180];
rot_adjust.foot_r   = [0, 90, 180]; % [0, 180, 180]

full_path = strcat(mouse_age,'/',mouse_name,'/CT_data');
cd(full_path)
% File definitions
bones.pelvis.ct = "CT_Right_Pelvis.stl"; bones.pelvis.osim = "../../../STL_opensim/Original_right_pelvis.stl";
bones.thigh_r.ct  = "CT_Right_Femur.stl";  bones.thigh_r.osim  = "../../../STL_opensim/femur_r.stl";
bones.leg_r.ct  = "CT_Right_Tibia.stl";  bones.leg_r.osim  = "../../../STL_opensim/tibfib_r.stl";
bones.foot_r.ct   = "CT_Right_Foot.stl";   bones.foot_r.osim   = "../../../STL_opensim/full_foot_r.stl";

markers.Right_Illiac.file = "CT_Marker_RASIS.stl";       markers.Right_Illiac.bone = "pelvis";
markers.Right_Hip.file   = "CT_Marker_Right_Hip.stl";   markers.Right_Hip.bone   = "pelvis";
markers.Left_Illiac.file = "CT_Marker_LASIS.stl";       markers.Left_Illiac.bone = "pelvis";
markers.Right_Knee.file  = "CT_Marker_Right_Knee.stl";  markers.Right_Knee.bone  = "thigh_r";
markers.Right_Ankle.file = "CT_Marker_Right_Ankle.stl"; markers.Right_Ankle.bone = "leg_r";
markers.Right_Met.file   = "CT_Marker_Right_Foot.stl";   markers.Right_Met.bone   = "foot_r";

bone_names = fieldnames(bones);
marker_names = fieldnames(markers);
final_tforms = struct();
bone_scales = struct();

%% ---- REGISTER BONES ----
figure('Color', 'w', 'Name', 'OpenSim (Blue) vs Registered CT (Red)');

for i = 1:length(bone_names)
    name = bone_names{i};
    fprintf("Registering %s...\n", name);
    
    % 1. Load and scale
    ct_tri = stlread(bones.(name).ct);
    os_tri = stlread(bones.(name).osim);
    ct_pts = ct_tri.Points * unit_conversion;
    os_pts = os_tri.Points;
    
    % 2. APPLY MANUAL PRE-ROTATION
    % Convert degrees to radians and create rotation matrix
    deg = rot_adjust.(name);
    R_manual = deg2rad(deg);
    RotM = eul2rotm(R_manual, 'XYZ'); 
    
    % Rotate CT points around their own center
    mu_ct_raw = mean(ct_pts, 1);
    ct_pts_centered = ct_pts - mu_ct_raw;
    ct_pts_rotated = (RotM * ct_pts_centered')' + mu_ct_raw;
    
    % 3. Centroid Pre-Alignment (Translation)
    mu_ct_rot = mean(ct_pts_rotated, 1);
    mu_os = mean(os_pts, 1);
    t_init = mu_os - mu_ct_rot;
    ct_init = ct_pts_rotated + t_init;
    
    % 4. ICP Registration
    source_pc = pointCloud(ct_init);
    target_pc = pointCloud(os_pts);
    [tform_icp, ~, rmse] = pcregistericp(source_pc, target_pc, ...
        'Metric', 'pointToPoint', 'MaxIterations', 500, 'Tolerance', [0.00000001,0.00000005], 'InlierRatio', 0.9);
    
    % 5. Build Combined 4x4 Matrix
    % This accounts for: [Raw CT] -> [Rotate] -> [Translate] -> [ICP]
    T_rot = eye(4); T_rot(1:3,1:3) = RotM;
    % Adjustment for rotating around center
    T_cent = eye(4); T_cent(1:3,4) = -mu_ct_raw';
    T_recenter = eye(4); T_recenter(1:3,4) = mu_ct_raw';
    
    T_trans = eye(4); T_trans(1:3, 4) = t_init';
    T_icp = tform_icp.T'; 
    
    % Order: Recentered * Translation * ICP * Rotation * Centering
    final_tforms.(name) = T_icp * T_trans * T_recenter * T_rot * T_cent;
    
    % --- PLOTTING ---
    hold on;
    patch('Faces', os_tri.ConnectivityList, 'Vertices', os_tri.Points, ...
        'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    ct_v_homog = [ct_pts, ones(size(ct_pts,1),1)] * final_tforms.(name)';
    patch('Faces', ct_tri.ConnectivityList, 'Vertices', ct_v_homog(:,1:3), ...
        'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % title(sprintf("%s (RMSE: %.4f)", name, rmse));
    axis equal; grid off; view(3); camlight; lighting gouraud;

    % Get dimensions (Bounding Box)
    ct_dim = max(ct_v_homog(:,1:3)) - min(ct_v_homog(:,1:3));
    os_dim = max(os_tri.Points) - min(os_tri.Points);
    
    % Scale Factor = CT / OpenSim
    % Using a safety check to avoid division by zero
    scale_xyz = ct_dim ./ os_dim;
    bone_scales.(name) = scale_xyz;
    
    fprintf("%s Scale Factors: [X:%.4f, Y:%.4f, Z:%.4f]\n", name, scale_xyz(1), scale_xyz(2), scale_xyz(3));

end

%% ---- TRANSFORM & PLOT MARKERS ----
results = zeros(length(marker_names), 3);
for i = 1:length(marker_names)
    m_name = marker_names{i};
    b_name = markers.(m_name).bone;
    
    m_tri = stlread(markers.(m_name).file);
    m_center_ct = mean(m_tri.Points * unit_conversion, 1);
    
    T = final_tforms.(b_name);
    m_pos_homog = [m_center_ct, 1] * T';
    results(i, :) = m_pos_homog(1:3);
    
    hold on
    plot3(results(i,1), results(i,2), results(i,3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [1.0, 0.41, 0.70]);
    % text(results(i,1), results(i,2), results(i,3), m_name, 'FontSize', 20);
end

%% ---- UPDATE OSIM FILE (MARKERS + SCALING) ----
% Load the XML document
cd('../../')
dom = xmlread(osim_input_file);
all_markers = dom.getElementsByTagName('Marker');

%% 1. Update Markers 
for k = 0:all_markers.getLength - 1
    this_marker = all_markers.item(k);
    
    % Get the name of the marker from the XML attribute
    xml_marker_name = char(this_marker.getAttribute('name'));
    
    % Check if this marker matches any in our 'results' table
    match_idx = find(strcmp(marker_names, xml_marker_name));
    
    if ~isempty(match_idx)
        % Get the new coordinates
        new_coords = results(match_idx, :)./100;
        new_coords_str = sprintf('%.8f %.8f %.8f', new_coords(1), new_coords(2), new_coords(3));
        
        % Find the <location> tag within this marker
        location_node = this_marker.getElementsByTagName('location').item(0);
        
        if ~isempty(location_node)
            % Update the text inside the <location> tag
            location_node.getFirstChild.setData(new_coords_str);
            fprintf("Updated %s in XML to: [%s]\n", xml_marker_name, new_coords_str);
        end
    end
end

%% 2 Update bones (Scaling)
if anatomicalScaling == true
    % File and folder names
    vtpNames = {"../STL_opensim/Pelvis_r.vtp",
                "../STL_opensim/femur_r.vtp",
                "../STL_opensim/patella_r.vtp",
                "../STL_opensim/tibfib_r.vtp",
                "../STL_opensim/foot_r.vtp",
                "../STL_opensim/phalanges_r.vtp"};
    
    vtpScaleFactors = [bone_scales.pelvis;
                       bone_scales.thigh_r;
                       bone_scales.leg_r;
                       bone_scales.leg_r;
                       bone_scales.foot_r;
                       bone_scales.foot_r];
    
    vtpOutputNames = {'pelvis','femur','patella','tibia','foot','toes'};
    
    output_vtp = strcat(output_folder,'Geometry/');
    scaled_mesh_folder = output_vtp;
    
    scaledMeshMatrix = cell(length(vtpNames),1);
    
    for ivtp = 1:length(vtpNames)
        % Output VTP file
        scaled_mesh_file = fullfile(scaled_mesh_folder,strcat(vtpOutputNames{ivtp},'_scaled.vtp'));
        % Read VTP using VTK
        %% 1. Read the original VTP mesh
        % vtkRead returns a MATLAB structure.
        mesh = vtkRead(char(vtpNames{ivtp}));
        
        % Check that the mesh contains points
        if ~isfield(mesh, 'points') || isempty(mesh.points)
            error('No mesh points found in: %s', vtpNames{ivtp});
        end
        
        vertices = mesh.points;
        
        %% 2. Scale vertices around their centroid
        centroid = mean(vertices, 1);
        
        vertices_scaled = ...
            (vertices - centroid) .* reshape(vtpScaleFactors(ivtp,:), 1, 3) + centroid;
        
        % Replace only the vertex coordinates.
        % Preserve connectivity and associated mesh data.
        mesh.points = vertices_scaled;
        
        %% 3. Write the scaled mesh
        vtkWrite(mesh, char(scaled_mesh_file), 'false', 'ascii');
    
        fprintf('Saved scaled VTP mesh:\n%s\n', scaled_mesh_file);
    
    
        
        % Read VTP file
        vtp_text = fileread(scaled_mesh_file);
    
        % Identify the first and last delimiters
        start_tag = '<InformationKey name="L2_NORM_RANGE" location="vtkDataArray" length="2">';
        end_tag   = '</InformationKey>';
    
        start_idx = strfind(vtp_text, start_tag);
        end_idx   = strfind(vtp_text, end_tag);
    
        if ~isempty(start_idx) && ~isempty(end_idx)
    
            % Select the first matching block
            start_idx = start_idx(1);
            end_candidates = end_idx(end_idx >= start_idx);
            end_idx = end_candidates(1) + length(end_tag) - 1;
    
            % Remove everything from the opening to closing tag
            vtp_text(start_idx:end_idx) = [];
    
            % Write the modified file
            fid = fopen(scaled_mesh_file, 'w');
    
            if fid == -1
                error('Could not open VTP file for writing.');
            end
    
            fwrite(fid, vtp_text, 'char');
            fclose(fid);
        end
    end
    
    % Change file names in the opensim model file
    all_bodies = dom.getElementsByTagName('Body');
    bone_names{5} = 'toes_r';
    for k = 0:all_bodies.getLength - 1
        this_body = all_bodies.item(k);
        % Get the XML body name
        xml_body_name = char(this_body.getAttribute('name'));
    
        % Match body name to our bone names (case-insensitive)
        for b = 1:length(bone_names)
            if strcmpi(xml_body_name, bone_names{b})
                % 1. Update FrameGeometry (the axes display)
                attached_geoms = this_body.getElementsByTagName('attached_geometry');
                if attached_geoms.getLength > 0
                    % Get the mesh_file XML element
                    mesh_file_node = attached_geoms.item(0).getElementsByTagName('mesh_file').item(0);
                    if ~isempty(mesh_file_node)
                        % Store a path relative to the output .osim file
                        if b>3
                        mesh_filename = strcat(vtpOutputNames{b+1},'_scaled.vtp');
                        else
                        mesh_filename = strcat(vtpOutputNames{b},'_scaled.vtp');
                        end
                        % Update the XML element
                        mesh_file_node.getFirstChild.setData(mesh_filename);
                    
                        fprintf('Updated mesh path: %s\n', mesh_filename);
                    end
                    if b == 3
                        mesh_file_node_2 = attached_geoms.item(0).getElementsByTagName('mesh_file').item(1);
                        mesh_filename_2 = strcat(vtpOutputNames{b+1},'_scaled.vtp');
                        mesh_file_node_2.getFirstChild.setData(mesh_filename_2);
                    end
                end
            end
        end
    end
    
    %% ---- SCALE MUSCLE PATH POINTS ----
    % Assumptions:
    % 1. Muscle path-point coordinates are in their parent body's local frame.
    % 2. body_centroids_local contains each mesh centroid expressed in
    %    the corresponding OpenSim body-local frame.
    % 3. Centroids and path-point coordinates use the same units.
    
    % ---- CALCULATE BODY CENTROIDS FROM VTP MESHES ----
    
    body_names = {'Pelvis', 'Thigh_r', 'Leg_r', 'Foot_r', 'Toes_r'};
    vtpBodyNames = {'pelvis','femur','tibia','foot','toes'};
    
    body_centroids_local = struct();
    
    for i = 1:numel(body_names)
    
        body_name = body_names{i};
    
        % Path to the corresponding VTP mesh
        vtp_file = fullfile(scaled_mesh_folder,strcat(vtpBodyNames{i},'_scaled.vtp'));
    
        % Read mesh using the MATLAB VTK toolbox
        mesh = vtkRead(char(vtp_file));
    
        if ~isfield(mesh, 'points') || isempty(mesh.points)
            error('No vertices found in VTP mesh: %s', vtp_file);
        end
    
        % Mesh vertices in their existing coordinate system
        vertices = double(mesh.points);
    
        % Calculate vertex centroid
        body_centroids_local.(body_name) = mean(vertices, 1);
    
        fprintf('%s centroid: [%.6f %.6f %.6f]\n', ...
            body_name, ...
            body_centroids_local.(body_name)(1), ...
            body_centroids_local.(body_name)(2), ...
            body_centroids_local.(body_name)(3));
    end
    
    % Map body names to their corresponding scaling factors.
    body_scale_map = struct();
    body_scale_map.Pelvis = bone_scales.pelvis;
    body_scale_map.Thigh_r = bone_scales.thigh_r;
    body_scale_map.Leg_r = bone_scales.leg_r;
    body_scale_map.Foot_r = bone_scales.foot_r;
    body_scale_map.Toes_r = bone_scales.foot_r;
    
    % These centroids must be calculated in the local body frames.
    % Example structure:
    % body_centroids_local.thigh_r = [cx, cy, cz];
    
    point_types = {'PathPoint', ...
                   'ConditionalPathPoint', ...
                   'MovingPathPoint'};
    
    for p = 1:numel(point_types)
    
        all_points = dom.getElementsByTagName(point_types{p});
    
        for k = 0:all_points.getLength-1
    
            point_node = all_points.item(k);
    
            % Identify the parent frame.
            frame_node = point_node.getElementsByTagName( ...
                'socket_parent_frame').item(0);
    
            if ~isempty(frame_node) && ~isempty(frame_node.getFirstChild)
                frame_name = char(frame_node.getFirstChild.getData());
    
                % Handle paths such as /bodyset/thigh_r
                parts = strsplit(strtrim(frame_name), '/');
                body_name = parts{end};
    
            else
                % Compatibility with older OpenSim XML formats.
                body_node = point_node.getElementsByTagName('body').item(0);
    
                if isempty(body_node) || isempty(body_node.getFirstChild)
                    continue
                end
    
                body_name = strtrim(char(body_node.getFirstChild.getData()));
            end
    
            % Only process points whose parent body has defined scaling.
            if ~isfield(body_scale_map, body_name)
                continue
            end
    
            if ~isfield(body_centroids_local, body_name)
                warning('Missing local centroid for body %s', body_name);
                continue
            end
    
            % Read original local coordinates.
            location_node = point_node.getElementsByTagName( ...
                'location').item(0);
    
            if isempty(location_node) || isempty(location_node.getFirstChild)
                continue
            end
    
            xyz = sscanf(char(location_node.getFirstChild.getData()), ...
                '%f', 3)';
    
            if numel(xyz) ~= 3
                continue
            end
    
            % Centroid and scale for this parent body.
            c = body_centroids_local.(body_name)*0.01;
            s = body_scale_map.(body_name);
    
            % Scale around the body-specific mesh centroid.
            xyz_scaled = (xyz - c) .* s + c;
    
            % Write the new coordinates.
            new_location = sprintf('%.10f %.10f %.10f', ...
                xyz_scaled(1), xyz_scaled(2), xyz_scaled(3));
    
            location_node.getFirstChild.setData(new_location);
    
            fprintf('Scaled %s point on %s\n', ...
                point_types{p}, body_name);
    
        end
    end
end

%% Save the modified XML
if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

xmlwrite(output_file, dom);
fprintf('Model saved to:\n%s\n', output_file);
