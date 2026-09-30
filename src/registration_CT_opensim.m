clc; clear; close all;

%% ==============================
% CONFIGURATION
% ==============================
% Mouse number
mouse_age = "22weeks";
mouse_name = "Cage51L";

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
rot_adjust.foot_r   = [0, 180, 180];

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
    if i==1
        scale_xyz = ct_dim ./ os_dim;
    end
   % if i==4
   %     scale_xyz(:) = 1;
   % elseif i==8
   %     scale_xyz(:) = 1;
   % end
    bone_scales.(name) = scale_xyz./100;
    
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
    plot3(results(i,1), results(i,2), results(i,3), 'ko', 'MarkerSize', 20, 'MarkerFaceColor', 'm');
    % text(results(i,1), results(i,2), results(i,3), m_name, 'FontSize', 20);
end

%% ---- UPDATE OSIM FILE ----
cd('../../')

%% ---- UPDATE OSIM FILE (MARKERS + SCALING) ----
% Load the XML document
dom = xmlread(osim_input_file);
all_markers = dom.getElementsByTagName('Marker');

% 1. Update Markers 
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

% 2. Update Body Scaling
% OpenSim stores scaling in <ScaleSet> usually, but for a simplified 
% direct mesh-to-model personalization, we target the <Body> scale or 
% simply report these for the Scale Tool. 
% If you want to modify the Body 'mass_center' or 'inertia' it's complex, 
% but here is how you find the Body to update its 'attached_geometry' scale:

% all_bodies = dom.getElementsByTagName('Body');
% for k = 0:all_bodies.getLength - 1
%     this_body = all_bodies.item(k);
%     % Get the XML body name
%     xml_body_name = char(this_body.getAttribute('name'));
% 
%     % Match body name to our bone names (case-insensitive)
%     for b = 1:length(bone_names)
%         if strcmpi(xml_body_name, bone_names{b})
%             s = bone_scales.(bone_names{b});
%             s_str = sprintf('%.10f %.10f %.10f', s(1), s(2), s(3));
% 
%             % 1. Update FrameGeometry (the axes display)
%             frame_geoms = this_body.getElementsByTagName('FrameGeometry');
%             if frame_geoms.getLength > 0
%                 sf_node = frame_geoms.item(0).getElementsByTagName('scale_factors').item(0);
%                 if ~isempty(sf_node)
%                     sf_node.getFirstChild.setData(s_str);
%                 end
%             end
% 
%             % 2. Update all Mesh geometries (Pelvis_l, Pelvis_r, etc.)
%             % We look for 'Mesh' tags specifically
%             meshes = this_body.getElementsByTagName('Mesh');
%             if meshes.getLength > 0
%                 fprintf("Body %s: Found %d meshes to scale.\n", xml_body_name, meshes.getLength);
%                 for m = 0:meshes.getLength - 1
%                     this_mesh = meshes.item(m);
%                     % Update <scale_factors>
%                     sf_nodes = this_mesh.getElementsByTagName('scale_factors');
%                     if sf_nodes.getLength > 0
%                         sf_nodes.item(0).getFirstChild.setData(s_str);
% 
%                         % Print verification
%                         m_file = char(this_mesh.getElementsByTagName('mesh_file').item(0).getFirstChild.getData());
%                         fprintf("  -> Updated scale for [%s] to [%s]\n", m_file, s_str);
%                     end
%                 end
%             else
%                 fprintf("Body %s: No <Mesh> tags found.\n", xml_body_name);
%             end
%         end
%     end
% end

% Save the modified XML
if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

xmlwrite(output_file, dom);
fprintf('Model saved to:\n%s\n', output_file);
