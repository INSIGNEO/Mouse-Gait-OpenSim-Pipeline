clear; close all; clc

import org.opensim.modeling.*

%% Inverse kinematic
% 1. Setup File Paths
mouse_age = "20weeks";
mouse_name = "Cage51L";
trcName = 'cage51L20weeks1_StanceNorm_02.trc';

modelFile  = strcat(mouse_age,'/',mouse_name,'/Model_data/Model_mouse_right_markers.osim');
trcFile    = strcat(mouse_age,'/',mouse_name,'/Model_data/',trcName);
resultsDir = strcat(mouse_age,'/',mouse_name,'/Model_data');

if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end

% 2. Initialize Model
model = Model(modelFile);
model.initSystem();

% 3. Initialize IK Tool
ikTool = InverseKinematicsTool();
ikTool.setModel(model);
ikTool.setMarkerDataFileName(trcFile);
ikTool.setResultsDir(resultsDir);

% 4. Custom Marker and Coordinate Weighting
% Define your weights in a Map or Cell array
markerWeights = {
    'Right_Illiac', 1.0; 
    'Right_Hip', 1.0;
    'Right_Knee', 1.0;
    'Right_Ankle', 1.0;
    'Right_Met', 1.0;
    'Left_Illiac', 1.0; 
};

% Clear existing tasks and add new ones
ikTaskSet = ikTool.getIKTaskSet();
ikTaskSet.clearAndDestroy();

for i = 1:size(markerWeights, 1)
    name = markerWeights{i,1};
    weight = markerWeights{i,2};

    newMarkerTask = IKMarkerTask();
    newMarkerTask.setName(name);
    newMarkerTask.setApply(true);
    newMarkerTask.setWeight(weight);

    ikTaskSet.adoptAndAppend(newMarkerTask);
end

% 5. Lock/Unlock Degrees of Freedom (Coordinate Tasks)
% Locking/Unlocking a coordinate forces IK to ignore or consider the marker data for that DOF
coordsToLock = {'zR_TorsoToPelvis', 'MTP_flexion_r', 'Hip_rotation_r', 'Knee_rotation_r'};
coordsToUnlock = {'PelvisToGround_Rx', 'PelvisToGround_Ry'};

for i = 1:length(coordsToLock)
    coordName = coordsToLock{i};
    if model.getCoordinateSet().contains(coordName)
        % Set the coordinate to locked in the model itself
        model.getCoordinateSet().get(coordName).set_locked(true);
    end
end

for i = 1:length(coordsToUnlock)
    coordName = coordsToUnlock{i};
    if model.getCoordinateSet().contains(coordName)
        % Set the coordinate to unlocked in the model itself
        model.getCoordinateSet().get(coordName).set_locked(false);
    end
end

% 6. Set Time Range and Output Name
markerData = MarkerData(trcFile);
ikTool.setStartTime(markerData.getStartFrameTime());
ikTool.setEndTime(markerData.getLastFrameTime());

% Specify the name of the output .mot file
outputMotionFile = strcat(mouse_age,'/',mouse_name,'/Model_data/ik_results_output.mot');
ikTool.setOutputMotionFileName(outputMotionFile);

% 7. Run IK
disp(['Running IK. Output will be saved to: ', outputMotionFile]);

% Save the setup file (.xml) for future record/OpenSim GUI use
ikTool.print(fullfile(resultsDir, 'ik_setup_settings.xml')); 

% Execute
ikTool.run();

disp('IK Processing Complete. The .mot file has been generated.');

%% --- 0. Settings ---
Scale = 0.01; 
InDegrees = true;
motionFile = strcat(mouse_age,'/',mouse_name,'/Model_data/ik_results_output.mot');
grfFile = strcat(mouse_age,'/',mouse_name,'/Model_data/', trcName(1:end-4),'_grf.mot');
grfXML = strcat(mouse_age,'/',mouse_name,'/Model_data/', trcName(1:end-4),'_grf.xml');
markerRadius = 0.0005; 

% Force plate visualization parameters
plateTName = strcat(mouse_age,'/',mouse_name,'/Model_data/', trcName(1:end-4),'_TranslationMatrix.mat');
plateRName = strcat(mouse_age,'/',mouse_name,'/Model_data/', trcName(1:end-4),'_RotationMatrix.mat');
plateT = load(plateTName).T/1000;
plateR = load(plateRName).R;
plateDimension = 152.4/1000; % 6 inch 
plateHeight = 0; % Distance between the wand markers and the plate
platePoint1 = [-0.5*plateDimension, 0.5*plateDimension, -plateHeight];
platePoint2 = [0.5*plateDimension, 0.5*plateDimension, -plateHeight];
platePoint3 = [0.5*plateDimension, -0.5*plateDimension, -plateHeight];
platePoint4 = [-0.5*plateDimension, -0.5*plateDimension, -plateHeight];
plateVertices = [platePoint1; platePoint2; platePoint3; platePoint4];

plateVerticesNew = plateVertices*plateR' + plateT';

% Read the GRF file
grfTable = readtable(grfFile,'FileType','text','Delimiter','\t','HeaderLines',4);

grfTime = grfTable.time;

GRF = [ ...
    grfTable.("x1_ground_force_vx"), ...
    grfTable.("x1_ground_force_vy"), ...
    grfTable.("x1_ground_force_vz")];

COP = [ ...
    grfTable.("x1_ground_force_px"), ...
    grfTable.("x1_ground_force_py"), ...
    grfTable.("x1_ground_force_pz")];

Torque = [ ...
    grfTable.("x1_ground_torque_x"), ...
    grfTable.("x1_ground_torque_y"), ...
    grfTable.("x1_ground_torque_z")];

% Remove some GRF components
% GRF(:,1) = 0; % anteroposterior
GRF(:,3) = 0; % mediolateral

% Correct CoP position
% COP(:,1) = COP(:,1)./10;
% COP(:,2) = COP(:,2) - 0.018;
% COP(:,3) = COP(:,3)./10;

% Set the moment at zero 
Torque = zeros(length(GRF),3);

% --- NEW: Video Settings ---
saveVideo = true; % Set to false if you just want to preview
videoName = 'Hindlimb_Animation.avi';
v = VideoWriter(videoName, 'Uncompressed AVI');
v.FrameRate = 10; % Adjust based on your pause(0.1)
if saveVideo, open(v); end

% --- 1. Load Model & Motion ---
m(1).model = Model(modelFile);
m(1).state = m(1).model.initSystem();
m(1).color = [0.8, 0.8, 0.8]; 
m(1).alpha = 0.8; 

motion = Storage(motionFile);
labels = motion.getColumnLabels();

trcData = MarkerData(trcFile);
nExpMarkers = trcData.getNumMarkers();
nFrames = trcData.getNumFrames();
expMarkerNames = {'Right_Illiac', 'Right_Hip', 'Right_Knee', 'Right_Ankle', 'Right_Met',...
                  'Left_Illiac'}; % Must match TRC column names

% Susbample GRF data
nKin = nFrames;
nGRF = size(GRF,1);

xOld = linspace(0,1,nGRF);
xNew = linspace(0,1,nKin);

GRF_sub = interp1(xOld, GRF, xNew, 'linear');
COP = interp1(xOld, COP, xNew, 'linear');

% --- 2. Setup Visualization ---
bodies = {'Pelvis', 'Thigh_r', 'Leg_r', 'Foot_r', 'Toes_r'}; 
stls   = {'Geometry/Pelvis_r.stl', 'Geometry/femur_r.stl', ...
          'Geometry/tibfib_r.stl', 'Geometry/foot_r.stl', 'Geometry/phalanges_r.stl'};

figure('Color', 'w'); hold on; axis equal; view(2); light; grid off;

% Initialize Meshes and Capture Neutral Transforms
m(1).model.realizePosition(m(1).state);
for i = 1:length(bodies)
    if exist(stls{i}, 'file')
        TR = stlread(stls{i});
        m(1).mesh(i).V_orig = TR.Points * Scale;
        m(1).mesh(i).bodyName = bodies{i};
        m(1).mesh(i).handle = patch('Vertices', m(1).mesh(i).V_orig, ...
            'Faces', TR.ConnectivityList, 'FaceColor', m(1).color, ...
            'EdgeColor', 'none', 'FaceAlpha', m(1).alpha, 'FaceLighting', 'gouraud');
        
        body = m(1).model.getBodySet().get(bodies{i});
        T_Neutral = body.getTransformInGround(m(1).state);
        R_n = T_Neutral.R().asMat33(); p_n = T_Neutral.p();
        M_Neutral = [R_n.get(0,0), R_n.get(0,1), R_n.get(0,2), p_n.get(0);
                     R_n.get(1,0), R_n.get(1,1), R_n.get(1,2), p_n.get(1);
                     R_n.get(2,0), R_n.get(2,1), R_n.get(2,2), p_n.get(2);
                     0, 0, 0, 1];
        m(1).mesh(i).M_Neutral_Inv = inv(M_Neutral);
    end
end

% --- 3. Setup Spheres ---
[uX, uY, uZ] = sphere(10);
uX = uX*markerRadius; uY = uY*markerRadius; uZ = uZ*markerRadius;

% Model Markers (Pink)
markerSet = m(1).model.getMarkerSet();
for k = 1:markerSet.getSize()
    m(1).marker(k).handle = surf(uX, uY, uZ, 'FaceColor', [1.0, 0.41, 0.70], 'EdgeColor', 'none');
end

% Experimental Markers (Blue)
for k = 1:length(expMarkerNames)
    m(1).expMarker(k).handle = surf(uX, uY, uZ, 'FaceColor', [0 0 1], 'EdgeColor', 'none');
end

%% Ground reaction force arrow
GRFscale = 0.1; % adjust visually
grfHandle = quiver3( 0,0,0,0,0,0,'LineWidth',2,'MaxHeadSize',0.5);

% Force plate
patch('Faces', [1, 2, 3, 4], 'Vertices', plateVerticesNew, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.5);

% --- 4. ANIMATION LOOP ---
CoP = zeros(nFrames,3);
for f = 0:nFrames-1
    dataVector = motion.getStateVector(f).getData();
    
    % [A through F: Your existing update logic remains exactly as is]
    
    % A. Global Translation
    globalOffset = [dataVector.get(4), dataVector.get(5), dataVector.get(6)];
    
    % B. Update Model Rotations (SKIP TRANSLATIONS)
    for j = 0:labels.size() - 1
        colName = char(labels.get(j));
        if m(1).model.getCoordinateSet().contains(colName)
            if ~contains(colName, '_T') && ~contains(colName, 'pos')
                val = dataVector.get(j - 1); 
                if InDegrees, val = val * (pi/180); end 
                m(1).model.getCoordinateSet().get(colName).setValue(m(1).state, val);
            else
                m(1).model.getCoordinateSet().get(colName).setValue(m(1).state, 0);
            end
        end
    end
    
    m(1).model.realizePosition(m(1).state);
    
    % D. Update Bone Meshes
    for i = 1:length(m(1).mesh)
        body = m(1).model.getBodySet().get(m(1).mesh(i).bodyName);
        T = body.getTransformInGround(m(1).state);
        R_curr = T.R().asMat33(); p_curr = T.p();
        M_Abs = [R_curr.get(0,0), R_curr.get(0,1), R_curr.get(0,2), p_curr.get(0);
                 R_curr.get(1,0), R_curr.get(1,1), R_curr.get(1,2), p_curr.get(1);
                 R_curr.get(2,0), R_curr.get(2,1), R_curr.get(2,2), p_curr.get(2);
                 0, 0, 0, 1];
        M_Final = M_Abs * m(1).mesh(i).M_Neutral_Inv;
        new_V = (m(1).mesh(i).V_orig * M_Final(1:3, 1:3)') + (M_Final(1:3, 4)' + globalOffset);
        if i==4
            CoP(f+1,:) = estimateFootCoP(new_V);
        end
        set(m(1).mesh(i).handle, 'Vertices', new_V);
    end
    
    % E. Update Model Markers
    for k = 1:markerSet.getSize()
        loc = markerSet.get(k-1).getLocationInGround(m(1).state);
        set(m(1).marker(k).handle, 'XData', uX + loc.get(0) + globalOffset(1), ...
                                   'YData', uY + loc.get(1) + globalOffset(2), ...
                                   'ZData', uZ + loc.get(2) + globalOffset(3));
    end

    % F. Update Experimental Markers
    frameData = trcData.getFrame(f);
    for k = 1:length(expMarkerNames)
        mIdx = trcData.getMarkerIndex(expMarkerNames{k});
        p_exp_raw = frameData.getMarker(mIdx);
        set(m(1).expMarker(k).handle, 'XData', uX + p_exp_raw.get(0), ...
                                      'YData', uY + p_exp_raw.get(1), ...
                                      'ZData', uZ + p_exp_raw.get(2));
    end

    %% G. Update GRF

    if f+1 <= size(GRF_sub,1)
    
        cop = CoP(f+1,:);
    
        force = GRF_sub(f+1,:);
    
        set(grfHandle,...
            'XData',cop(1),...
            'YData',cop(2),...
            'ZData',cop(3),...
            'UData',force(1)*GRFscale,...
            'VData',force(2)*GRFscale,...
            'WData',force(3)*GRFscale);
    
    end

    drawnow limitrate;
    
    % --- NEW: Capture Frame for Video ---
    if saveVideo
        frame = getframe(gcf); % Captures the current figure
        writeVideo(v, frame);
    end
    
    pause(0.1);
end

% --- NEW: Finalize Video ---
if saveVideo
    close(v);
    fprintf('Video saved as %s\n', videoName);
end

%% Calculate the RMSE from the inverse kinematic
% Preallocate
nMarkers = length(expMarkerNames);
markerErrors = zeros(nFrames,nMarkers);

% Loop through all frames again
for f = 0:nFrames-1

    % Get IK state
    dataVector = motion.getStateVector(f).getData();

    % Apply coordinates
    for j = 0:labels.size()-1
        colName = char(labels.get(j));

        if m(1).model.getCoordinateSet().contains(colName)

            if ~contains(colName,'_T') && ~contains(colName,'pos')
                val = dataVector.get(j-1);

                if InDegrees
                    val = val*pi/180;
                end

                m(1).model.getCoordinateSet().get(colName).setValue(m(1).state,val);

            else
                m(1).model.getCoordinateSet().get(colName).setValue(m(1).state,0);
            end
        end
    end

    m(1).model.realizePosition(m(1).state);

    % Global translation used in animation
    globalOffset = [dataVector.get(4), ...
                    dataVector.get(5), ...
                    dataVector.get(6)];

    % Experimental marker positions
    frameData = trcData.getFrame(f);

    for k = 1:nMarkers

        % Experimental marker
        mIdx = trcData.getMarkerIndex(expMarkerNames{k});
        pExp = frameData.getMarker(mIdx);

        expPos = [pExp.get(0), ...
                  pExp.get(1), ...
                  pExp.get(2)];

        % Corresponding model marker
        modelMarker = markerSet.get(expMarkerNames{k});

        pModel = modelMarker.getLocationInGround(m(1).state);

        modelPos = [pModel.get(0), ...
                    pModel.get(1), ...
                    pModel.get(2)] + globalOffset;

        % Euclidean distance error
        markerErrors(f+1,k) = norm(expPos - modelPos);

    end
end

%% RMSE per marker

RMSE_marker = sqrt(mean(markerErrors.^2,1));

fprintf('\nRMSE per marker:\n');
for k = 1:nMarkers
    fprintf('%s : %.4f mm\n', ...
        expMarkerNames{k}, RMSE_marker(k)*1000);
end

%% Global RMSE

Global_RMSE = sqrt(mean(markerErrors(:).^2));

fprintf('\nGlobal IK RMSE = %.4f mm\n',Global_RMSE*1000);

%% Plot RMSE per marker

figure('Color','w');

subplot(2,1,1)
plot(markerErrors*1000)
ylabel('RMSE (mm)')
legend(expMarkerNames)
title('RMSE of each marker')
grid on

%% Plot RMS error through time

frameRMSE = sqrt(mean(markerErrors.^2,2));

subplot(2,1,2)
plot(frameRMSE*1000,'LineWidth',2)
xlabel('Frame')
ylabel('RMSE (mm)')
title(sprintf('Global RMSE through time (mean = %.2f mm)',...
      Global_RMSE*1000))
grid on

%% Update GRF .mot file with new CoP

% Interpolate COP data
CoP = interp1(xNew, CoP, xOld, 'linear');

% Replace Force columns
grfTable.("x1_ground_force_vx") = GRF(:,1);
grfTable.("x1_ground_force_vy") = GRF(:,2);
grfTable.("x1_ground_force_vz") = GRF(:,3);

% Replace CoP columns
grfTable.("x1_ground_force_px") = CoP(:,1);
grfTable.("x1_ground_force_py") = CoP(:,2);
grfTable.("x1_ground_force_pz") = CoP(:,3);

% Replace Torque columns
grfTable.("x1_ground_torque_x") = Torque(:,1);
grfTable.("x1_ground_torque_y") = Torque(:,2);
grfTable.("x1_ground_torque_z") = Torque(:,3);

COP(:,1) = COP(:,1)./10;
COP(:,2) = COP(:,2) - 0.018;
COP(:,3) = COP(:,3)./10;

%% Read original header

fid = fopen(grfFile,'r');

headerLines = cell(4,1);

for i = 1:4
    headerLines{i} = fgetl(fid);
end

fclose(fid);

%% Write updated file

fid = fopen(grfFile,'w');

for i = 1:length(headerLines)
    fprintf(fid,'%s\n',headerLines{i});
end

% Write column names
varNames = grfTable.Properties.VariableNames;

for i = 1:length(varNames)-1
    fprintf(fid,'%s\t',varNames{i});
end
fprintf(fid,'%s\n',varNames{end});

% Write data
formatSpec = [repmat('%.8f\t',1,width(grfTable)-1) '%.8f\n'];

for r = 1:height(grfTable)
    fprintf(fid,formatSpec,grfTable{r,:});
end

fclose(fid);

fprintf('GRF file updated with new CoP locations:\n%s\n',grfFile);

% Copy the grf file
newGRFMotName = strcat(mouse_age,'/',mouse_name,'/Model_data/grf_output.mot');
newGRFXMLName = strcat(mouse_age,'/',mouse_name,'/Model_data/grf_output.xml');
copyfile(grfFile, newGRFMotName);
copyfile(grfXML, newGRFXMLName);
