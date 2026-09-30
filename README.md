# Mouse-Gait-OpenSim-Pipeline
A MATLAB/OpenSim pipeline for processing mouse gait data from Vicon motion capture and estimating muscle and joint forces using multibody dynamics simulations.

<img width="1887" height="683" alt="image" src="https://github.com/user-attachments/assets/557e77a1-970b-4ab4-be65-3a535cf2ba55" />

## Features

Register segmented mouse hindlimb bones and experimental marker positions to a template OpenSim model

Update marker locations in the OpenSim model

Import Vicon C3D files

Clean and filter marker trajectories

Detect stance phases automatically

Align experimental markers to a mouse OpenSim model

Export normalized .trc files

Process force plate data

Export OpenSim-compatible .mot and .xml files

Prepare inputs for inverse kinematics and inverse dynamics

Run inverse kinematics

Run OpenSim multibody dynamics simulations to estimate muscle and joint forces

## Requirements

MATLAB R2024a or newer

OpenSim 4.4 or newer

BTK Toolbox for MATLAB 0.3.0 or newer

## Installation

Clone the repository

git clone https://github.com/YOURNAME/Mouse-Gait-OpenSim-Pipeline.git

Add BTK and OpenSim to the MATLAB path.

Open MATLAB in the repository folder.

## Running the pipeline

registration_CT_opensim

c3d_to_trc

inverse_kinematic

MBD_simulations

## Outputs

*.osim — updated opensim model with personalised marker coordinates 

*.trc — normalized marker trajectories

*.mot — ground reaction forces

*.xml — OpenSim external loads

*.mot — joint angles from inverse kinematic

*.sto — muscle and joint forces

## Citation

If you use this pipeline, please cite the repository DOI (generated via Zenodo).

## License

Released under the MIT License.
