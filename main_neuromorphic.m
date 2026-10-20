 %% =========================================================
%  NEUROMORPHIC PHOTONIC DEVICE SIMULATION
%  Based on RTD-VCSEL PRL Node + Photonic SNN
%  Papers: Hejda et al. 2022, Owen-Newns et al. 2023
%
%  Run this file to execute all simulations in sequence.
%  Each module can also be run independently.
%% =========================================================

clear; clc; close all;

fprintf('==============================================\n');
fprintf('  NEUROMORPHIC RTD-VCSEL PHOTONIC SIMULATION\n');
fprintf('==============================================\n\n');

%% --- MODULE 1: RTD I-V Characteristic ---
fprintf('[1/5] Plotting RTD I-V characteristic...\n');
rtd_iv_characteristic;
pause(1);

%% --- MODULE 2: RTD Spiking Neuron (Single Node) ---
fprintf('[2/5] Simulating single RTD spiking neuron...\n');
rtd_spiking_neuron;
pause(1);

%% --- MODULE 3: Excitability Threshold Sweep ---
fprintf('[3/5] Sweeping bias voltage (excitability threshold)...\n');
rtd_threshold_sweep;
pause(1);

%% --- MODULE 4: PRL Node (Photodetector + RTD + VCSEL) ---
fprintf('[4/5] Simulating full PRL optoelectronic node...\n');
prl_node_simulation;
pause(1);

%% --- MODULE 5: Photonic SNN - Iris Classification ---
fprintf('[5/5] Running photonic SNN (reservoir computing)...\n');
photonic_snn_iris;

fprintf('\n==============================================\n');
fprintf('  ALL SIMULATIONS COMPLETE\n');
fprintf('==============================================\n');
