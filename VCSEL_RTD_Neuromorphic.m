clear; clc; close all;

fprintf('==============================================\n');
fprintf('  VCSEL-RTD Neuromorphic Device Simulator\n');
fprintf('  3-Layer VCSEL | 3-Barrier RTD | 100 Hz\n');
fprintf('==============================================\n\n');
q    = 1.602e-19;   
h    = 6.626e-34;   
hbar = h/(2*pi);     
kb   = 1.381e-23;   
c0   = 3e8;         
m0   = 9.109e-31;   
eps0 = 8.854e-12;    % Permittivity of free space [F/m]
VCSEL.lambda    = 850e-9;          % Lasing wavelength [m]
VCSEL.nu        = c0/VCSEL.lambda; % Optical frequency [Hz]
VCSEL.n_layers  = 3;               % Number of QW layers
VCSEL.Lz        = 8e-9;            % QW thickness per layer [m]
VCSEL.L_total   = VCSEL.n_layers * VCSEL.Lz; % Total active length
VCSEL.V_active  = pi*(3e-6)^2 * VCSEL.L_total; % Active volume [m^3]
VCSEL.n_r       = 3.5;             % Refractive index (GaAs)
VCSEL.Gamma     = 0.04 * VCSEL.n_layers; % Confinement factor (×layers)
VCSEL.g0        = 1200;            % Material gain coeff [1/m]
VCSEL.N_tr      = 1.8e24;          % Transparency carrier density [m^-3]
VCSEL.epsilon   = 1e-23;           % Gain compression factor [m^3]
VCSEL.tau_n     = 2e-9;            % Carrier lifetime [s]
VCSEL.tau_p     = 1.5e-12;         % Photon lifetime [s]
VCSEL.beta_sp   = 1e-4;            % Spontaneous emission factor
VCSEL.vg        = c0/VCSEL.n_r;    % Group velocity [m/s]
VCSEL.I_th_est  = 1.5e-3;          % Estimated threshold current [A]
VCSEL.eta_i     = 0.85;            % Internal quantum efficiency
VCSEL.T         = 300;             % Operating temperature [K]

% DBR Mirror parameters
VCSEL.R_top     = 0.998;           % Top DBR reflectivity
VCSEL.R_bot     = 0.9998;          % Bottom DBR reflectivity
VCSEL.alpha_i   = 10;              % Internal loss [1/cm] → [1/m]
VCSEL.L_cav     = 1e-6;            % Cavity length [m]

% Compute photon lifetime from mirror losses
alpha_m = (1/(2*VCSEL.L_cav)) * log(1/(VCSEL.R_top*VCSEL.R_bot));
VCSEL.tau_p = 1 / (VCSEL.vg * (alpha_m + VCSEL.alpha_i*100));

fprintf('[VCSEL] 3-Layer QW Structure:\n');
fprintf('  Wavelength     : %.0f nm\n', VCSEL.lambda*1e9);
fprintf('  QW Layers      : %d × %.0f nm\n', VCSEL.n_layers, VCSEL.Lz*1e9);
fprintf('  Confinement Γ  : %.4f\n', VCSEL.Gamma);
fprintf('  Photon lifetime: %.2f ps\n', VCSEL.tau_p*1e12);
fprintf('  Active Volume  : %.2e m³\n\n', VCSEL.V_active);

%% ===== RTD STRUCTURE PARAMETERS (3-Barrier) =====
% InGaAs/AlAs triple-barrier RTD
RTD.n_barriers  = 3;               % Number of tunneling barriers
RTD.t_barrier   = 1.7e-9;          % Barrier thickness [m] (AlAs)
RTD.t_well      = 5e-9;            % Well thickness [m] (InGaAs)
RTD.m_eff       = 0.042 * m0;      % Effective mass InGaAs [kg]
RTD.m_bar       = 0.15 * m0;       % Effective mass AlAs barrier [kg]
RTD.E_b         = 0.55;            % Barrier height [eV] → [J]
RTD.E_b_J       = RTD.E_b * q;
RTD.area        = (3e-6)^2 * pi;   % RTD mesa area [m²]

% NDR (Negative Differential Resistance) parameters
RTD.V_peak      = 0.18;            % Peak voltage [V]
RTD.V_valley    = 0.42;            % Valley voltage [V]
RTD.J_peak      = 8e4;             % Peak current density [A/m²]
RTD.J_valley    = 1.2e4;           % Valley current density [A/m²]
RTD.PVR         = RTD.J_peak/RTD.J_valley; % Peak-to-valley ratio
RTD.C_rtd       = 15e-15;          % RTD capacitance [F]
RTD.R_series    = 5;               % Series resistance [Ω]
RTD.V_bias      = 0.28;            % Operating bias point [V] (in NDR)

% Compute resonant level energies (particle in box approx)
for i = 1:RTD.n_barriers-1
    RTD.E_res(i) = (hbar^2 * pi^2 * i^2) / (2 * RTD.m_eff * RTD.t_well^2 * q);
end

fprintf('[RTD] 3-Barrier Tunnel Structure:\n');
fprintf('  Barriers       : %d × %.1f nm AlAs\n', RTD.n_barriers, RTD.t_barrier*1e9);
fprintf('  Wells          : %d × %.1f nm InGaAs\n', RTD.n_barriers-1, RTD.t_well*1e9);
fprintf('  Peak Voltage   : %.2f V\n', RTD.V_peak);
fprintf('  Valley Voltage : %.2f V\n', RTD.V_valley);
fprintf('  Peak-to-Valley : %.1f\n', RTD.PVR);
fprintf('  Resonant E1    : %.4f eV\n', RTD.E_res(1));
fprintf('  Resonant E2    : %.4f eV\n\n', RTD.E_res(2));

%% ===== NEUROMORPHIC CHIP PARAMETERS =====
NEURO.f_chip     = 100;            % Chip frequency [Hz]
NEURO.T_chip     = 1/NEURO.f_chip; % Clock period [s]
NEURO.N_neurons  = 64;             % Neurons in network (8×8 array)
NEURO.N_syn      = NEURO.N_neurons^2; % Synaptic connections
NEURO.tau_mem    = 5e-3;           % Membrane time constant [s]
NEURO.V_thresh   = 0.35;           % Firing threshold voltage [V]
NEURO.V_reset    = 0.05;           % Reset voltage [V]
NEURO.V_rest     = 0.10;           % Resting voltage [V]
NEURO.t_refrac   = 2e-3;           % Refractory period [s]
NEURO.w_syn      = 0.02;           % Synaptic weight [V]
NEURO.leak       = 0.15;           % Leak conductance [normalized]

fprintf('[NEURO] 100 Hz Neuromorphic Chip:\n');
fprintf('  Clock frequency : %d Hz\n', NEURO.f_chip);
fprintf('  Neurons (8×8)   : %d\n', NEURO.N_neurons);
fprintf('  Membrane τ      : %.1f ms\n', NEURO.tau_mem*1e3);
fprintf('  Firing threshold: %.2f V\n', NEURO.V_thresh);
fprintf('  Refractory time : %.1f ms\n\n', NEURO.t_refrac*1e3);

%% ===== SIMULATION TIME SETUP =====
dt      = 0.5e-12;     % Time step [s] — sub-ps for VCSEL dynamics
t_end   = 50e-9;       % Total simulation time [s]
t       = 0:dt:t_end;
Nt      = length(t);

fprintf('[SIM] Time setup:\n');
fprintf('  dt = %.1f ps | t_end = %.0f ns | Steps = %d\n\n', ...
        dt*1e12, t_end*1e9, Nt);

%% ===== RTD I-V CHARACTERISTIC =====
fprintf('[COMPUTE] Building RTD I-V characteristic...\n');
V_rtd = linspace(0, 0.7, 1000);
I_rtd = RTD_IV(V_rtd, RTD);

%% ===== VCSEL RATE EQUATIONS SIMULATION =====
fprintf('[COMPUTE] Running VCSEL rate equations...\n');

% State variables: N = carrier density, S = photon density
N = zeros(1, Nt);   % Carrier density [m^-3]
S = zeros(1, Nt);   % Photon density [m^-3]
I_drive = zeros(1, Nt);

% Initial conditions
N(1) = VCSEL.N_tr * 0.5;
S(1) = 1e15;

% Drive current: step + RTD modulation
I_bias   = 2.5e-3;    % DC bias [A]
I_mod_amp = 0.8e-3;   % Modulation amplitude [A]
f_mod    = 5e9;        % Modulation frequency [Hz]

for i = 1:Nt
    % Drive current with RTD-like NDR pulsing
    I_drive(i) = I_bias + I_mod_amp * sin(2*pi*f_mod*t(i)) .* ...
                 (1 + 0.3*sin(2*pi*2e9*t(i)));
end

% Integrate rate equations (Euler — fast for demo; RK4 block below)
for i = 1:Nt-1
    [dN,dS] = VCSEL_RateEqs(N(i), S(i), I_drive(i), VCSEL);
    N(i+1) = max(N(i) + dt*dN, 0);
    S(i+1) = max(S(i) + dt*dS, 1e10);
end

% Optical power output
P_out = (h * VCSEL.nu * VCSEL.V_active .* S) / (2 * VCSEL.tau_p);

fprintf('  Peak optical power : %.2f mW\n', max(P_out)*1e3);
fprintf('  Avg optical power  : %.2f mW\n', mean(P_out)*1e3);

%% ===== RTD OSCILLATION SIMULATION =====
fprintf('[COMPUTE] Simulating RTD NDR oscillations...\n');

V_rtd_t = zeros(1, Nt);
I_rtd_t = zeros(1, Nt);
V_rtd_t(1) = RTD.V_bias;

% LC tank circuit oscillation in NDR region
L_ind  = 0.1e-9;    % Parasitic inductance [H]
I_L    = zeros(1, Nt);
I_L(1) = RTD_IV(RTD.V_bias, RTD) / RTD.area;

for i = 1:Nt-1
    J_t = RTD_IV(V_rtd_t(i), RTD) / RTD.area;
    I_rtd_t(i) = J_t * RTD.area;
    V_src = RTD.V_bias + 0.02*sin(2*pi*50e9*t(i));
    dV = (I_L(i) - I_rtd_t(i)) / RTD.C_rtd;
    dI = (V_src - V_rtd_t(i) - I_L(i)*RTD.R_series) / L_ind;
    V_rtd_t(i+1) = V_rtd_t(i) + dt*dV;
    V_rtd_t(i+1) = max(0, min(V_rtd_t(i+1), 0.7));
    I_L(i+1)     = I_L(i) + dt*dI;
end
I_rtd_t(Nt) = RTD_IV(V_rtd_t(Nt), RTD) / RTD.area;

%% ===== LIF NEURON (VCSEL-RTD COUPLED) SIMULATION =====
fprintf('[COMPUTE] Running LIF neuron network (100 Hz chip)...\n');

% Simulate 10 representative neurons
N_show  = 8;
t_neuro = 0 : NEURO.T_chip/100 : 10*NEURO.T_chip;
Nt_n    = length(t_neuro);
dt_n    = t_neuro(2)-t_neuro(1);

V_mem   = NEURO.V_rest * ones(N_show, Nt_n);
spikes  = zeros(N_show, Nt_n);
t_last  = -inf(N_show, 1);

% Synaptic weight matrix (random sparse)
rng(42);
W = 0.03 * rand(N_show) .* (rand(N_show) > 0.6);
W(1:N_show+1:end) = 0;  % No self-connections

% Input stimulus: Poisson-like spike train
I_stim = zeros(N_show, Nt_n);
for n = 1:N_show
    fr = 20 + 10*rand();   % 20–30 Hz input
    for k = 1:Nt_n
        if rand() < fr * dt_n
            I_stim(n,k) = NEURO.w_syn * (2 + rand());
        end
    end
end

for k = 1:Nt_n-1
    for n = 1:N_show
        refrac = (t_neuro(k) - t_last(n)) < NEURO.t_refrac;
        if refrac
            V_mem(n,k+1) = NEURO.V_reset;
        else
            % Synaptic input from other neurons
            syn_in = sum(W(n,:)' .* spikes(:,k));
            % LIF dynamics
            dV = (-(V_mem(n,k) - NEURO.V_rest) / NEURO.tau_mem + ...
                   I_stim(n,k) + syn_in) / NEURO.tau_mem;
            V_mem(n,k+1) = V_mem(n,k) + dt_n * dV;
            if V_mem(n,k+1) >= NEURO.V_thresh
                spikes(n,k+1) = 1;
                V_mem(n,k+1)  = NEURO.V_reset;
                t_last(n)     = t_neuro(k+1);
            end
        end
    end
end

% Compute firing rates
fire_rates = sum(spikes, 2) / (t_neuro(end));
fprintf('  Mean firing rate : %.1f Hz\n', mean(fire_rates));
fprintf('  Max  firing rate : %.1f Hz\n', max(fire_rates));

%% ===== ENERGY EFFICIENCY ANALYSIS =====
fprintf('[COMPUTE] Energy efficiency metrics...\n');
E_per_spike  = 0.5 * RTD.C_rtd * VCSEL.I_th_est^2 * VCSEL.tau_n;
E_photon     = h * VCSEL.nu;
eta_wall     = mean(P_out) / (I_bias * 1.8);   % Wall-plug efficiency
fprintf('  Energy per spike : %.2f fJ\n', E_per_spike*1e15);
fprintf('  Wall-plug η      : %.1f%%\n', eta_wall*100);

%% ===== PLOTTING =====
fprintf('\n[PLOT] Generating figures...\n');

fig = figure('Name','VCSEL-RTD Neuromorphic Device','Color','#0d1117',...
             'Position',[50 30 1400 900],'NumberTitle','off');

% Color palette
col_laser = '#00FFAA';
col_rtd   = '#FF6B6B';
col_neuro = '#4ECDC4';
col_drive = '#FFE66D';
col_bg    = '#0d1117';
col_panel = '#161b22';
col_text  = '#e6edf3';
col_grid  = '#30363d';

% --- Panel 1: VCSEL Photon & Carrier Density ---
ax1 = subplot(3,3,1);
t_ns = t * 1e9;
yyaxis left
plot(t_ns, S*1e-21, 'Color', col_laser, 'LineWidth', 1.4);
ylabel('Photon Density [×10²¹ m⁻³]','Color',col_laser,'FontSize',9);
yyaxis right
plot(t_ns, N*1e-24, 'Color', col_drive, 'LineWidth', 1.2, 'LineStyle','--');
ylabel('Carrier Density [×10²⁴ m⁻³]','Color',col_drive,'FontSize',9);
xlabel('Time [ns]','Color',col_text);
title('VCSEL: Photon & Carrier Dynamics','Color',col_text,'FontWeight','bold');
style_ax(ax1, col_panel, col_text, col_grid);

% --- Panel 2: Optical Power Output ---
ax2 = subplot(3,3,2);
plot(t_ns, P_out*1e3, 'Color', col_laser, 'LineWidth', 1.5);
xlabel('Time [ns]','Color',col_text);
ylabel('Optical Power [mW]','Color',col_text);
title('VCSEL Optical Output (3-QW)','Color',col_text,'FontWeight','bold');
yline(mean(P_out)*1e3,'--','Color','#888888','LineWidth',1,'Label','Mean');
style_ax(ax2, col_panel, col_text, col_grid);

% --- Panel 3: RTD I-V Characteristic ---
ax3 = subplot(3,3,3);
plot(V_rtd, I_rtd*1e3, 'Color', col_rtd, 'LineWidth', 2.2);
xline(RTD.V_peak,'--','Color','#AAAAFF','LineWidth',1.2,'Label','V_{peak}');
xline(RTD.V_valley,'--','Color','#FF9999','LineWidth',1.2,'Label','V_{valley}');
xline(RTD.V_bias,'-','Color','#FFFF00','LineWidth',1.5,'Label','Bias');
xlabel('Voltage [V]','Color',col_text);
ylabel('Current [mA]','Color',col_text);
title('RTD I-V (3-Barrier NDR)','Color',col_text,'FontWeight','bold');
style_ax(ax3, col_panel, col_text, col_grid);
text(0.05,max(I_rtd)*1e3*0.9,sprintf('PVR = %.1f',RTD.PVR),...
     'Color',col_text,'FontSize',9,'FontWeight','bold','Parent',ax3);

% --- Panel 4: RTD Voltage Oscillation ---
ax4 = subplot(3,3,4);
t_ps = t * 1e12;
plot(t_ps, V_rtd_t, 'Color', col_rtd, 'LineWidth', 1.3);
yline(RTD.V_peak,'--','Color','#AAAAFF','LineWidth',1,'Label','V_p');
yline(RTD.V_valley,'--','Color','#FF9999','LineWidth',1,'Label','V_v');
xlabel('Time [ps]','Color',col_text);
ylabel('RTD Voltage [V]','Color',col_text);
title('RTD NDR Oscillation','Color',col_text,'FontWeight','bold');
style_ax(ax4, col_panel, col_text, col_grid);

% --- Panel 5: Neuron Membrane Potentials ---
ax5 = subplot(3,3,5);
cmap = lines(N_show);
for n = 1:N_show
    plot(t_neuro*1e3, V_mem(n,:), 'Color', cmap(n,:), 'LineWidth', 0.9); hold on;
end
yline(NEURO.V_thresh,'--','Color','#FF4444','LineWidth',1.5,'Label','V_{th}');
xlabel('Time [ms]','Color',col_text);
ylabel('Membrane Voltage [V]','Color',col_text);
title('LIF Neurons (100 Hz Chip)','Color',col_text,'FontWeight','bold');
xlim([0 t_neuro(end)*1e3]); ylim([0 0.45]);
style_ax(ax5, col_panel, col_text, col_grid);

% --- Panel 6: Spike Raster Plot ---
ax6 = subplot(3,3,6);
for n = 1:N_show
    spike_times = t_neuro(logical(spikes(n,:))) * 1e3;
    if ~isempty(spike_times)
        stem(spike_times, n*ones(size(spike_times)), 'Color', cmap(n,:), ...
             'Marker','none','LineWidth',1.5); hold on;
    end
end
xlabel('Time [ms]','Color',col_text);
ylabel('Neuron Index','Color',col_text);
title('Spike Raster (100 Hz Chip)','Color',col_text,'FontWeight','bold');
xlim([0 t_neuro(end)*1e3]); ylim([0 N_show+1]);
yticks(1:N_show);
style_ax(ax6, col_panel, col_text, col_grid);

% --- Panel 7: Drive Current ---
ax7 = subplot(3,3,7);
plot(t_ns, I_drive*1e3, 'Color', col_drive, 'LineWidth', 1.2);
xlabel('Time [ns]','Color',col_text);
ylabel('Drive Current [mA]','Color',col_text);
title('VCSEL Drive Current','Color',col_text,'FontWeight','bold');
yline(VCSEL.I_th_est*1e3,'--','Color','#FF4444','LineWidth',1,'Label','I_{th}');
style_ax(ax7, col_panel, col_text, col_grid);

% --- Panel 8: Firing Rate Histogram ---
ax8 = subplot(3,3,8);
bar(1:N_show, fire_rates, 'FaceColor', col_neuro, 'EdgeColor','none');
xlabel('Neuron Index','Color',col_text);
ylabel('Firing Rate [Hz]','Color',col_text);
title('Neuron Firing Rates','Color',col_text,'FontWeight','bold');
yline(NEURO.f_chip,'--','Color','#FF4444','LineWidth',1.5,'Label','100 Hz');
style_ax(ax8, col_panel, col_text, col_grid);

% --- Panel 9: Device Architecture Diagram ---
ax9 = subplot(3,3,9);
draw_device_diagram(ax9, VCSEL, RTD, NEURO, col_text, col_laser, col_rtd, col_neuro);

% Super title
annotation(fig,'textbox',[0.02 0.96 0.96 0.04],...
    'String','VCSEL–RTD Neuromorphic Device  |  3-Layer VCSEL + 3-Barrier RTD  |  100 Hz Chip',...
    'Color','#58a6ff','FontSize',13,'FontWeight','bold',...
    'EdgeColor','none','BackgroundColor','none','HorizontalAlignment','center');

set(fig,'Color',col_bg);
drawnow;

fprintf('\n[DONE] Simulation complete!\n');
fprintf('==============================================\n');
print_summary(VCSEL, RTD, NEURO, max(P_out), mean(fire_rates), E_per_spike, eta_wall);

%% ===== SAVE FIGURE =====
saveas(fig, '/mnt/user-data/outputs/VCSEL_RTD_Neuromorphic.png');
fprintf('\n[SAVED] Figure → VCSEL_RTD_Neuromorphic.png\n');

%% ============================================================
%                     HELPER FUNCTIONS
%% ============================================================

function [dN, dS] = VCSEL_RateEqs(N, S, I, VC)
% VCSEL rate equations: carrier and photon densities
    g   = VC.vg * VC.g0 * log(max(N/VC.N_tr, 1e-3)) / (1 + VC.epsilon*S);
    dN  = (VC.eta_i * I)/(q * VC.V_active) - N/VC.tau_n - g*S;
    dS  = VC.Gamma*g*S - S/VC.tau_p + VC.Gamma*VC.beta_sp*N/VC.tau_n;
end

function I = RTD_IV(V, RTD)
% 3-Barrier RTD I-V model: Gaussian peak + valley + linear background
    q    = 1.602e-19;
    % First resonance (main NDR peak)
    dV1  = 0.045;
    J1   = RTD.J_peak * exp(-((V - RTD.V_peak).^2)/(2*dV1^2));
    % Second resonance (weak)
    dV2  = 0.03;
    J2   = 0.25*RTD.J_peak * exp(-((V - RTD.V_peak - 0.18).^2)/(2*dV2^2));
    % Valley current (thermionic + tunneling background)
    J_val = RTD.J_valley * (1 - exp(-V/0.05)) .* (1 + 0.5*(V/RTD.V_valley).^2);
    % Total current density
    J    = J1 + J2 + J_val;
    % Smooth NDR region
    in_ndr = (V >= RTD.V_peak) & (V <= RTD.V_valley);
    J_ndr  = RTD.J_peak - (RTD.J_peak - RTD.J_valley) * ...
             ((V - RTD.V_peak)/(RTD.V_valley - RTD.V_peak)).^0.7;
    J(in_ndr) = J_ndr(in_ndr);
    I = J * RTD.area;
end

function style_ax(ax, bg, fg, gc)
    set(ax, 'Color', bg, 'XColor', fg, 'YColor', fg, ...
        'GridColor', gc, 'MinorGridColor', gc, ...
        'FontSize', 8, 'FontName', 'Consolas', ...
        'Box','on','GridAlpha',0.4);
    grid(ax,'on');
end

function draw_device_diagram(ax, VCSEL, RTD, NEURO, ct, cl, cr, cn)
    cla(ax); hold(ax,'on');
    set(ax,'Color','#0d1117','XColor',ct,'YColor',ct,...
        'XTick',[],'YTick',[],'Box','on');
    
    % Top DBR
    for i=1:5
        y = 0.85 + i*0.022;
        rectangle(ax,'Position',[0.1 y 0.8 0.018],...
            'FaceColor',hex2rgb('#1f6feb'),'EdgeColor','none');
    end
    text(ax,0.5,0.97,'Top DBR','Color',ct,'HorizontalAlignment','center','FontSize',7);
    
    % 3 QW layers
    qw_cols = {cl,'#88FF88','#44FFFF'};
    for i=1:3
        y = 0.62 - (i-1)*0.075;
        rectangle(ax,'Position',[0.1 y 0.8 0.05],...
            'FaceColor',hex2rgb(qw_cols{i}),'EdgeColor','none','FaceAlpha',0.7);
        text(ax,0.5,y+0.025,sprintf('QW Layer %d',i),...
            'Color','#000000','HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
    end
    text(ax,0.5,0.66,'3-Layer VCSEL Active Region','Color',ct,...
        'HorizontalAlignment','center','FontSize',7.5,'FontWeight','bold');
    
    % RTD section
    rectangle(ax,'Position',[0.1 0.30 0.8 0.08],...
        'FaceColor',hex2rgb('#21262d'),'EdgeColor',hex2rgb(cr),'LineWidth',1.5);
    for i=1:3
        xb = 0.15 + (i-1)*0.22;
        rectangle(ax,'Position',[xb 0.31 0.05 0.06],...
            'FaceColor',hex2rgb(cr),'EdgeColor','none');
    end
    text(ax,0.5,0.34,'RTD: 3-Barrier (AlAs/InGaAs)','Color',ct,...
        'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
    
    % Bottom DBR
    for i=1:5
        y = 0.20 - i*0.022;
        rectangle(ax,'Position',[0.1 y 0.8 0.018],...
            'FaceColor',hex2rgb('#1f6feb'),'EdgeColor','none');
    end
    text(ax,0.5,0.24,'Bottom DBR','Color',ct,'HorizontalAlignment','center','FontSize',7);
    
    % Neuromorphic chip label
    rectangle(ax,'Position',[0.05 0.02 0.9 0.10],...
        'FaceColor',hex2rgb('#161b22'),'EdgeColor',hex2rgb(cn),...
        'LineWidth',2,'Curvature',0.3);
    text(ax,0.5,0.07,sprintf('100 Hz Neuromorphic Chip  [%d neurons]',NEURO.N_neurons),...
        'Color',hex2rgb(cn),'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
    
    title(ax,'Device Architecture','Color',ct,'FontWeight','bold','FontSize',8);
    xlim(ax,[0 1]); ylim(ax,[0 1.05]);
end

function rgb = hex2rgb(hex)
    hex = strrep(hex,'#','');
    rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))]/255;
end

function print_summary(VCSEL, RTD, NEURO, P_peak, fr_mean, E_spike, eta)
    fprintf('\n========= DEVICE SUMMARY =========\n');
    fprintf('VCSEL (3-QW):\n');
    fprintf('  λ = %.0f nm | Γ = %.4f | τ_p = %.2f ps\n',...
            VCSEL.lambda*1e9, VCSEL.Gamma, VCSEL.tau_p*1e12);
    fprintf('  Peak P_out = %.2f mW\n', P_peak*1e3);
    fprintf('RTD (3-Barrier):\n');
    fprintf('  PVR = %.1f | V_peak = %.2f V | V_NDR = %.2f V\n',...
            RTD.PVR, RTD.V_peak, RTD.V_bias);
    fprintf('Neuromorphic Chip:\n');
    fprintf('  f_chip = %d Hz | N_neurons = %d\n', NEURO.f_chip, NEURO.N_neurons);
    fprintf('  Mean firing rate = %.1f Hz\n', fr_mean);
    fprintf('Energy:\n');
    fprintf('  E/spike = %.2f fJ | η_wall = %.1f%%\n', E_spike*1e15, eta*100);
    fprintf('===================================\n');
end
