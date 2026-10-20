%% =========================================================
%  MODULE 4: Full PRL Node (Photodetector + RTD + VCSEL)
%  O/E/O architecture from Paper 1:
%  Optical in → PD → RTD spike → VCSEL → Optical out
%  Also demonstrates coincidence detection (XOR logic)
%% =========================================================

function prl_node_simulation()

%% --- System Parameters ---
% Photodetector
PD_responsivity = 0.8;    % A/W  (750-1650 nm detector, Paper 1)
PD_bandwidth    = 10e9;   % Hz

% RTD
C_rtd  = 15e-15;
Rs_rtd = 4.0;
L_rtd  = 0.5e-9;
Ibias = 0.880;  % biased at 880 mV

% VCSEL (modulation model — Paper 5, 1300 nm device)
I_th_vcsel    = 1.4e-3;   % threshold current [A]
eta_vcsel     = 0.5;      % differential efficiency [W/A]
I_vcsel_bias  = 3.5e-3;   % VCSEL bias [A]

%% --- Time setup ---
dt  = 0.05e-9;
T   = 1200e-9;
t   = 0:dt:T;
N   = length(t);
t_ns = t * 1e9;

%% --- Optical input: two channels (coincidence detection) ---
% Channel A: single pulses
% Channel B: paired pulses offset in time
% XOR: fire only when exactly one channel active

P_optical_A = zeros(1,N);   % optical power [W]
P_optical_B = zeros(1,N);

pw = 3e-9;  % pulse width

% Channel A events
evA = [100 300 500 700 900] * 1e-9;
for e = evA
    P_optical_A((t>=e)&(t<e+pw)) = 0.8e-3;  % 0.8 mW
end

% Channel B events (some overlap, some not)
evB = [200 300 600 700 1000] * 1e-9;
for e = evB
    P_optical_B((t>=e)&(t<e+pw)) = 0.8e-3;
end

% Combined optical input to PD
P_combined = P_optical_A + P_optical_B;

%% --- PD converts optical → electrical current ---
% Simple RC lowpass for bandwidth limit
tau_pd = 1/(2*pi*PD_bandwidth);
I_pd = zeros(1,N);
for i = 2:N
 alpha_pd = 1 - exp(-dt/tau_pd);
% line 63 becomes:
I_pd(i) = I_pd(i-1) + alpha_pd * (PD_responsivity*P_combined(i) - I_pd(i-1));
end

%% --- RTD spiking dynamics ---
V_rtd = zeros(1,N);
I_L   = zeros(1,N);
V_rtd(1) = 0.82;

for i = 2:N
    V   = V_rtd(i-1);
    IL  = I_L(i-1);
    Iin = I_pd(i-1);

    dV  = (Ibias + Iin - rtd_current(V) - IL) / C_rtd;
    dIL = (V - Rs_rtd*IL) / L_rtd;

    V_rtd(i) = max(0, min(V + dt*dV, 1.4));
    I_L(i)   = IL + dt*dIL;
end

%% --- VCSEL output (rate model — Paper 5) ---
% VCSEL modulated by RTD voltage above threshold
V_spike_thresh = 0.90;
I_vcsel_mod = I_vcsel_bias + max(V_rtd - V_spike_thresh, 0) * 8e-3;

% VCSEL power above threshold
P_vcsel = max(eta_vcsel * (I_vcsel_mod - I_th_vcsel), 0);

% Add 100 ps VCSEL spiking timescale (lowpass)
tau_vcsel = 100e-12;
P_vcsel_filtered = zeros(1,N);
for i = 2:N
    P_vcsel_filtered(i) = P_vcsel_filtered(i-1) + ...
        dt/tau_vcsel * (P_vcsel(i) - P_vcsel_filtered(i-1));
end

%% --- Detect spikes and compute coincidence ---
spike_thresh_W = 0.5e-3;  % detection threshold [W]
spikes_out = find(diff(P_vcsel_filtered > spike_thresh_W) == 1);

% Classify: coincidence (both A&B) vs single channel
coin_thresh = 1.4e-3;  % 2x power = both channels
is_coincidence = P_combined(spikes_out) > coin_thresh;

%% --- Plot ---
figure('Name','PRL Node — O/E/O Neuromorphic Neuron','Position',[100 100 950 750]);

subplot(5,1,1);
plot(t_ns, P_optical_A*1e3,'b-','LineWidth',1.5); hold on;
plot(t_ns, P_optical_B*1e3,'r-','LineWidth',1.5);
ylabel('Optical P (mW)','FontSize',10);
title('PRL Node: Photodetector → RTD → VCSEL (O/E/O Architecture)','FontSize',13);
legend('Channel A','Channel B','Location','northeast','FontSize',9);
grid on; xlim([0 T*1e9]);

subplot(5,1,2);
plot(t_ns, I_pd*1e3,'Color',[0.2 0.7 0.4],'LineWidth',1.5);
ylabel('I_{PD} (mA)','FontSize',10);
title('Photodetector output current','FontSize',10);
grid on; xlim([0 T*1e9]);

subplot(5,1,3);
plot(t_ns, V_rtd,'b-','LineWidth',1.8); hold on;
yline(V_spike_thresh,'--r','Spike threshold','LineWidth',1.2,...
      'LabelHorizontalAlignment','left','FontSize',8);
ylabel('V_{RTD} (V)','FontSize',10);
title('RTD voltage (spiking dynamics)','FontSize',10);
ylim([0.5 1.3]); grid on; xlim([0 T*1e9]);

subplot(5,1,4);
plot(t_ns, P_vcsel_filtered*1e3,'Color',[0.8 0.2 0.8],'LineWidth',1.8); hold on;
if ~isempty(spikes_out)
    plot(t_ns(spikes_out(is_coincidence)), ...
         P_vcsel_filtered(spikes_out(is_coincidence))*1e3, ...
         'ko','MarkerSize',9,'MarkerFaceColor','k','DisplayName','Coincidence');
    plot(t_ns(spikes_out(~is_coincidence)), ...
         P_vcsel_filtered(spikes_out(~is_coincidence))*1e3, ...
         'r^','MarkerSize',8,'MarkerFaceColor','r','DisplayName','Single channel');
end
ylabel('VCSEL P_{out} (mW)','FontSize',10);
title('VCSEL optical output spikes','FontSize',10);
legend('VCSEL output','Coincidence spike','Single-ch spike',...
       'Location','northeast','FontSize',8);
grid on; xlim([0 T*1e9]);

subplot(5,1,5);
% XOR truth: coincidence=0 (both active, suppress), single=1 (fire)
xor_signal = zeros(1,N);
for s = 1:length(spikes_out)
    idx = spikes_out(s);
    rng_idx = max(1,idx-20):min(N,idx+40);
    xor_signal(rng_idx) = ~is_coincidence(s);
end
plot(t_ns, xor_signal, 'Color',[0.8 0.5 0.0],'LineWidth',2.0);
ylabel('XOR output','FontSize',10);
xlabel('Time (ns)','FontSize',12);
title('Logical XOR output (single-channel only)','FontSize',10);
ylim([-0.1 1.4]); grid on; xlim([0 T*1e9]);

fprintf('   Total output spikes: %d\n', length(spikes_out));
if ~isempty(spikes_out)
    fprintf('   Coincidence spikes: %d\n', sum(is_coincidence));
    fprintf('   Single-channel spikes: %d\n', sum(~is_coincidence));
end

end

%% --- RTD I-V helper ---
function I = rtd_current(V)
    Ip = 7.5e-3; Vp = 0.55; Iv = 1.2e-3; Vv = 0.90;
    I_res = Ip.*exp(-((V-Vp).^2)/(2*0.06^2));
    I_val = Iv + 3.5e-3.*(V-Vv).^2.*(V>Vv);
    I_bg  = 0.3e-3.*V + 0.8e-3.*exp(3.5*(V-1.1)).*(V>1.0);
    I = max(I_res + I_val + I_bg, 0.05e-3);
end
