%% =========================================================
%  MODULE 3: Excitability Threshold Sweep
%  Replicates Paper 1, Figs 1(d-f):
%  RTD bias = 860, 870, 880 mV → threshold shifts
%  Shows all-or-nothing spiking with adjustable threshold
%% =========================================================

function rtd_threshold_sweep()

%% --- Sweep parameters (from Paper 1) ---
bias_voltages = [0.860, 0.870, 0.880];   % [V]
bias_labels   = {'860 mV', '870 mV', '880 mV'};
colors        = {[0.1 0.4 0.9], [0.9 0.5 0.1], [0.8 0.1 0.1]};

dt  = 0.05e-9;
T   = 700e-9;
t   = 0:dt:T;
N   = length(t);
t_us = t * 1e6;   % convert to microseconds (Paper 1 x-axis)

% Input: graded pulse train — 7 pulses of increasing amplitude
n_pulses = 7;
pulse_start = linspace(50e-9, 600e-9, n_pulses);
pulse_width = 3e-9;
pulse_amps  = linspace(0.05e-3, 0.75e-3, n_pulses);  % increasing

I_input = zeros(1, N);
for k = 1:n_pulses
    mask = (t >= pulse_start(k)) & (t < pulse_start(k)+pulse_width);
    I_input(mask) = pulse_amps(k);
end

%% --- Simulate for each bias ---
figure('Name','Excitability Threshold Sweep','Position',[100 100 900 700]);

for b = 1:3
    Vbias = bias_voltages(b);

    % Convert bias voltage to equivalent bias current
    Ibias = rtd_current(Vbias);

    % State variables
    C  = 15e-15; Rs = 4.0; L = 0.5e-9;
    V_rtd = zeros(1,N);
    I_L   = zeros(1,N);
    V_rtd(1) = Vbias - 0.05;
    I_L(1)   = 0;

    for i = 2:N
        V   = V_rtd(i-1);
        IL  = I_L(i-1);
        Iin = I_input(i-1);

        dV  = (Ibias + Iin - rtd_current(V) - IL) / C;
        dIL = (V - Rs*IL) / L;

        V_rtd(i) = max(0, min(V + dt*dV, 1.4));
        I_L(i)   = IL + dt*dIL;
    end

    % Normalize to optical output (VCSEL modulation — Paper 1)
    % Optical output intensity proportional to V_rtd deviation
    V_spike_thresh = 0.92;
    optical_out = max(V_rtd - V_spike_thresh, 0) * 60;  % scaled to mV (Fig 1d-f)

    subplot(4,1,b+1);
    plot(t_us, optical_out*1e3, 'Color', colors{b}, 'LineWidth', 1.8);
    hold on;

    % Mark spikes
    spikes = find(diff(optical_out > 0.002) == 1);
    if ~isempty(spikes)
        plot(t_us(spikes), optical_out(spikes)*1e3, 'kv', ...
             'MarkerSize',7,'MarkerFaceColor','k');
    end

    ylabel('Optical out (mV)','FontSize',10);
    text(0.62, max(optical_out*1e3)*0.65, ...
         sprintf('U_S = %s | Spikes: %d', bias_labels{b}, length(spikes)), ...
         'FontSize',10,'Color',colors{b},'FontWeight','bold', ...
         'Units','normalized');
    grid on; xlim([0 T*1e6]); ylim([-2 60]);
    set(gca,'FontSize',10);
end

% Input panel (top)
subplot(4,1,1);
stairs(t_us, I_input*1e3, 'k-', 'LineWidth', 1.5);
ylabel('I_{input} (mA)','FontSize',10);
title('Excitability Threshold Sweep — Bias Voltage Effect (Paper 1, Fig 1d-f)', ...
      'FontSize',13);
grid on; xlim([0 T*1e6]);
text(0.02, 0.75, 'Graded input pulses (increasing amplitude →)', ...
     'FontSize',9,'Units','normalized','Color',[0.4 0.4 0.4]);

xlabel('Time (\mus)','FontSize',12);

fprintf('   Bias sweep complete: 860, 870, 880 mV\n');
fprintf('   Closer to NDR peak → lower threshold → more spikes\n');

end

%% --- RTD I-V helper ---
function I = rtd_current(V)
    Ip = 7.5e-3; Vp = 0.55;
    Iv = 1.2e-3; Vv = 0.90;
    I_res = Ip .* exp(-((V-Vp).^2)/(2*0.06^2));
    I_val = Iv + 3.5e-3.*(V-Vv).^2.*(V>Vv);
    I_bg  = 0.3e-3.*V + 0.8e-3.*exp(3.5*(V-1.1)).*(V>1.0);
    I = max(I_res + I_val + I_bg, 0.05e-3);
end
