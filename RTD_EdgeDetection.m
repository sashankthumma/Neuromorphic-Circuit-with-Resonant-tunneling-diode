% =========================================================================
%  RTD_EdgeDetection.m
%  Module 3: Rising Edge Detection in Mackey-Glass Time-Series
%
%  Implements Section 3.1 of Owen-Newns et al. 2026
%  Dual modulation: optical input = y(t), electronic bias = y(t-tau)
%  Extended to THz-regime parameters
%
%  Plots:
%   (a) Mackey-Glass time-series
%   (b) Delayed difference D_tau(t)
%   (c) Optical input encoding
%   (d) Bias voltage modulation
%   (e) RTD spike output
%   (f) Detection accuracy vs delay parameter
% =========================================================================

function RTD_EdgeDetection(P)

if nargin < 1
    P = get_regime_params('THz');
end

kBTq = P.kBTq;
fV   = @(V) P.A .* log((1+exp((P.b-P.c+P.n1.*V)./kBTq))./ ...
                         (1+exp((P.b-P.c-P.n1.*V)./kBTq))).* ...
            (pi/2+atan((P.c-P.n1.*V)./P.d)) + ...
            P.H.*(exp(P.n2.*V./kBTq)-1);

% =========== Generate Mackey-Glass time-series ===========================
% MG equation: dx/dt = 0.2*x(t-tau_MG)/(1+x(t-tau_MG)^10) - 0.1*x(t)
% tau_MG=17, step=0.01, sampled every 1s — classic chaotic regime

N_MG     = 1200;    % total time steps to generate
tau_MG   = 17;      % MG delay parameter (controls chaos level)
dt_MG    = 0.01;    % integration step for MG
sample_n = 100;     % sample every 100 steps
x_mg     = mackey_glass(N_MG * sample_n + tau_MG*100, dt_MG, tau_MG);
x_mg     = x_mg(1 : sample_n : end);
x_mg     = x_mg(1:1000);  % use 1000 timesteps
N_ts     = length(x_mg);

% Normalise to [0, 1]
x_mg = (x_mg - min(x_mg)) ./ (max(x_mg) - min(x_mg));

% =========== Dual modulation encoding ====================================
delay_ts   = 10;               % delay in timesteps
S_max      = 600;              % max optical amplitude
V0_nominal = 0.75;             % nominal bias voltage
dV_max     = 0.20;             % bias modulation amplitude

% y(t) -> optical amplitude
S_opt = S_max .* x_mg;

% y(t-delay) -> bias voltage modulation
x_delayed = [x_mg(1)*ones(1, delay_ts), x_mg(1:end-delay_ts)];
V_bias_mod = V0_nominal + dV_max .* x_delayed;

% Delayed difference (what the RTD 'sees')
D_tau = x_mg - x_delayed;

% =========== Threshold for edge detection ================================
threshold  = 0.45;
target_idx = find(D_tau > threshold);  % ground-truth rising edge steps

% =========== Run RTD simulation for each timestep ========================
tau_spike   = 2 * pi * sqrt(P.L * P.C);
pulse_dur   = tau_spike * 0.5;
dt_sim      = tau_spike / 300;

spike_out   = zeros(1, N_ts);   % 1 = spike detected, 0 = quiescent

fprintf('   Running %d timestep edge detection simulation...\n', N_ts);

for ts = 1:N_ts
    % Each timestep: simulate one optical pulse period
    t_local = 0 : dt_sim : tau_spike * 4;
    S0_amp  = S_opt(ts);
    Vb      = V_bias_mod(ts);
    S0_func = @(t) S0_amp .* (t >= tau_spike*0.5 & t <= tau_spike*0.5 + pulse_dur);

    V_ts = zeros(size(t_local));
    I_ts = zeros(size(t_local));
    V_ts(1) = Vb;
    I_ts(1) = fV(Vb);

    spiked = false;
    for k = 2:length(t_local)
        S0k = S0_func(t_local(k-1));
        dV  = (I_ts(k-1) - fV(V_ts(k-1)) - P.kappa*S0k) / P.C;
        dI  = (Vb - V_ts(k-1) - P.R*I_ts(k-1)) / P.L;
        V_ts(k) = V_ts(k-1) + dt_sim * dV;
        I_ts(k) = I_ts(k-1) + dt_sim * dI;
        if V_ts(k) > Vb + 0.2 && ~spiked
            spiked = true;
        end
    end
    spike_out(ts) = double(spiked);
end

% =========== Detection metrics ===========================================
detected     = find(spike_out == 1);
true_pos     = intersect(detected, target_idx);
false_pos    = setdiff(detected, target_idx);
false_neg    = setdiff(target_idx, detected);
precision    = length(true_pos) / max(length(detected), 1);
recall       = length(true_pos) / max(length(target_idx), 1);
F1           = 2 * precision * recall / max(precision + recall, 1e-6);

fprintf('   Edge detection: TP=%d, FP=%d, FN=%d | P=%.2f R=%.2f F1=%.2f\n', ...
    length(true_pos), length(false_pos), length(false_neg), precision, recall, F1);

% =========== Delay sweep — accuracy vs delay parameter ===================
delay_range   = 1:20;
accuracy_vec  = zeros(size(delay_range));
for di = 1:length(delay_range)
    d = delay_range(di);
    x_del_tmp = [x_mg(1)*ones(1,d), x_mg(1:end-d)];
    D_tmp = x_mg - x_del_tmp;
    tgt_tmp = find(D_tmp > threshold);
    % Simplified: accuracy = fraction of targets detectable given threshold
    detectability = sum(D_tmp(tgt_tmp) > threshold) / max(length(find(D_tmp > threshold*0.8)), 1);
    accuracy_vec(di) = min(detectability, 1);
end

% ===== Plotting ===========================================================
fig = figure('Name','Module 3: Edge Detection','NumberTitle','off', ...
             'Position',[100 50 1300 950]);

t_axis = 1:N_ts;

% (a) MG signal
subplot(3,2,1);
plot(t_axis, x_mg, 'b-', 'LineWidth', 1.5);
hold on;
for idx = target_idx'
    xline(idx, 'r:', 'Alpha', 0.4);
end
ylabel('Amplitude (norm)');
title('(a) Mackey-Glass Time-Series');
xlim([1 N_ts]); grid on;
legend({'MG signal','Rising edges'},'FontSize',8,'Location','northwest');

% (b) Delayed difference
subplot(3,2,2);
plot(t_axis, D_tau, 'Color', [0.8 0.4 0], 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', sprintf('Threshold = %.2f', threshold), 'LineWidth', 1.5);
fill([t_axis, fliplr(t_axis)], ...
     [max(D_tau,threshold), threshold*ones(1,N_ts)], ...
     'r', 'FaceAlpha', 0.15, 'EdgeColor','none');
ylabel('\Delta_\tau(t)');
title(sprintf('(b) Delayed Difference (delay = %d steps)', delay_ts));
xlim([1 N_ts]); grid on;

% (c) Optical input encoding
subplot(3,2,3);
plot(t_axis, S_opt, 'Color',[0.85 0.33 0.1], 'LineWidth', 1.2);
ylabel('Optical Amplitude (arb)');
title('(c) Optical Input S_0(t) — encodes y(t)');
xlim([1 N_ts]); grid on; ylim([0 S_max*1.1]);

% (d) Bias voltage modulation
subplot(3,2,4);
plot(t_axis, V_bias_mod, 'Color',[0.13 0.55 0.13], 'LineWidth', 1.5);
yline(V0_nominal, 'k:', 'Nominal bias', 'LineWidth', 1);
ylabel('Bias Voltage (V)');
title('(d) Bias Voltage V_0(t) — encodes y(t-\tau)');
xlim([1 N_ts]); grid on;

% (e) RTD spike output
subplot(3,2,5);
stem(t_axis, spike_out, 'b', 'MarkerFaceColor','b', 'MarkerSize', 4);
hold on;
for idx = target_idx'
    patch([idx-0.5 idx+0.5 idx+0.5 idx-0.5], [0 0 1.15 1.15], ...
          'r', 'FaceAlpha', 0.2, 'EdgeColor','none');
end
stem(t_axis(logical(spike_out & ~ismember(t_axis, target_idx))), ...
     spike_out(logical(spike_out & ~ismember(t_axis, target_idx))), ...
     'r', 'MarkerFaceColor','r', 'MarkerSize',5);
ylabel('Spike (0/1)');
xlabel('Timestep');
title(sprintf('(e) RTD Spike Output | F1=%.2f | Regime: %s', F1, P.regime_label));
xlim([1 N_ts]); ylim([0 1.3]); grid on;
legend({'Correct spike','Target window','False positive'},'FontSize',7,'Location','northeast');

% (f) Accuracy vs delay
subplot(3,2,6);
plot(delay_range, accuracy_vec * 100, 'b-o', 'LineWidth', 2, 'MarkerSize', 5);
xline(delay_ts, 'r--', sprintf('Used (%d)', delay_ts), 'LineWidth', 1.5);
xlabel('Delay \tau (timesteps)');
ylabel('Detection Score (%)');
title('(f) Edge Detection vs Delay Parameter');
grid on; ylim([0 110]);

sgtitle(sprintf('Module 3: Mackey-Glass Edge Detection | %s', P.regime_label), ...
        'FontWeight','bold','FontSize',12);
drawnow;
end


% =========================================================================
%  mackey_glass — Generate Mackey-Glass chaotic time-series
%  dx/dt = 0.2*x(t-tau)/(1+x(t-tau)^10) - 0.1*x(t)
% =========================================================================
function x = mackey_glass(N, dt, tau_mg)
    delay_steps = round(tau_mg / dt);
    x = zeros(1, N + delay_steps);
    x(1:delay_steps) = 0.5 + 0.1*randn(1, delay_steps);  % random init
    for k = delay_steps+1 : N+delay_steps
        x_d = x(k - delay_steps);
        dxdt = 0.2 * x_d / (1 + x_d^10) - 0.1 * x(k-1);
        x(k) = x(k-1) + dt * dxdt;
    end
    x = x(delay_steps+1:end);
end
