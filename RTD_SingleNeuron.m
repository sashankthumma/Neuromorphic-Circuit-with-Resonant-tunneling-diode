% =========================================================================
%  RTD_SingleNeuron.m
%  Module 2: Single RTD Neuron — Spike Dynamics
%
%  Implements Equations (1–5) from Owen-Newns et al. 2026
%  Gauss-Legendre 4th order integrator (as in paper) + ode15s fallback
%
%  Plots:
%   (a) Single spike response to optical pulse — V(t) and I(t)
%   (b) Multi-pulse response — refractory period demonstration
%   (c) Excitability threshold: sub- vs supra-threshold inputs
%   (d) RTD-Laser coupled system: electrical spike -> optical output
%   (e) Spike width vs regime (THz scaling)
% =========================================================================

function RTD_SingleNeuron(P)

if nargin < 1
    P = get_regime_params('THz');
end

kBTq = P.kBTq;
fV   = @(V) P.A .* log((1+exp((P.b-P.c+P.n1.*V)./kBTq))./ ...
                         (1+exp((P.b-P.c-P.n1.*V)./kBTq))).* ...
            (pi/2+atan((P.c-P.n1.*V)./P.d)) + ...
            P.H.*(exp(P.n2.*V./kBTq)-1);

% ----- Time scale from circuit parameters --------------------------------
tau_spike = 2 * pi * sqrt(P.L * P.C);   % natural oscillation period
dt        = tau_spike / 500;              % integration step
t_end     = tau_spike * 20;              % simulate 20 spike periods
t_vec     = 0 : dt : t_end;
N_steps   = length(t_vec);

V0_bias   = 0.75;    % bias in valley region
kappa     = P.kappa;

fprintf('   Spike duration ~ %s  |  Max freq ~ %s\n', ...
    format_time(tau_spike), format_freq(1/(3*tau_spike)));

% =========== (a) Single spike — Gauss-Legendre integrator ================
% 4th-order implicit Gauss-Legendre (as used in the paper)
S_pulse_amp = 600;   % optical pulse amplitude (arb units)
pulse_dur   = tau_spike * 0.5;
pulse_start = t_vec(round(N_steps*0.1));

S0_single = @(t) S_pulse_amp * (t >= pulse_start & t <= pulse_start + pulse_dur);

[V_single, I_single] = gauss_legendre_rtd(t_vec, V0_bias, S0_single, P, fV);

% =========== (b) Multi-pulse — refractory period =========================
n_pulses   = 5;
pulse_sep  = tau_spike * 3;
t_multi    = 0 : dt : pulse_sep * n_pulses * 1.5;
S0_multi   = @(t) S_pulse_amp * any(abs(t - (1:n_pulses)*pulse_sep) < pulse_dur/2);

V_multi = zeros(size(t_multi));
I_multi = zeros(size(t_multi));
V_multi(1) = V0_bias;
I_multi(1) = fV(V0_bias);
for k = 2:length(t_multi)
    t_k = t_multi(k-1);
    S0_k = S_pulse_amp * any(abs(t_k - (1:n_pulses)*pulse_sep) < pulse_dur);
    [dV, dI] = rtd_derivatives(V_multi(k-1), I_multi(k-1), V0_bias, S0_k, P, fV);
    V_multi(k) = V_multi(k-1) + dt * dV;
    I_multi(k) = I_multi(k-1) + dt * dI;
end

% =========== (c) Threshold sweep — sub vs supra ==========================
amps    = [100, 300, 450, 600, 750];   % optical amplitudes to test
colors  = lines(length(amps));
V_thresh = zeros(length(amps), N_steps);

for ia = 1:length(amps)
    S0_a = @(t) amps(ia) * (t >= pulse_start & t <= pulse_start + pulse_dur);
    V_thresh(ia,:) = gauss_legendre_rtd(t_vec, V0_bias, S0_a, P, fV);
end

% =========== (d) RTD-Laser coupled system ================================
% Couple RTD output to laser rate equations
[V_laser, I_laser, N_laser, S_laser] = solve_rtd_laser(t_vec, V0_bias, S0_single, P, fV);

% =========== (e) Regime scaling ------------------------------------------
regimes_list = {'Experiment','GHz_Sim','NearTHz','THz'};
spike_durs   = zeros(1,4);
max_freqs    = zeros(1,4);
for kr = 1:4
    Pr = get_regime_params(regimes_list{kr});
    spike_durs(kr) = 2*pi*sqrt(Pr.L*Pr.C);
    max_freqs(kr)  = 1/(3*spike_durs(kr));
end

% ===== Plotting ===========================================================
fig = figure('Name','Module 2: Single RTD Neuron','NumberTitle','off', ...
             'Position',[80 50 1300 900]);

t_scale  = get_tscale(tau_spike);
t_unit   = get_tunit(tau_spike);
t_plot   = t_vec .* t_scale;

% (a) Single spike
ax1 = subplot(3,2,1);
yyaxis left;
plot(t_plot, V_single, 'b-', 'LineWidth', 2);
ylabel('RTD Voltage (V)');
yyaxis right;
plot(t_plot, I_single.*1000, 'r--', 'LineWidth', 1.5);
ylabel('Current (mA)');
xlabel(sprintf('Time (%s)', t_unit));
title('(a) Single Spike Response');
xline(pulse_start*t_scale, 'g:', 'Pulse', 'LabelOrientation','horizontal');
grid on;

% (b) Multi-pulse
ax2 = subplot(3,2,2);
t_multi_plot = t_multi .* t_scale;
plot(t_multi_plot, V_multi, 'b-', 'LineWidth', 1.5); hold on;
for np = 1:n_pulses
    xline(np*pulse_sep*t_scale, 'g:', 'LineWidth', 1);
end
xlabel(sprintf('Time (%s)', t_unit));
ylabel('RTD Voltage (V)');
title(sprintf('(b) %d-Pulse Train — Refractory Period', n_pulses));
grid on;

% (c) Threshold sweep
ax3 = subplot(3,2,3);
hold on;
for ia = 1:length(amps)
    plot(t_plot, V_thresh(ia,:), 'Color', colors(ia,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('S_0 = %d', amps(ia)));
end
xlabel(sprintf('Time (%s)', t_unit));
ylabel('RTD Voltage (V)');
title('(c) Excitability — Sub vs Supra-Threshold');
legend('Location','northeast','FontSize',8);
yline(fV(V0_bias), 'k:', 'Bias', 'LineWidth', 1);
grid on;

% (d) RTD-Laser output
ax4 = subplot(3,2,4);
yyaxis left;
plot(t_plot, V_laser, 'b-', 'LineWidth', 2);
ylabel('RTD Voltage (V)');
yyaxis right;
plot(t_plot, S_laser, 'r-', 'LineWidth', 1.5);
ylabel('Photon Number S');
xlabel(sprintf('Time (%s)', t_unit));
title('(d) RTD–Laser Coupled System (Optical I/O)');
grid on;

% (e) Regime comparison
ax5 = subplot(3,2,5);
regime_names = {'Experiment','GHz Sim','Near-THz','THz'};
regime_colors = [0.89 0.15 0.15; 0.11 0.62 0.46; 0.60 0.40 0.80; 0.21 0.47 0.75];
b = bar(spike_durs * 1e12, 'FaceColor','flat');
for kr=1:4, b.CData(kr,:) = regime_colors(kr,:); end
set(gca,'XTickLabel', regime_names, 'YScale','log');
ylabel('Spike Duration (ps)');
title('(e) Spike Duration vs Regime');
grid on;
yline(1, 'k--', '1 ps (THz threshold)', 'LabelVerticalAlignment','bottom');

% (f) Max frequency bar
ax6 = subplot(3,2,6);
b2 = bar(max_freqs / 1e9, 'FaceColor','flat');
for kr=1:4, b2.CData(kr,:) = regime_colors(kr,:); end
set(gca,'XTickLabel', regime_names, 'YScale','log');
ylabel('Max Spike Frequency (GHz)');
title('(f) Max Spike Frequency vs Regime');
grid on;
yline(1000, 'r--', '1 THz', 'LabelVerticalAlignment','bottom');

sgtitle(sprintf('Module 2: Single RTD Neuron Dynamics | %s', P.regime_label), ...
        'FontWeight','bold','FontSize',12);
drawnow;
end


% =========================================================================
%  gauss_legendre_rtd  —  4th-order implicit Gauss-Legendre integrator
%  As used in the paper (Section 2)
% =========================================================================
function [V_out, I_out] = gauss_legendre_rtd(t_vec, V0_bias, S0_func, P, fV)
    N  = length(t_vec);
    V  = zeros(1,N);
    I  = zeros(1,N);
    V(1) = V0_bias;
    I(1) = fV(V0_bias);

    % Gauss-Legendre nodes and weights (2-point, 4th order)
    c1 = 0.5 - sqrt(3)/6;
    c2 = 0.5 + sqrt(3)/6;
    w1 = 0.5; w2 = 0.5;

    for k = 1:N-1
        dt  = t_vec(k+1) - t_vec(k);
        tk  = t_vec(k);
        S1  = S0_func(tk + c1*dt);
        S2  = S0_func(tk + c2*dt);

        % Explicit Euler predictor for initial guess
        [dV_k, dI_k] = rtd_derivatives(V(k), I(k), V0_bias, S0_func(tk), P, fV);
        V1 = V(k) + c1*dt*dV_k;
        I1 = I(k) + c1*dt*dI_k;
        V2 = V(k) + c2*dt*dV_k;
        I2 = I(k) + c2*dt*dI_k;

        % 2 Newton iterations to solve implicit stage equations
        for iter = 1:2
            [dV1, dI1] = rtd_derivatives(V1, I1, V0_bias, S1, P, fV);
            [dV2, dI2] = rtd_derivatives(V2, I2, V0_bias, S2, P, fV);
            V1 = V(k) + dt*(0.25*dV1 + (0.25-sqrt(3)/6)*dV2);
            I1 = I(k) + dt*(0.25*dI1 + (0.25-sqrt(3)/6)*dI2);
            V2 = V(k) + dt*((0.25+sqrt(3)/6)*dV1 + 0.25*dV2);
            I2 = I(k) + dt*((0.25+sqrt(3)/6)*dI1 + 0.25*dI2);
        end

        [dV1, dI1] = rtd_derivatives(V1, I1, V0_bias, S1, P, fV);
        [dV2, dI2] = rtd_derivatives(V2, I2, V0_bias, S2, P, fV);
        V(k+1) = V(k) + dt*(w1*dV1 + w2*dV2);
        I(k+1) = I(k) + dt*(w1*dI1 + w2*dI2);
    end
    V_out = V;
    I_out = I;
end


% =========================================================================
%  rtd_derivatives  —  RHS of RTD ODEs (Eqs 1-2)
% =========================================================================
function [dV, dI] = rtd_derivatives(V, I, V0, S0, P, fV)
    dV = (I - fV(V) - P.kappa * S0) / P.C;
    dI = (V0 - V - P.R * I) / P.L;
end


% =========================================================================
%  solve_rtd_laser  —  Solve RTD + laser rate equations (Eqs 1-5)
% =========================================================================
function [V_out, I_out, N_out, S_out] = solve_rtd_laser(t_vec, V0_bias, S0_func, P, fV)
    N_steps = length(t_vec);
    V_out = zeros(1,N_steps); I_out = zeros(1,N_steps);
    N_out = zeros(1,N_steps); S_out = zeros(1,N_steps);
    V_out(1) = V0_bias;
    I_out(1) = fV(V0_bias);
    N_out(1) = P.N0 * 1.1;
    S_out(1) = 0;

    for k = 1:N_steps-1
        dt = t_vec(k+1) - t_vec(k);
        V  = V_out(k); I = I_out(k);
        N  = N_out(k); S = S_out(k);
        S0k = S0_func(t_vec(k));

        [dV, dI] = rtd_derivatives(V, I, V0_bias, S0k, P, fV);
        dN = P.J_laser - (P.gamma_m + P.gamma_l + P.gamma_nr)*N ...
             - P.gamma_m*(N - P.N0)*S + P.eta*I/1.6e-19;
        dS = P.gamma_m*(N - P.N0)*S - S/P.tau_ph + P.gamma_m*N;

        V_out(k+1) = V + dt*dV;
        I_out(k+1) = max(I + dt*dI, 0);
        N_out(k+1) = max(N + dt*dN, 0);
        S_out(k+1) = max(S + dt*dS, 0);
    end
end


% =========================================================================
%  Time axis helpers
% =========================================================================
function s = get_tscale(tau)
    if tau < 1e-12, s = 1e15;
    elseif tau < 1e-9, s = 1e12;
    elseif tau < 1e-6, s = 1e9;
    else, s = 1e6; end
end

function s = get_tunit(tau)
    if tau < 1e-12, s = 'fs';
    elseif tau < 1e-9, s = 'ps';
    elseif tau < 1e-6, s = 'ns';
    else, s = 'us'; end
end

function s = format_time(t)
    if t < 1e-12,     s = sprintf('%.2f fs', t*1e15);
    elseif t < 1e-9,  s = sprintf('%.2f ps', t*1e12);
    elseif t < 1e-6,  s = sprintf('%.2f ns', t*1e9);
    else,             s = sprintf('%.2f us', t*1e6); end
end

function s = format_freq(f)
    if f >= 1e12,     s = sprintf('%.2f THz', f/1e12);
    elseif f >= 1e9,  s = sprintf('%.2f GHz', f/1e9);
    elseif f >= 1e6,  s = sprintf('%.2f MHz', f/1e6);
    else,             s = sprintf('%.1f kHz', f/1e3); end
end
