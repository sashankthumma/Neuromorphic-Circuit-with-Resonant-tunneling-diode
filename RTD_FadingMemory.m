% =========================================================================
%  RTD_FadingMemory.m
%  Module 5: Neuromorphic Photonic Fading Memory Cell
%
%  Implements Section 3.2.2 of Owen-Newns et al. 2026
%  10 coupled RTD-laser neurons, all-to-all via attenuated feedback loop
%  Bias voltage formula: Bn = 0.9 - 0.185*n^2
%  Tunable memory depth via optical attenuation
%
%  Plots:
%   (a) Autaptic single-neuron memory (self-feedback loop)
%   (b) 10-neuron fading memory — 0% attenuation (permanent)
%   (c) 10-neuron fading memory — 10% attenuation
%   (d) 10-neuron fading memory — 20% attenuation
%   (e) Memory depth vs attenuation curve
%   (f) Per-cycle heatmap
% =========================================================================

function RTD_FadingMemory(P)

if nargin < 1
    P = get_regime_params('THz');
end

kBTq = P.kBTq;
fV   = @(V) P.A .* log((1+exp((P.b-P.c+P.n1.*V)./kBTq))./ ...
                         (1+exp((P.b-P.c-P.n1.*V)./kBTq))).* ...
            (pi/2+atan((P.c-P.n1.*V)./P.d)) + ...
            P.H.*(exp(P.n2.*V./kBTq)-1);

% =========== Parameters ==================================================
N_cells    = 10;    % number of RTD-laser neurons in memory cell
N_cycles   = 10;    % number of memory cycles to simulate
tau_spike  = 2 * pi * sqrt(P.L * P.C);
tau_delay  = tau_spike * 5;    % feedback delay = 5 spike periods
dt         = tau_spike / 200;

% Bias voltages: parabolic distribution (Eq from paper: Bn = 0.9 - 0.185*n^2)
n_vec    = (1:N_cells);
V_biases = 0.9 - 0.185 .* (n_vec ./ N_cells).^2;
% Clamp to reasonable range
V_biases = max(min(V_biases, 1.05), 0.75);

fprintf('   Bias voltages: %.3f to %.3f V\n', min(V_biases), max(V_biases));
fprintf('   Delay: %s  |  Cycle duration: %s\n', ...
    format_time(tau_delay), format_time(tau_delay));

attenuations = [0, 10, 20, 40, 60, 80, 90, 95];

% =========== (a) Autaptic single-neuron memory ===========================
fprintf('   Simulating autaptic memory...\n');
t_auto    = 0 : dt : tau_delay * N_cycles * 1.2;
N_auto    = length(t_auto);
V_auto    = zeros(1, N_auto);
I_auto    = zeros(1, N_auto);
S_fb_auto = zeros(1, N_auto);   % feedback signal (delayed laser output)

V_auto(1) = V_biases(1);
I_auto(1) = fV(V_biases(1));

delay_steps_auto = round(tau_delay / dt);
S_laser_auto     = zeros(1, N_auto);  % laser output history

% Initial trigger at t=10*dt
for k = 2:N_auto
    % Feedback: delayed laser output
    fb_idx = max(1, k - delay_steps_auto);
    S_in   = S_laser_auto(fb_idx) + ...
             (k <= 15 ? 600 : 0);   % external trigger for first few steps

    dV = (I_auto(k-1) - fV(V_auto(k-1)) - P.kappa * S_in) / P.C;
    dI = (V_biases(1) - V_auto(k-1) - P.R * I_auto(k-1)) / P.L;
    V_auto(k) = V_auto(k-1) + dt * dV;
    I_auto(k) = I_auto(k-1) + dt * dI;

    % Simple laser model: output proportional to RTD voltage spike
    if V_auto(k) > V_biases(1) + 0.15
        S_laser_auto(k) = 500;
    else
        S_laser_auto(k) = 0;
    end
end

% =========== (b-d) Fading memory at different attenuations ===============
fprintf('   Simulating fading memory (3 attenuation levels)...\n');
att_to_plot  = [0, 10, 20];
V_all        = cell(1,3);
cycle_times  = zeros(1, N_cycles+1);
for cy = 0:N_cycles
    cycle_times(cy+1) = cy * tau_delay;
end

for ai = 1:3
    att    = att_to_plot(ai) / 100;
    [V_mem, ~] = simulate_memory_cell(V_biases, N_cells, N_cycles, ...
                                       tau_delay, dt, att, P, fV);
    V_all{ai} = V_mem;
end

% =========== (e) Memory depth vs attenuation =============================
fprintf('   Computing memory depth vs attenuation...\n');
mem_depth = zeros(size(attenuations));
for ai = 1:length(attenuations)
    att = attenuations(ai) / 100;
    [~, depth] = simulate_memory_cell(V_biases, N_cells, N_cycles, ...
                                       tau_delay, dt, att, P, fV);
    mem_depth(ai) = depth;
end

% ===== Plotting ===========================================================
fig = figure('Name','Module 5: Fading Memory','NumberTitle','off', ...
             'Position',[140 50 1300 950]);

t_scale = get_tscale(tau_delay);
t_unit  = get_tunit(tau_delay);

% (a) Autaptic memory
subplot(3,3,1);
t_auto_plot = (0:N_auto-1) * dt * t_scale;
plot(t_auto_plot, V_auto, 'b-', 'LineWidth', 1.5);
for cy = 1:N_cycles
    xline(cy*tau_delay*t_scale, 'r:', 'Alpha',0.5);
end
xlabel(sprintf('Time (%s)',t_unit)); ylabel('RTD Voltage (V)');
title('(a) Autaptic Single-Neuron Memory');
ylim([0.5 1.3]); grid on;
text(0.05, 0.95, sprintf('\\tau = %s', format_time(tau_delay)), ...
     'Units','normalized','FontSize',8);

% (b-d) Fading memory traces
att_labels = {'(b) 0% Attenuation (Permanent)', ...
              '(c) 10% Attenuation', ...
              '(d) 20% Attenuation'};
cmap_cells = lines(N_cells);

for ai = 1:3
    subplot(3,3, ai+1);
    V_mem = V_all{ai};
    N_t   = size(V_mem, 2);
    t_mem_plot = (0:N_t-1) * dt * t_scale;
    hold on;
    for nc = 1:N_cells
        plot(t_mem_plot, V_mem(nc,:), 'Color', cmap_cells(nc,:), 'LineWidth', 1.2);
    end
    for cy = 1:N_cycles
        xline(cy*tau_delay*t_scale, 'k:', 'Alpha', 0.4);
    end
    xlabel(sprintf('Time (%s)',t_unit)); ylabel('RTD Voltage (V)');
    title(att_labels{ai});
    ylim([0.5 1.3]); grid on;
    xline(tau_delay*0.2*t_scale, 'g--', 'Trigger', 'LineWidth', 1.5);
end

% (e) Per-cycle spike heatmap — show which neurons fire each cycle
subplot(3,3,[5 6]);
att_heat   = 10/100;
spike_heat = zeros(N_cells, N_cycles);
V_heat     = V_all{2};   % 10% attenuation
N_t_heat   = size(V_heat, 2);

for cy = 1:N_cycles
    t_start = round((cy-0.9)*tau_delay/dt);
    t_end   = round(cy*tau_delay/dt);
    t_start = max(1, t_start);
    t_end   = min(N_t_heat, t_end);
    for nc = 1:N_cells
        window = V_heat(nc, t_start:t_end);
        spike_heat(nc, cy) = any(window > V_biases(nc) + 0.15);
    end
end

imagesc(1:N_cycles, 1:N_cells, spike_heat);
colormap(gca, [0.95 0.95 0.95; 0.2 0.6 0.9]);
colorbar('Ticks',[0.25 0.75],'TickLabels',{'Quiescent','Spiking'});
xlabel('Memory cycle'); ylabel('RTD neuron');
title('(e) Spike Activity per Cycle (10% attenuation)');
set(gca,'XTick',1:N_cycles,'YTick',1:N_cells);
grid on;

% (f) Memory depth vs attenuation
subplot(3,3,[7 8 9]);
plot(attenuations, mem_depth, 'b-o', 'LineWidth', 2.5, ...
     'MarkerFaceColor','b', 'MarkerSize', 8);
hold on;
fill([attenuations, fliplr(attenuations)], ...
     [mem_depth, zeros(1,length(attenuations))], ...
     'b', 'FaceAlpha', 0.1, 'EdgeColor','none');
xlabel('Optical Attenuation (%)');
ylabel('Memory Depth (cycles)');
title('(f) Memory Depth vs Attenuation — Tunable Fading Memory');
xlim([0 100]); ylim([0 N_cycles+1]);
grid on;
yline(N_cycles, 'r--', 'Permanent memory', 'LineWidth', 1.5, ...
      'LabelVerticalAlignment','bottom');
yline(0, 'k--', 'Zero memory', 'LineWidth', 1);
text(5, N_cycles*0.5, sprintf('THz regime: \\tau_{delay} = %s', ...
     format_time(tau_delay)), 'FontSize', 9);
legend({'Memory depth','Area'},'Location','northeast','FontSize',8);

sgtitle(sprintf('Module 5: Neuromorphic Fading Memory Cell | %s', P.regime_label), ...
        'FontWeight','bold','FontSize',12);
drawnow;
end


% =========================================================================
%  simulate_memory_cell — run 10-neuron all-to-all memory cell
% =========================================================================
function [V_all, memory_depth] = simulate_memory_cell(V_biases, N_cells, ...
                                   N_cycles, tau_delay, dt, attenuation, P, fV)
    N_t = round(N_cycles * tau_delay / dt) + 1;
    V_all = zeros(N_cells, N_t);
    I_all = zeros(N_cells, N_t);
    S_laser_all = zeros(N_cells, N_t);

    delay_steps = round(tau_delay / dt);
    trigger_dur = round(0.1 * tau_delay / dt);

    for nc = 1:N_cells
        V_all(nc,1) = V_biases(nc);
        I_all(nc,1) = fV(V_biases(nc));
    end

    for k = 2:N_t
        % All-to-all feedback: mean of all delayed laser outputs
        fb_idx = max(1, k - delay_steps);
        mean_fb = mean(S_laser_all(:, fb_idx)) * (1 - attenuation);

        % External trigger: first 10% of first cycle
        ext_trigger = (k <= trigger_dur) * 600;
        S_total_in  = mean_fb + ext_trigger;

        for nc = 1:N_cells
            S_in = S_total_in;
            dV = (I_all(nc,k-1) - fV(V_all(nc,k-1)) - P.kappa*S_in) / P.C;
            dI = (V_biases(nc)  - V_all(nc,k-1) - P.R*I_all(nc,k-1)) / P.L;
            V_all(nc,k) = V_all(nc,k-1) + dt*dV;
            I_all(nc,k) = I_all(nc,k-1) + dt*dI;

            if V_all(nc,k) > V_biases(nc) + 0.15
                S_laser_all(nc,k) = 400 + nc*20;
            else
                S_laser_all(nc,k) = 0;
            end
        end
    end

    % Count memory depth: how many cycles have at least one spiking neuron
    memory_depth = 0;
    for cy = 2:N_cycles
        t_start = round((cy-0.9)*tau_delay/dt);
        t_end   = round(cy*tau_delay/dt);
        t_start = max(1,t_start); t_end = min(N_t,t_end);
        any_spike = false;
        for nc = 1:N_cells
            if any(V_all(nc,t_start:t_end) > V_biases(nc) + 0.15)
                any_spike = true; break;
            end
        end
        if any_spike, memory_depth = memory_depth + 1;
        else, break; end
    end
end


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
    if t<1e-12, s=sprintf('%.2f fs',t*1e15);
    elseif t<1e-9, s=sprintf('%.2f ps',t*1e12);
    elseif t<1e-6, s=sprintf('%.2f ns',t*1e9);
    else, s=sprintf('%.2f us',t*1e6); end
end
