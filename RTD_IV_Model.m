% =========================================================================
%  RTD_IV_Model.m
%  Module 1: RTD I-V Curve and NDR Characterisation
%
%  Plots:
%   (a) Full I-V curve with NDR region highlighted
%   (b) f(V) conductivity function components
%   (c) NDR parameter table across all regimes
%   (d) Threshold map: optical amplitude vs bias voltage
% =========================================================================

function RTD_IV_Model(P)

if nargin < 1
    P = get_regime_params('THz');
end

% Helper functions defined in get_regime_params.m
% (MATLAB requires them in same file or on path — inline here for safety)

kBTq = P.kBTq;

fV = @(V) P.A .* log((1 + exp((P.b - P.c + P.n1.*V)./kBTq)) ./ ...
                      (1 + exp((P.b - P.c - P.n1.*V)./kBTq))) .* ...
          (pi/2 + atan((P.c - P.n1.*V) ./ P.d)) + ...
          P.H .* (exp(P.n2.*V./kBTq) - 1);

V_range = linspace(0, 1.2, 2000);
I_range = fV(V_range) .* 1000;   % convert to mA

% ----- Find NDR region ---------------------------------------------------
[Ipk, ipk] = max(fV(V_range));
sub = fV(V_range(ipk:end));
[Ival, ival_rel] = min(sub);
ival = ipk + ival_rel - 1;
Vpk  = V_range(ipk);
Vval = V_range(ival);
PVCR = Ipk / max(Ival, 1e-12);

fprintf('   I-V: V_peak=%.3fV, V_valley=%.3fV, NDR=%.0f mV, PVCR=%.2f\n', ...
    Vpk, Vval, (Vval-Vpk)*1000, PVCR);

fig = figure('Name','Module 1: RTD I-V Model','NumberTitle','off', ...
             'Position',[50 50 1200 800]);

% ---- (a) Full I-V curve --------------------------------------------------
ax1 = subplot(2,2,1);
plot(V_range, I_range, 'b-', 'LineWidth', 2); hold on;

% Shade NDR region
V_ndr = V_range(ipk:ival);
I_ndr = I_range(ipk:ival);
fill([V_ndr, fliplr(V_ndr)], [I_ndr, min(I_ndr)*ones(1,length(V_ndr))], ...
     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Mark peak and valley
plot(Vpk,  Ipk*1000,  'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(Vval, Ival*1000, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Mark bias point
V0 = 0.75;
plot(V0, fV(V0)*1000, 'g*', 'MarkerSize', 12, 'LineWidth', 2);

xlabel('Bias Voltage (V)');
ylabel('Current (mA)');
title('RTD I–V Curve — Schulman Model');
legend({'I–V curve', 'NDR region', 'Peak', 'Valley', 'Bias point (0.75V)'}, ...
       'Location', 'northwest', 'FontSize', 8);
grid on;
text(Vpk+0.02, Ipk*1000*1.05, sprintf('PVCR = %.2f', PVCR), ...
     'FontSize', 9, 'Color', 'r');

% ---- (b) I-V components breakdown ----------------------------------------
ax2 = subplot(2,2,2);
term1 = P.A .* log((1 + exp((P.b - P.c + P.n1.*V_range)./kBTq)) ./ ...
                    (1 + exp((P.b - P.c - P.n1.*V_range)./kBTq)));
term2 = (pi/2 + atan((P.c - P.n1.*V_range) ./ P.d));
term3 = P.H .* (exp(P.n2.*V_range./kBTq) - 1);

plot(V_range, term1.*1000, 'b--', 'LineWidth', 1.5); hold on;
plot(V_range, (term1.*term2).*1000, 'm-', 'LineWidth', 1.5);
plot(V_range, term3.*1000, 'g-.', 'LineWidth', 1.5);
plot(V_range, I_range, 'k-', 'LineWidth', 2);
xlabel('Voltage (V)'); ylabel('Current (mA)');
title('f(V) Component Breakdown');
legend({'Term1 (tunnelling log)', 'Term1×Term2 (NDR core)', ...
        'Term3 (thermal/leakage)', 'Total f(V)'}, 'FontSize', 8);
grid on; ylim([-0.5 2]);

% ---- (c) Threshold map: spike firing regions ----------------------------
ax3 = subplot(2,2,3);
V_bias_vec  = linspace(0.70, 1.05, 60);
S_amp_vec   = linspace(0, 800, 60);
spike_map   = zeros(length(S_amp_vec), length(V_bias_vec));

for iv = 1:length(V_bias_vec)
    for is = 1:length(S_amp_vec)
        % Spike fires if optical perturbation pushes V below NDR entry
        % Simplified threshold: I_optical > I_peak - I(V_bias)
        I_bias = fV(V_bias_vec(iv));
        I_opt  = P.kappa * S_amp_vec(is);
        % RTD spikes if net current exceeds peak by pushing into NDR
        delta  = I_bias - fV(Vpk) + I_opt;
        spike_map(is, iv) = (delta < 0) && (V_bias_vec(iv) <= Vval + 0.05);
    end
end

imagesc(V_bias_vec, S_amp_vec, spike_map);
colormap(ax3, [0.3 0.3 0.7; 1 0.9 0.4]);
xlabel('RTD Bias Voltage (V)');
ylabel('Optical Input Amplitude (arb)');
title('Spike Firing Threshold Map');
cb = colorbar; cb.Ticks = [0.25 0.75];
cb.TickLabels = {'Non-spiking','Spiking'};
grid on;

% ---- (d) NDR sensitivity to barrier thickness (THz design) --------------
ax4 = subplot(2,2,4);
% Vary n1 parameter (proxy for barrier thickness — thinner barrier = higher n1)
n1_vec = linspace(0.10, 0.35, 200);
PVCR_vec = zeros(size(n1_vec));
NDR_width_vec = zeros(size(n1_vec));

for k = 1:length(n1_vec)
    P_tmp = P; P_tmp.n1 = n1_vec(k);
    fV_tmp = @(V) P_tmp.A .* log((1+exp((P_tmp.b-P_tmp.c+P_tmp.n1.*V)./kBTq))./ ...
                                   (1+exp((P_tmp.b-P_tmp.c-P_tmp.n1.*V)./kBTq))).* ...
                  (pi/2+atan((P_tmp.c-P_tmp.n1.*V)./P_tmp.d)) + ...
                  P_tmp.H.*(exp(P_tmp.n2.*V./kBTq)-1);
    IV_tmp = fV_tmp(V_range);
    [pk, ip] = max(IV_tmp);
    sub_tmp = IV_tmp(ip:end);
    [vl, iv_rel] = min(sub_tmp);
    if pk > 0 && vl > 0 && ip+iv_rel < length(IV_tmp)
        PVCR_vec(k) = pk / max(vl, 1e-12);
        NDR_width_vec(k) = (V_range(ip+iv_rel-1) - V_range(ip)) * 1000;
    end
end

yyaxis left;
plot(n1_vec, PVCR_vec, 'b-', 'LineWidth', 2);
ylabel('PVCR'); ylim([0 10]);
yyaxis right;
plot(n1_vec, NDR_width_vec, 'r--', 'LineWidth', 2);
ylabel('NDR width (mV)');
xlabel('n_1 parameter (barrier thickness proxy)');
title('PVCR and NDR Width vs Barrier Design');
xline(P.n1, 'k:', 'LineWidth', 1.5, 'Label', 'Paper value');
legend({'PVCR','NDR width'}, 'Location', 'northwest', 'FontSize', 8);
grid on;

sgtitle(sprintf('Module 1: RTD I–V Model | Regime: %s', P.regime_label), ...
        'FontWeight', 'bold', 'FontSize', 12);

drawnow;
end
