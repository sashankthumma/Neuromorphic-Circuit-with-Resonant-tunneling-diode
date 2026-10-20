% =========================================================================
%  RTD_RegimeComparison.m
%  Module 6: Regime Comparison — Experiment / GHz / Near-THz / THz
%
%  Compares spike dynamics across all four parameter regimes
%  Shows the path from current experimental reality to THz target
%
%  Plots:
%   (a) I-V curves (identical — device physics unchanged)
%   (b) Spike waveforms overlaid (time-normalised)
%   (c) Key metrics: RC, spike duration, max freq — bar charts
%   (d) Capacitance vs max frequency (design space)
%   (e) Inductance vs max frequency (parasitic elimination)
%   (f) THz design roadmap (C, L targets)
% =========================================================================

function RTD_RegimeComparison()

regime_names = {'Experiment', 'GHz_Sim', 'NearTHz', 'THz'};
labels       = {'Experiment\n(~MHz)', 'GHz Sim\n(Paper)', 'Near-THz\n(~100 GHz)', 'THz Target'};
colors       = [0.89 0.15 0.15;
                0.11 0.62 0.46;
                0.60 0.40 0.80;
                0.21 0.47 0.75];

Ps = cell(1,4);
for k = 1:4
    Ps{k} = get_regime_params(regime_names{k});
end

kBTq = Ps{1}.kBTq;
fV   = @(P_,V) P_.A .* log((1+exp((P_.b-P_.c+P_.n1.*V)./kBTq))./ ...
                              (1+exp((P_.b-P_.c-P_.n1.*V)./kBTq))).* ...
               (pi/2+atan((P_.c-P_.n1.*V)./P_.d)) + ...
               P_.H.*(exp(P_.n2.*V./kBTq)-1);

% ----- Derived metrics ---------------------------------------------------
tau_spikes = zeros(1,4);
max_freqs  = zeros(1,4);
RC_times   = zeros(1,4);
for k=1:4
    P_ = Ps{k};
    tau_spikes(k) = 2*pi*sqrt(P_.L*P_.C);
    max_freqs(k)  = 1/(3*tau_spikes(k));
    RC_times(k)   = P_.R * P_.C;
end

fprintf('   Regime comparison:\n');
for k=1:4
    fprintf('   %-12s  tau=%-10s  fmax=%-12s  RC=%-10s\n', ...
        regime_names{k}, format_time(tau_spikes(k)), ...
        format_freq(max_freqs(k)), format_time(RC_times(k)));
end

% ===== Plotting ===========================================================
fig = figure('Name','Module 6: Regime Comparison','NumberTitle','off', ...
             'Position',[160 50 1400 1000]);

% ---- (a) I-V curves (same for all regimes) ------------------------------
subplot(3,3,1);
V_range = linspace(0, 1.2, 2000);
I_range = fV(Ps{1}, V_range) * 1000;
plot(V_range, I_range, 'k-', 'LineWidth', 2.5);
hold on;
% Find and shade NDR
[~, ipk] = max(fV(Ps{1}, V_range));
sub_ = fV(Ps{1}, V_range(ipk:end));
[~, iv_rel] = min(sub_);
ival = ipk + iv_rel - 1;
fill([V_range(ipk:ival), fliplr(V_range(ipk:ival))], ...
     [I_range(ipk:ival), min(I_range(ipk:ival))*ones(1,ival-ipk+1)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor','none');
xlabel('Voltage (V)'); ylabel('Current (mA)');
title('(a) I–V Curve (same all regimes)');
text(0.7, max(I_range)*0.8, 'NDR region', 'Color','r','FontSize',9);
grid on;

% ---- (b) Spike waveforms — time-normalised per regime -------------------
subplot(3,3,2);
hold on;
legend_entries = {};
for k=1:4
    P_   = Ps{k};
    ts   = 2*pi*sqrt(P_.L*P_.C);
    dt   = ts/300;
    t    = 0:dt:ts*8;
    N    = length(t);
    V    = zeros(1,N); I_ = zeros(1,N);
    V(1) = 0.75; I_(1) = fV(P_, 0.75);
    pulse_dur = ts*0.3;

    for step=2:N
        S0k = 600*(t(step-1) >= ts*0.3 & t(step-1) <= ts*0.3+pulse_dur);
        dV  = (I_(step-1) - fV(P_,V(step-1)) - P_.kappa*S0k)/P_.C;
        dI  = (0.75 - V(step-1) - P_.R*I_(step-1))/P_.L;
        V(step)  = V(step-1) + dt*dV;
        I_(step) = I_(step-1) + dt*dI;
    end

    t_norm = t ./ ts;  % normalise to spike period
    plot(t_norm, V, 'Color', colors(k,:), 'LineWidth', 2.0);
    legend_entries{k} = strrep(regime_names{k},'_',' ');
end
xlabel('Time / \tau_{spike} (normalised)');
ylabel('RTD Voltage (V)');
title('(b) Spike Waveforms — Normalised Time');
legend(legend_entries, 'Location','northeast','FontSize',8);
grid on;

% ---- (c) Spike duration bar chart ---------------------------------------
subplot(3,3,3);
spike_ps = tau_spikes * 1e12;
b1 = bar(spike_ps, 'FaceColor','flat');
for k=1:4, b1.CData(k,:) = colors(k,:); end
set(gca,'XTickLabel',{'Exp','GHz','~THz','THz'},'YScale','log');
ylabel('Spike Duration (ps)');
title('(c) Spike Duration per Regime');
grid on;
yline(1, 'r--', '1 ps (THz)', 'LineWidth',1.5,'LabelVerticalAlignment','bottom');
for k=1:4
    text(k, spike_ps(k)*1.5, format_time(tau_spikes(k)), ...
         'HorizontalAlignment','center','FontSize',7,'Rotation',0);
end

% ---- (d) Max frequency bar chart ----------------------------------------
subplot(3,3,4);
freq_GHz = max_freqs / 1e9;
b2 = bar(freq_GHz, 'FaceColor','flat');
for k=1:4, b2.CData(k,:) = colors(k,:); end
set(gca,'XTickLabel',{'Exp','GHz','~THz','THz'},'YScale','log');
ylabel('Max Spike Frequency (GHz)');
title('(d) Max Frequency per Regime');
grid on;
yline(1e3, 'r--', '1 THz', 'LineWidth',1.5,'LabelVerticalAlignment','bottom');
for k=1:4
    text(k, freq_GHz(k)*2, format_freq(max_freqs(k)), ...
         'HorizontalAlignment','center','FontSize',7);
end

% ---- (e) Design space: C vs max frequency -------------------------------
subplot(3,3,5);
C_range = logspace(-19, -12, 200);    % 1 aF to 1 pF
L_fixed = 1e-12;                       % 1 pH (THz target inductance)
f_vs_C  = 1 ./ (3 * 2*pi*sqrt(L_fixed .* C_range));
semilogx(C_range*1e15, f_vs_C/1e9, 'k-', 'LineWidth', 2); hold on;

% Mark regimes on this curve
for k=1:4
    P_ = Ps{k};
    f_k = 1/(3*2*pi*sqrt(L_fixed*P_.C));
    plot(P_.C*1e15, f_k/1e9, 'o', 'Color', colors(k,:), ...
         'MarkerFaceColor', colors(k,:), 'MarkerSize', 10);
end
xlabel('Capacitance (fF)');
ylabel('Max Frequency (GHz)');
title('(e) Design Space: C vs f_{max} (L=1pH fixed)');
set(gca, 'XScale','log','YScale','log');
yline(1e3,'r--','1 THz','LineWidth',1.5,'LabelVerticalAlignment','bottom');
legend([{'f_{max}(C)'},strrep(regime_names,'_',' ')], 'FontSize',7,'Location','southwest');
grid on;

% ---- (f) Design space: L vs max frequency -------------------------------
subplot(3,3,6);
L_range = logspace(-13, -6, 200);    % 0.1 pH to 1 uH
C_fixed = 1e-18;                      % 1 aF (THz target capacitance)
f_vs_L  = 1 ./ (3 * 2*pi*sqrt(L_range .* C_fixed));
semilogx(L_range*1e12, f_vs_L/1e9, 'k-', 'LineWidth', 2); hold on;

for k=1:4
    P_ = Ps{k};
    f_k = 1/(3*2*pi*sqrt(P_.L*C_fixed));
    plot(P_.L*1e12, f_k/1e9, 'o', 'Color', colors(k,:), ...
         'MarkerFaceColor', colors(k,:), 'MarkerSize', 10);
end
xlabel('Inductance (pH)');
ylabel('Max Frequency (GHz)');
title('(f) Design Space: L vs f_{max} (C=1aF fixed)');
set(gca, 'XScale','log','YScale','log');
yline(1e3,'r--','1 THz','LineWidth',1.5,'LabelVerticalAlignment','bottom');
grid on;

% ---- (g) THz design roadmap — 2D contour map ---------------------------
subplot(3,3,[7 8 9]);
C_2d = logspace(-19, -12, 80);
L_2d = logspace(-13, -6, 80);
[C_grid, L_grid] = meshgrid(C_2d, L_2d);
F_grid = 1 ./ (3 * 2*pi*sqrt(L_grid .* C_grid));

contourf(log10(C_2d*1e15), log10(L_2d*1e12), log10(F_grid/1e9), 20, ...
         'LineColor','none');
colorbar('Label','log_{10}(f_{max} / GHz)');
colormap(parula);
hold on;

% Regime markers
for k=1:4
    P_ = Ps{k};
    plot(log10(P_.C*1e15), log10(P_.L*1e12), 's', ...
         'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
         'MarkerSize', 14, 'LineWidth', 2);
    text(log10(P_.C*1e15)+0.1, log10(P_.L*1e12), ...
         strrep(regime_names{k},'_',' '), 'FontSize',9,'Color','w','FontWeight','bold');
end

% THz target line
[C_thz_line, L_thz_line] = deal(logspace(-19,-12,100), logspace(-13,-6,100));
% f = 1THz line: L*C = (1/(3*2*pi*1e12))^2
LC_1THz = (1/(3*2*pi*1e12))^2;
C_thz_curve = logspace(-19,-12,200);
L_thz_curve = LC_1THz ./ C_thz_curve;
valid = L_thz_curve >= min(L_2d) & L_thz_curve <= max(L_2d);
plot(log10(C_thz_curve(valid)*1e15), log10(L_thz_curve(valid)*1e12), ...
     'r--', 'LineWidth', 2.5);
text(-2.5, -2, '1 THz boundary', 'Color','r','FontSize',10,'FontWeight','bold');

xlabel('log_{10}(C / fF)');
ylabel('log_{10}(L / pH)');
title('(g) THz Design Roadmap: C–L Space | Red dashed = 1 THz boundary');
grid on;

sgtitle('Module 6: Regime Comparison — Experiment to THz Target', ...
        'FontWeight','bold','FontSize',12);
drawnow;
end


function s = format_time(t)
    if t<1e-12, s=sprintf('%.2f fs',t*1e15);
    elseif t<1e-9, s=sprintf('%.2f ps',t*1e12);
    elseif t<1e-6, s=sprintf('%.2f ns',t*1e9);
    else, s=sprintf('%.2f us',t*1e6); end
end
function s = format_freq(f)
    if f>=1e12, s=sprintf('%.2f THz',f/1e12);
    elseif f>=1e9, s=sprintf('%.2f GHz',f/1e9);
    elseif f>=1e6, s=sprintf('%.2f MHz',f/1e6);
    else, s=sprintf('%.1f kHz',f/1e3); end
end
