% =========================================================================
%  get_regime_params.m
%  Returns RTD circuit parameters for a given operating regime.
%
%  Regimes:
%   'Experiment' — Real experimental circuit (paper Section 5)
%                  Large parasitics, ~150 ns spike duration
%   'GHz_Sim'    — Paper numerical simulation (Table 1)
%                  ~300 ps spike duration, ~GHz max rate
%   'NearTHz'    — Optimised integrated circuit (projected)
%                  ~10 ps spike duration, ~100 GHz
%   'THz'        — Ideal monolithic device (target)
%                  ~0.5 ps spike duration, ~THz
% =========================================================================

function P = get_regime_params(regime)

% ----- Shared Schulman f(V) parameters (same for all regimes) ------------
P.A   = -5.5e-5;   % [A]
P.b   =  0.033;    % [V]
P.c   =  0.113;    % [V]
P.d   = -2.8e-6;   % [V]
P.n1  =  0.185;
P.n2  =  0.045;
P.H   =  18e-5;    % [A]
P.kappa = 0.1e-6;  % [A / arb optical unit]
P.kBTq  = 0.02585; % kT/q at 300K [V]

% ----- Laser parameters (two-level rate equation model) ------------------
P.gamma_m  = 1e10;   % spontaneous emission rate [s^-1]
P.gamma_l  = 1e9;    % radiative leakage rate [s^-1]
P.gamma_nr = 1e9;    % non-radiative leakage rate [s^-1]
P.N0       = 1e8;    % threshold inversion number
P.tau_ph   = 3e-12;  % photon lifetime [s]
P.eta      = 0.5;    % coupling efficiency RTD->laser
P.J_laser  = 2e10;   % laser bias current (normalised)

% ----- Circuit parameters — regime dependent -----------------------------
switch regime
    case 'Experiment'
        % Real lab circuit — bond wires, bias-tee, connectors
        P.R = 10;           % [Ohm]
        P.L = 126e-3;       % [H]  — 126 mH equivalent parasitics
        P.C = 2e-12;        % [F]  — 2 pF
        P.regime_label = 'Experimental (~MHz)';
        P.color = [0.89 0.15 0.15];

    case 'GHz_Sim'
        % Paper Table 1 — numerical simulation
        P.R = 10;           % [Ohm]
        P.L = 126e-9;       % [H]  — 126 nH
        P.C = 0.002e-12;    % [F]  — 0.002 pF = 2 fF
        P.regime_label = 'GHz Simulation (Paper)';
        P.color = [0.11 0.62 0.46];

    case 'NearTHz'
        % Projected optimised integrated circuit
        P.R = 10;           % [Ohm]
        P.L = 10e-12;       % [H]  — 10 pH
        P.C = 0.05e-15;     % [F]  — 0.05 fF
        P.regime_label = 'Near-THz (~100 GHz)';
        P.color = [0.33 0.29 0.72];

    case 'THz'
        % Ideal monolithic InGaAs/AlAs on InP — THz target
        % Mesa ~80 nm, all-on-chip, no external parasitics
        P.R = 10;           % [Ohm]
        P.L = 1e-12;        % [H]  — 1 pH  (intrinsic device only)
        P.C = 1e-18;        % [F]  — 0.001 fF = 1 aF
        P.regime_label = 'THz Target (~THz)';
        P.color = [0.21 0.47 0.75];

    otherwise
        error('Unknown regime: %s', regime);
end

P.regime = regime;
end


% =========================================================================
%  fV_RTD.m  —  Schulman RTD I-V conductivity model
%  Equation (3) from Owen-Newns et al. 2026
% =========================================================================
function I = fV_RTD(V, P)
    kBTq = P.kBTq;
    num  = 1 + exp((P.b - P.c + P.n1 .* V) ./ kBTq);
    den  = 1 + exp((P.b - P.c - P.n1 .* V) ./ kBTq);
    term1 = P.A .* log(num ./ den);
    term2 = (pi/2 + atan((P.c - P.n1 .* V) ./ P.d));
    term3 = P.H .* (exp(P.n2 .* V ./ kBTq) - 1);
    I = term1 .* term2 + term3;
end


% =========================================================================
%  rtd_ode.m  —  ODE right-hand side for RTD circuit
%  State: x = [V; I]
%  Equations (1) and (2) from paper
% =========================================================================
function dxdt = rtd_ode(~, x, V0, S0, P)
    V = x(1);
    I = x(2);
    dV = (I - fV_RTD(V, P) - P.kappa * S0) / P.C;
    dI = (V0 - V - P.R * I) / P.L;
    dxdt = [dV; dI];
end


% =========================================================================
%  rtd_laser_ode.m  —  ODE for RTD + laser coupled system
%  State: x = [V; I_rtd; N; S]
%  Equations (1-5) from paper
% =========================================================================
function dxdt = rtd_laser_ode(~, x, V0, S0_ext, P)
    V    = x(1);
    I    = x(2);
    N    = x(3);
    S    = x(4);

    % RTD equations
    dV = (I - fV_RTD(V, P) - P.kappa * S0_ext) / P.C;
    dI = (V0 - V - P.R * I) / P.L;

    % Laser rate equations — driven by RTD current I
    J_eff = P.J_laser + P.eta * I / 1.6e-19;
    dN = J_eff - (P.gamma_m + P.gamma_l + P.gamma_nr) * N ...
         - P.gamma_m * (N - P.N0) * S;
    dS = P.gamma_m * (N - P.N0) * S - S / P.tau_ph + P.gamma_m * N;

    dxdt = [dV; dI; dN; dS];
end


% =========================================================================
%  solve_rtd.m  —  Wrapper: solve RTD ODE over time span
% =========================================================================
function [t, V, I] = solve_rtd(tspan, V0_bias, S0_func, P, x0)
    if nargin < 5
        x0 = [V0_bias; fV_RTD(V0_bias, P)];
    end
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-12, 'MaxStep', P.L/P.R * 0.1);
    [t, x] = ode15s(@(t,x) rtd_ode(t, x, V0_bias, S0_func(t), P), ...
                    tspan, x0, opts);
    V = x(:,1);
    I = x(:,2);
end


% =========================================================================
%  get_ndr_info.m  —  Find peak/valley of I-V, compute PVCR
% =========================================================================
function ndr = get_ndr_info(P)
    V_scan = linspace(0.3, 1.1, 5000);
    I_scan = fV_RTD(V_scan, P);

    [~, ipk] = max(I_scan);
    peak_region = I_scan(ipk:end);
    [~, ival_rel] = min(peak_region);
    ival = ipk + ival_rel - 1;

    ndr.V_peak   = V_scan(ipk);
    ndr.I_peak   = I_scan(ipk);
    ndr.V_valley = V_scan(ival);
    ndr.I_valley = I_scan(ival);
    ndr.PVCR     = I_scan(ipk) / max(I_scan(ival), 1e-12);
    ndr.width_V  = ndr.V_valley - ndr.V_peak;
end


% =========================================================================
%  Utility: format_time, format_freq
% =========================================================================
function s = format_time(t)
    if t < 1e-12
        s = sprintf('%.2f fs', t * 1e15);
    elseif t < 1e-9
        s = sprintf('%.2f ps', t * 1e12);
    elseif t < 1e-6
        s = sprintf('%.2f ns', t * 1e9);
    elseif t < 1e-3
        s = sprintf('%.2f us', t * 1e6);
    else
        s = sprintf('%.2f ms', t * 1e3);
    end
end

function s = format_freq(f)
    if f >= 1e12
        s = sprintf('%.2f THz', f / 1e12);
    elseif f >= 1e9
        s = sprintf('%.2f GHz', f / 1e9);
    elseif f >= 1e6
        s = sprintf('%.2f MHz', f / 1e6);
    else
        s = sprintf('%.2f kHz', f / 1e3);
    end
end
