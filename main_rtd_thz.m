clc; clear; close all;

fprintf('=================================================================\n');
fprintf('  RTD THz NEUROMORPHIC SIMULATION SUITE\n');
fprintf('  Owen-Newns et al. 2026 — Extended to THz Regime\n');
fprintf('=================================================================\n\n');

% ----- Global parameter sets ----------------------------------------------
REGIME = 'THz';   % Options: 'Experiment' | 'GHz_Sim' | 'NearTHz' | 'THz'
P = get_regime_params(REGIME);

fprintf('Active regime : %s\n', REGIME);
fprintf('  C = %.4e F\n', P.C);
fprintf('  L = %.4e H\n', P.L);
fprintf('  R = %.1f Ohm\n', P.R);
fprintf('  Vp = %.4f V  (peak voltage)\n',  P.Vp);
fprintf('  Ip = %.4e A  (peak current)\n',  P.Ip);
fprintf('  Vv = %.4f V  (valley voltage)\n', P.Vv);
fprintf('  Iv = %.4e A  (valley current)\n', P.Iv);
fprintf('  Spike duration (approx): %s\n', format_time(pi * sqrt(P.L * P.C)));
fprintf('  Max spike freq  (approx): %s\n', format_freq(1 / (6 * pi * sqrt(P.L * P.C))));
fprintf('\n');

% ----- Run all modules ----------------------------------------------------
% Note: Ensure these .m files exist in your current MATLAB path.

modules = {
    'RTD_IV_Model',      'I-V curve and NDR characterisation';
    'RTD_SingleNeuron',  'Single neuron spike dynamics';
};

nModules = size(modules, 1);

for k = 1:nModules
    modName = modules{k, 1};
    modDesc = modules{k, 2};
    fprintf('[%d/%d] %s...\n', k, nModules, modDesc);
    try
        feval(modName, P);
    catch ME
        fprintf('  WARNING: Module "%s" failed.\n', modName);
        fprintf('  Reason : %s\n', ME.message);
        fprintf('  Ensure %s.m exists and is on the MATLAB path.\n\n', modName);
    end
end

fprintf('\n=== All modules executed. ===\n');

% =========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function P = get_regime_params(regime)
% Returns circuit and NDR I-V parameters for the chosen regime.
% NDR parameters (Vp, Ip, Vv, Iv) scale with the LC time constant so that
% the peak/valley voltages remain physically meaningful at each frequency.

    switch regime
        case 'Experiment'
            P.C  = 10e-12;   % F
            P.L  = 50e-9;    % H
            P.R  = 10;       % Ohm
            P.Vp = 0.20;     % V  — peak voltage
            P.Ip = 5e-3;     % A  — peak current
            P.Vv = 0.45;     % V  — valley voltage
            P.Iv = 1e-3;     % A  — valley current

        case 'GHz_Sim'
            P.C  = 1e-12;
            P.L  = 10e-9;
            P.R  = 5;
            P.Vp = 0.18;
            P.Ip = 2e-3;
            P.Vv = 0.40;
            P.Iv = 0.4e-3;

        case 'NearTHz'
            P.C  = 100e-15;
            P.L  = 1e-9;
            P.R  = 2;
            P.Vp = 0.15;
            P.Ip = 1e-3;
            P.Vv = 0.35;
            P.Iv = 0.2e-3;

        case 'THz'
            P.C  = 10e-15;
            P.L  = 100e-12;
            P.R  = 1;
            P.Vp = 0.12;
            P.Ip = 0.5e-3;
            P.Vv = 0.30;
            P.Iv = 0.1e-3;

        otherwise
            error('Unknown regime "%s". Choose: Experiment | GHz_Sim | NearTHz | THz', regime);
    end
end

% -------------------------------------------------------------------------

function str = format_time(t)
    if t < 1e-12
        str = sprintf('%.2f fs', t * 1e15);
    elseif t < 1e-9
        str = sprintf('%.2f ps', t * 1e12);
    else
        str = sprintf('%.2f ns', t * 1e9);
    end
end

% -------------------------------------------------------------------------

function str = format_freq(f)
    if f >= 1e12
        str = sprintf('%.2f THz', f / 1e12);
    elseif f >= 1e9
        str = sprintf('%.2f GHz', f / 1e9);
    else
        str = sprintf('%.2f MHz', f / 1e6);
    end
end
