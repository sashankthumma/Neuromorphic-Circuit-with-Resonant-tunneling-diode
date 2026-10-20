% =========================================================================
%  RTD_pSNN_ELM.m
%  Module 4: Photonic Spiking Extreme Learning Machine
%
%  Implements Section 3.2.1 of Owen-Newns et al. 2026
%  20 uncoupled RTD neurons, Iris flower dataset classification
%  Heaviside step activation, least squares + node-significance training
%
%  Plots:
%   (a) Network architecture diagram (text-based)
%   (b) Mean spike activity maps per class
%   (c) Binary spike pattern matrix (all 150 datapoints)
%   (d) Accuracy vs training set size — least squares
%   (e) Accuracy vs training size & nodes — node significance
%   (f) Confusion matrix
% =========================================================================

function RTD_pSNN_ELM(P)

if nargin < 1
    P = get_regime_params('THz');
end

kBTq = P.kBTq;
fV   = @(V) P.A .* log((1+exp((P.b-P.c+P.n1.*V)./kBTq))./ ...
                         (1+exp((P.b-P.c-P.n1.*V)./kBTq))).* ...
            (pi/2+atan((P.c-P.n1.*V)./P.d)) + ...
            P.H.*(exp(P.n2.*V./kBTq)-1);

rng(42);  % reproducible results

% =========== Iris Dataset (embedded) =====================================
% 4 features: sepal length, sepal width, petal length, petal width
% 3 classes: Setosa(0), Versicolor(1), Virginica(2) — 50 each
[X_iris, y_iris] = load_iris();
N_data    = size(X_iris, 1);   % 150
N_feat    = size(X_iris, 2);   % 4
N_class   = 3;
N_neurons = 20;                % RTD array size

% Normalise features to [0,1]
X_norm = (X_iris - min(X_iris)) ./ (max(X_iris) - min(X_iris) + 1e-8);

% =========== Random input mixing (Nf x N weight matrix) ==================
W_in = randn(N_feat, N_neurons);   % random mixing matrix

% =========== RTD threshold (Heaviside activation) ========================
% Each RTD fires if weighted input exceeds its threshold T_n
% Threshold is set by bias-dependent NDR position
V0_array  = 0.75 * ones(1, N_neurons);  % all biased at valley
S_max_arr = 600;  % max optical amplitude

% Compute weighted inputs for all datapoints
X_weighted = X_norm * W_in;  % [150 x 20]

% Normalise to optical amplitude range [0, S_max_arr]
X_scaled = (X_weighted - min(X_weighted(:))) ./ ...
           (max(X_weighted(:)) - min(X_weighted(:))) * S_max_arr;

% =========== RTD spike firing (Heaviside step function) ==================
% For each RTD neuron n, threshold T_n determined by bias and NDR
% Spike if optical input > T_n
T_neurons = rand(1, N_neurons) * S_max_arr * 0.6 + S_max_arr * 0.2;

spike_matrix = double(X_scaled > repmat(T_neurons, N_data, 1));  % [150 x 20]

fprintf('   Mean spike rate: %.2f\n', mean(spike_matrix(:)));

% =========== Per-class mean spike activity ================================
mean_activity = zeros(N_class, N_neurons);
for c = 0:N_class-1
    idx_c = find(y_iris == c);
    mean_activity(c+1, :) = mean(spike_matrix(idx_c, :), 1);
end

% =========== Training: Least Squares =====================================
train_sizes  = 5:5:45;
acc_ls       = zeros(size(train_sizes));

for ti = 1:length(train_sizes)
    n_train = train_sizes(ti);
    [acc_ls(ti)] = train_test_ls(spike_matrix, y_iris, n_train, N_class);
end

fprintf('   Peak LS accuracy: %.1f%%\n', max(acc_ls)*100);

% =========== Training: Node Significance Algorithm =======================
node_range = 2:2:20;
ns_acc     = zeros(length(train_sizes), length(node_range));

for ti = 1:length(train_sizes)
    for ni = 1:length(node_range)
        ns_acc(ti,ni) = train_test_node_sig(spike_matrix, y_iris, ...
                         train_sizes(ti), N_class, node_range(ni));
    end
end

fprintf('   Peak Node-Sig accuracy: %.1f%%\n', max(ns_acc(:))*100);

% =========== Confusion matrix (best LS model) ============================
[~, conf_mat] = train_test_ls(spike_matrix, y_iris, 39, N_class);

% ===== Plotting ===========================================================
fig = figure('Name','Module 4: Photonic Spiking ELM','NumberTitle','off', ...
             'Position',[120 50 1350 950]);

class_names = {'Setosa','Versicolor','Virginica'};
cmap_spike  = [0.27 0.18 0.50; 1 0.90 0.35];  % purple (spike) / yellow (no spike)

% (a) Mean spike activity maps
for c = 1:N_class
    ax = subplot(3, 4, c);
    activity_grid = reshape(mean_activity(c,:), 4, 5);
    imagesc(activity_grid);
    colormap(ax, hot);
    colorbar;
    caxis([0 1]);
    title(sprintf('(a) Mean activity\n%s', class_names{c}));
    xlabel('RTD column'); ylabel('RTD row');
    set(gca,'XTick',1:5,'YTick',1:4);
end

% (b) Spike pattern matrix — all 150 datapoints
subplot(3, 4, 4);
imagesc(spike_matrix');
colormap(gca, cmap_spike);
xlabel('Datapoint index');
ylabel('RTD neuron');
title('(b) Binary spike pattern matrix');
hold on;
xline(50, 'w-', 'LineWidth', 2);
xline(100, 'w-', 'LineWidth', 2);
text(25,  21, 'Setosa',     'Color','w','FontSize',8,'HorizontalAlignment','center');
text(75,  21, 'Versicolor', 'Color','w','FontSize',8,'HorizontalAlignment','center');
text(125, 21, 'Virginica',  'Color','w','FontSize',8,'HorizontalAlignment','center');

% (c) LS accuracy
subplot(3, 4, [5 6]);
plot(train_sizes, acc_ls*100, 'b-o', 'LineWidth', 2, 'MarkerFaceColor','b', 'MarkerSize',6);
hold on;
yline(96.5, 'r--', '96.5% (paper)', 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
xlabel('Training samples per class');
ylabel('Test Accuracy (%)');
title('(c) Least Squares Training Accuracy');
ylim([60 100]); grid on;
text(train_sizes(end), max(acc_ls)*100+1, sprintf('Peak: %.1f%%', max(acc_ls)*100), ...
     'FontSize', 9, 'Color','b');

% (d) Node significance heatmap
subplot(3, 4, [7 8]);
imagesc(train_sizes, node_range, ns_acc'*100);
colorbar;
xlabel('Training samples per class');
ylabel('Number of significant nodes');
title('(d) Node Significance Algorithm Accuracy (%)');
colormap(gca, parula);
caxis([60 100]);

% (e) Comparison LS vs NS at best conditions
subplot(3, 4, [9 10]);
best_ns = max(ns_acc, [], 2);
plot(train_sizes, acc_ls*100, 'b-o', 'LineWidth', 2, 'MarkerFaceColor','b', 'MarkerSize',5);
hold on;
plot(train_sizes, best_ns*100, 'r-s', 'LineWidth', 2, 'MarkerFaceColor','r', 'MarkerSize',5);
yline(96.5, 'b:', 'LS Paper', 'LineWidth',1);
yline(93.0, 'r:', 'NS Paper', 'LineWidth',1);
xlabel('Training samples per class');
ylabel('Accuracy (%)');
title('(e) LS vs Node-Significance Comparison');
legend({'Least Squares','Node Significance','Paper LS','Paper NS'}, ...
       'Location','southeast','FontSize',8);
ylim([60 100]); grid on;

% (f) Confusion matrix
subplot(3, 4, [11 12]);
imagesc(conf_mat);
colormap(gca, flipud(hot));
colorbar;
threshold = max(conf_mat(:)) / 2 * 0.7;
for ri = 1:N_class
    for ci = 1:N_class
        if conf_mat(ri,ci) > threshold
            txtColor = 'w';
        else
            txtColor = 'k';
        end
        text(ci, ri, num2str(conf_mat(ri,ci)), ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', ...
             'Color', txtColor);
    end
end

% Ensure class_names is a cell array or string array with 3 entries, e.g.:
% class_names = {'Setosa','Versicolor','Virginica'};

set(gca, ...
    'XTick', 1:3, ...
    'XTickLabel', class_names, ...
    'YTick', 1:3, ...
    'YTickLabel', class_names);

xlabel('Predicted');
ylabel('True');
title('(f) Confusion Matrix (39 train/class)');

% Ensure P.regime_label is defined (char or string)
sgtitle(sprintf('Module 4: Photonic Spiking ELM — Iris Classification | %s', P.regime_label), ...
       'FontWeight', 'bold', 'FontSize', 12);

drawnow;
% end  % only include if this matches an open block (function/if/for/while)

% =========================================================================
%  Iris dataset (built-in — no file needed)
% =========================================================================
function [X, y] = load_iris()
    % Iris dataset: 150x4 features, 150x1 labels
    % Using Fisher's original 1936 data
    data = [
    5.1,3.5,1.4,0.2; 4.9,3.0,1.4,0.2; 4.7,3.2,1.3,0.2; 4.6,3.1,1.5,0.2;
    5.0,3.6,1.4,0.2; 5.4,3.9,1.7,0.4; 4.6,3.4,1.4,0.3; 5.0,3.4,1.5,0.2;
    4.4,2.9,1.4,0.2; 4.9,3.1,1.5,0.1; 5.4,3.7,1.5,0.2; 4.8,3.4,1.6,0.2;
    4.8,3.0,1.4,0.1; 4.3,3.0,1.1,0.1; 5.8,4.0,1.2,0.2; 5.7,4.4,1.5,0.4;
    5.4,3.9,1.3,0.4; 5.1,3.5,1.4,0.3; 5.7,3.8,1.7,0.3; 5.1,3.8,1.5,0.3;
    5.4,3.4,1.7,0.2; 5.1,3.7,1.5,0.4; 4.6,3.6,1.0,0.2; 5.1,3.3,1.7,0.5;
    4.8,3.4,1.9,0.2; 5.0,3.0,1.6,0.2; 5.0,3.4,1.6,0.4; 5.2,3.5,1.5,0.2;
    5.2,3.4,1.4,0.2; 4.7,3.2,1.6,0.2; 4.8,3.1,1.6,0.2; 5.4,3.4,1.5,0.4;
    5.2,4.1,1.5,0.1; 5.5,4.2,1.4,0.2; 4.9,3.1,1.5,0.1; 5.0,3.2,1.2,0.2;
    5.5,3.5,1.3,0.2; 4.9,3.1,1.5,0.1; 4.4,3.0,1.3,0.2; 5.1,3.4,1.5,0.2;
    5.0,3.5,1.3,0.3; 4.5,2.3,1.3,0.3; 4.4,3.2,1.3,0.2; 5.0,3.5,1.6,0.6;
    5.1,3.8,1.9,0.4; 4.8,3.0,1.4,0.3; 5.1,3.8,1.6,0.2; 4.6,3.2,1.4,0.2;
    5.3,3.7,1.5,0.2; 5.0,3.3,1.4,0.2;
    7.0,3.2,4.7,1.4; 6.4,3.2,4.5,1.5; 6.9,3.1,4.9,1.5; 5.5,2.3,4.0,1.3;
    6.5,2.8,4.6,1.5; 5.7,2.8,4.5,1.3; 6.3,3.3,4.7,1.6; 4.9,2.4,3.3,1.0;
    6.6,2.9,4.6,1.3; 5.2,2.7,3.9,1.4; 5.0,2.0,3.5,1.0; 5.9,3.0,4.2,1.5;
    6.0,2.2,4.0,1.0; 6.1,2.9,4.7,1.4; 5.6,2.9,3.6,1.3; 6.7,3.1,4.4,1.4;
    5.6,3.0,4.5,1.5; 5.8,2.7,4.1,1.0; 6.2,2.2,4.5,1.5; 5.6,2.5,3.9,1.1;
    5.9,3.2,4.8,1.8; 6.1,2.8,4.0,1.3; 6.3,2.5,4.9,1.5; 6.1,2.8,4.7,1.2;
    6.4,2.9,4.3,1.3; 6.6,3.0,4.4,1.4; 6.8,2.8,4.8,1.4; 6.7,3.0,5.0,1.7;
    6.0,2.9,4.5,1.5; 5.7,2.6,3.5,1.0; 5.5,2.4,3.8,1.1; 5.5,2.4,3.7,1.0;
    5.8,2.7,3.9,1.2; 6.0,2.7,5.1,1.6; 5.4,3.0,4.5,1.5; 6.0,3.4,4.5,1.6;
    6.7,3.1,4.7,1.5; 6.3,2.3,4.4,1.3; 5.6,3.0,4.1,1.3; 5.5,2.5,4.0,1.3;
    5.5,2.6,4.4,1.2; 6.1,3.0,4.6,1.4; 5.8,2.6,4.0,1.2; 5.0,2.3,3.3,1.0;
    5.6,2.7,4.2,1.3; 5.7,3.0,4.2,1.2; 5.7,2.9,4.2,1.3; 6.2,2.9,4.3,1.3;
    5.1,2.5,3.0,1.1; 5.7,2.8,4.1,1.3;
    6.3,3.3,6.0,2.5; 5.8,2.7,5.1,1.9; 7.1,3.0,5.9,2.1; 6.3,2.9,5.6,1.8;
    6.5,3.0,5.8,2.2; 7.6,3.0,6.6,2.1; 4.9,2.5,4.5,1.7; 7.3,2.9,6.3,1.8;
    6.7,2.5,5.8,1.8; 7.2,3.6,6.1,2.5; 6.5,3.2,5.1,2.0; 6.4,2.7,5.3,1.9;
    6.8,3.0,5.5,2.1; 5.7,2.5,5.0,2.0; 5.8,2.8,5.1,2.4; 6.4,3.2,5.3,2.3;
    6.5,3.0,5.5,1.8; 7.7,3.8,6.7,2.2; 7.7,2.6,6.9,2.3; 6.0,2.2,5.0,1.5;
    6.9,3.2,5.7,2.3; 5.6,2.8,4.9,2.0; 7.7,2.8,6.7,2.0; 6.3,2.7,4.9,1.8;
    6.7,3.3,5.7,2.1; 7.2,3.2,6.0,1.8; 6.2,2.8,4.8,1.8; 6.1,3.0,4.9,1.8;
    6.4,2.8,5.6,2.1; 7.2,3.0,5.8,1.6; 7.4,2.8,6.1,1.9; 7.9,3.8,6.4,2.0;
    6.4,2.8,5.6,2.2; 6.3,2.8,5.1,1.5; 6.1,2.6,5.6,1.4; 7.7,3.0,6.1,2.3;
    6.3,3.4,5.6,2.4; 6.4,3.1,5.5,1.8; 6.0,3.0,4.8,1.8; 6.9,3.1,5.4,2.1;
    6.7,3.1,5.6,2.4; 6.9,3.1,5.1,2.3; 5.8,2.7,5.1,1.9; 6.8,3.2,5.9,2.3;
    6.7,3.3,5.7,2.5; 6.7,3.0,5.2,2.3; 6.3,2.5,5.0,1.9; 6.5,3.0,5.2,2.0;
    6.2,3.4,5.4,2.3; 5.9,3.0,5.1,1.8];

    y = [zeros(50,1); ones(50,1); 2*ones(50,1)];
    X = data;
end


% =========================================================================
%  Least-squares training and testing
% =========================================================================
function [acc, conf] = train_test_ls(spike_mat, y, n_train, n_class)
    N = size(spike_mat, 1);
    idx_train = []; idx_test = [];
    for c = 0:n_class-1
        class_idx = find(y == c);
        perm = randperm(length(class_idx));
        idx_train = [idx_train; class_idx(perm(1:n_train))];
        idx_test  = [idx_test;  class_idx(perm(n_train+1:end))];
    end

    X_tr = spike_mat(idx_train, :);
    X_te = spike_mat(idx_test, :);
    y_tr = y(idx_train);
    y_te = y(idx_test);

    % One-hot encode
    Y_oh = (y_tr == (0:n_class-1));
    W_out = pinv(X_tr) * Y_oh;

    % Predict
    Y_pred = X_te * W_out;
    [~, pred_class] = max(Y_pred, [], 2);
    pred_class = pred_class - 1;

    acc = mean(pred_class == y_te);

    conf = zeros(n_class);
    for i = 1:length(y_te)
        conf(y_te(i)+1, pred_class(i)+1) = conf(y_te(i)+1, pred_class(i)+1) + 1;
    end
end


% =========================================================================
%  Node significance training
% =========================================================================
function acc = train_test_node_sig(spike_mat, y, n_train, n_class, n_nodes)
    N_neurons = size(spike_mat, 2);
    idx_train = []; idx_test = [];
    for c = 0:n_class-1
        class_idx = find(y == c);
        perm = randperm(length(class_idx));
        idx_train = [idx_train; class_idx(perm(1:min(n_train,length(class_idx))))];
        idx_test  = [idx_test;  class_idx(perm(min(n_train,length(class_idx))+1:end))];
    end
    if isempty(idx_test), acc = 0; return; end

    X_tr = spike_mat(idx_train, :);
    X_te = spike_mat(idx_test, :);
    y_tr = y(idx_train);
    y_te = y(idx_test);

    % Node significance: per-class spike rate difference
    sig_scores = zeros(N_neurons, 1);
    for n = 1:N_neurons
        rates = zeros(n_class,1);
        for c = 0:n_class-1
            rates(c+1) = mean(X_tr(y_tr==c, n));
        end
        sig_scores(n) = std(rates);  % high std = more discriminative
    end

    [~, sorted_nodes] = sort(sig_scores, 'descend');
    sel_nodes = sorted_nodes(1:min(n_nodes, N_neurons));

    % Binary weights: 1 if node indicates this class most
    W_bin = zeros(N_neurons, n_class);
    for n = sel_nodes'
        rates = zeros(n_class,1);
        for c = 0:n_class-1
            rates(c+1) = mean(X_tr(y_tr==c, n));
        end
        [~, best_class] = max(rates);
        W_bin(n, best_class) = 1;
    end

    Y_pred = X_te * W_bin;
    [~, pred_class] = max(Y_pred, [], 2);
    pred_class = pred_class - 1;
    acc = mean(pred_class == y_te);
end