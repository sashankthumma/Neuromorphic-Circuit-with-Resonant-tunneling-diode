%% =========================================================
%  MODULE 5: Photonic SNN — Iris Flower Classification
%  Replicates Owen-Newns et al. 2023 (Paper 5)
%  VCSEL-based spiking reservoir computer
%  Time-multiplexed virtual nodes, θ = 250 ps
%  Task: Fisher's Iris dataset (150 flowers, 3 species)
%% =========================================================

function photonic_snn_iris()

rng(42);  % reproducibility

%% =========================================================
%  SECTION A: VCSEL Spiking Model
%% =========================================================

% VCSEL spiking neuron: leaky integrate-and-fire with
% injection locking and refractory period (Paper 5, Sec III)
theta_node  = 250e-12;   % node duration [s] = 250 ps
tau_refrac  = 1e-9;      % refractory period [s] = ~1 ns
tau_leak    = 400e-12;   % membrane leak [s]
spike_thresh = 1.0;       % normalized threshold

%% =========================================================
%  SECTION B: Load Iris Dataset
%% =========================================================
fprintf('   Loading Fisher Iris dataset...\n');

% Iris data: [SepalLength, SepalWidth, PetalLength, PetalWidth, Class]
% Class: 1=Setosa, 2=Versicolor, 3=Virginica
iris_data = load_iris_data();

X = iris_data(:, 1:4);   % 150x4
y = iris_data(:, 5);     % 150x1 labels (1,2,3)

% Normalize features to [0,1]
X_norm = (X - min(X)) ./ (max(X) - min(X));

n_samples  = size(X_norm, 1);   % 150
n_features = size(X_norm, 2);   % 4
n_classes  = 3;

%% =========================================================
%  SECTION C: SNN Parameters
%% =========================================================
Nv_list = [64, 128, 256, 512];   % virtual node counts to test
train_size = 30;                  % 10 per class (Paper 5)

%% =========================================================
%  SECTION D: Random Mask (input encoding)
%  Paper 5, Fig 6: random matrix Nv × 4 multiplied by data
%% =========================================================

results = struct();

for ni = 1:length(Nv_list)
    Nv = Nv_list(ni);
    fprintf('   Running SNN with Nv = %d virtual nodes...\n', Nv);

    % Random mask (fixed for all data points)
    mask = randn(Nv, n_features);

    %% --- Encode all data points through VCSEL SNN ---
    % Each flower → Nv-length binary spike pattern
    spike_patterns = zeros(n_samples, Nv);

    for s = 1:n_samples
        % Masked input vector
        u = mask * X_norm(s,:)';   % Nv × 1

        % Normalize input to reasonable drive range
        u_norm = (u - mean(u)) / (std(u) + 1e-8);

        % Simulate VCSEL spiking reservoir (leaky-IF)
        % Each element of u_norm drives one virtual node
        v_mem = zeros(Nv, 1);      % membrane voltage
        refrac = zeros(Nv, 1);     % refractory counter

        for node = 1:Nv
            drive = 0.5 + 0.4 * u_norm(node);  % bias + input

            % Couple from previous node (refractory creates memory)
            if node > 1
                coupling = 0.15 * spike_patterns(s, node-1);
            else
                coupling = 0;
            end

            % Leaky integrate
            v_mem(node) = (1 - theta_node/tau_leak) * 0 + ...
                          theta_node/tau_leak * (drive + coupling);

            % Threshold and fire
            if v_mem(node) >= spike_thresh && refrac(node) <= 0
                spike_patterns(s, node) = 1;
                refrac(node) = round(tau_refrac / theta_node);
            end

            % Propagate refractory
            if node < Nv
                refrac(node+1) = max(0, refrac(node) - 1);
            end
        end
    end

    %% --- Train/test split ---
    % Stratified: 10 per class for training (Paper 5)
    train_idx = [];
    test_idx  = [];
    for c = 1:n_classes
        cls_idx = find(y == c);
        cls_idx = cls_idx(randperm(length(cls_idx)));
        n_train_c = train_size / n_classes;
        train_idx = [train_idx; cls_idx(1:n_train_c)];
        test_idx  = [test_idx;  cls_idx(n_train_c+1:end)];
    end

    S_train = spike_patterns(train_idx, :);   % Nv-dim spike patterns
    S_test  = spike_patterns(test_idx,  :);
    y_train = y(train_idx);
    y_test  = y(test_idx);

    %% --- One-hot encode labels ---
    L_train = zeros(length(y_train), n_classes);
    for c = 1:n_classes
        L_train(:,c) = (y_train == c);
    end

    %% --- Train output layer (linear least squares, Paper 5 eq.1) ---
    % S * W = L  →  W = pinv(S) * L
    W = pinv(S_train) * L_train;   % Nv × 3

    %% --- Test ---
    out_test = S_test * W;          % n_test × 3
    [~, y_pred] = max(out_test, [], 2);

    accuracy = mean(y_pred == y_test) * 100;
    results(ni).Nv = Nv;
    results(ni).accuracy = accuracy;
    results(ni).spike_patterns = spike_patterns;
    results(ni).y = y;
    results(ni).W = W;

    fprintf('      Accuracy: %.1f%%\n', accuracy);
end

%% =========================================================
%  SECTION E: Plots
%% =========================================================

figure('Name','Photonic SNN — Iris Classification','Position',[50 50 1100 850]);

%% --- Plot 1: Spike pattern map (Nv=512 equivalent, use 256) ---
[~, best_idx] = max([results.accuracy]);
best = results(best_idx);
Nv_best = best.Nv;

subplot(3,3,[1 2]);
imagesc(best.spike_patterns);
colormap(gca, [1 1 1; 0.2 0.7 0.2]);  % white=no spike, green=spike
ylabel('Flower #','FontSize',11);
xlabel('Virtual node','FontSize',11);
title(sprintf('Spike patterns — Nv=%d nodes (best accuracy)', Nv_best),'FontSize',12);
hold on;
yline(50.5,'r-','LineWidth',2); yline(100.5,'b-','LineWidth',2);
text(Nv_best*0.02, 25, 'Setosa','Color','r','FontWeight','bold','FontSize',10);
text(Nv_best*0.02, 75, 'Versicolor','Color','b','FontWeight','bold','FontSize',10);
text(Nv_best*0.02, 125,'Virginica','Color','k','FontWeight','bold','FontSize',10);
colorbar('Ticks',[0.25 0.75],'TickLabels',{'No spike','Spike'},'FontSize',9);

%% --- Plot 2: Accuracy vs Nv ---
subplot(3,3,3);
Nv_vals  = [results.Nv];
acc_vals = [results.accuracy];
bar(Nv_vals, acc_vals, 0.5, 'FaceColor',[0.3 0.5 0.9],'EdgeColor','none');
hold on;
yline(97, '--r', '>97% (Paper 5)', 'LineWidth',1.5,'FontSize',9,...
      'LabelHorizontalAlignment','left');
xlabel('Virtual nodes N_v','FontSize',11);
ylabel('Classification accuracy (%)','FontSize',11);
title('Accuracy vs node count','FontSize',12);
ylim([50 102]); grid on;
xticks(Nv_vals);
for i=1:length(Nv_vals)
    text(Nv_vals(i), acc_vals(i)+1.5, sprintf('%.0f%%',acc_vals(i)), ...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end

%% --- Plot 3: Confusion matrix (best result) ---
subplot(3,3,[4 5]);
[~, best_idx2] = max([results.accuracy]);
best2 = results(best_idx2);
Nv2   = best2.Nv;

% Re-run to get confusion matrix
mask2 = randn(Nv2, n_features);
rng(42);
sp2 = zeros(n_samples, Nv2);
for s = 1:n_samples
    u = mask2 * X_norm(s,:)';
    u_norm = (u - mean(u)) / (std(u) + 1e-8);
    for node = 1:Nv2
        drive = 0.5 + 0.4*u_norm(node);
        coupling = 0;
        if node > 1, coupling = 0.15*sp2(s,node-1); end
        vm = theta_node/tau_leak * (drive+coupling);
        if vm >= spike_thresh, sp2(s,node) = 1; end
    end
end

train_idx2 = []; test_idx2 = [];
for c=1:n_classes
    ci = find(y==c); ci=ci(randperm(length(ci)));
    nc = train_size/n_classes;
    train_idx2=[train_idx2;ci(1:nc)]; test_idx2=[test_idx2;ci(nc+1:end)];
end
L2=zeros(length(train_idx2),n_classes);
for c=1:n_classes, L2(:,c)=(y(train_idx2)==c); end
W2=pinv(sp2(train_idx2,:))*L2;
out2=sp2(test_idx2,:)*W2;
[~,ypred2]=max(out2,[],2);
ytest2=y(test_idx2);

CM = confusionmat(ytest2, ypred2);
imagesc(CM); colormap(gca, parula);
colorbar;
xticks(1:3); yticks(1:3);
xticklabels({'Setosa','Versicolor','Virginica'});
yticklabels({'Setosa','Versicolor','Virginica'});
xlabel('Predicted','FontSize',11); ylabel('True','FontSize',11);
title(sprintf('Confusion matrix (Nv=%d, acc=%.1f%%)', Nv2, best2.accuracy),'FontSize',12);
for r=1:3, for c=1:3
    text(c,r,num2str(CM(r,c)),'HorizontalAlignment','center','FontSize',14,...
         'FontWeight','bold','Color','w');
end; end

%% --- Plot 4: Node firing rate per species ---
subplot(3,3,6);
sp_rates = zeros(3, Nv_best);
for c = 1:3
    sp_rates(c,:) = mean(best.spike_patterns(best.y==c,:), 1);
end
plot(sp_rates(1,:), 'r-', 'LineWidth',1.2,'DisplayName','Setosa'); hold on;
plot(sp_rates(2,:), 'b-', 'LineWidth',1.2,'DisplayName','Versicolor');
plot(sp_rates(3,:), 'k-', 'LineWidth',1.2,'DisplayName','Virginica');
xlabel('Virtual node index','FontSize',11);
ylabel('Mean firing rate','FontSize',11);
title('Per-species node firing rates','FontSize',12);
legend('FontSize',9,'Location','best'); grid on;

%% --- Plot 5: VCSEL SNN timing schematic ---
subplot(3,3,[7 8 9]);
% Show timing of virtual nodes for one example flower
flower_idx = 5;  % pick a Setosa
sp_flower = best.spike_patterns(flower_idx, :);
t_nodes = (0:Nv_best-1) * 250;   % time in ps

stem(t_nodes, sp_flower, 'Marker','none','Color',[0.2 0.6 0.9],'LineWidth',0.8);
hold on;
% Mark spike events
fired = t_nodes(sp_flower == 1);
plot(fired, ones(1,length(fired)), 'rv','MarkerSize',5,'MarkerFaceColor','r');
xlabel('Time (ps)', 'FontSize',11);
ylabel('Spike (1/0)','FontSize',11);
title(sprintf('Virtual node spike train — Flower #%d (Setosa) | θ=250 ps/node | Nv=%d', ...
              flower_idx, Nv_best),'FontSize',12);
xlim([0 Nv_best*250]); ylim([-0.1 1.5]);
text(0.7, 0.88, sprintf('Spike rate: %.1f%%', mean(sp_flower)*100), ...
     'Units','normalized','FontSize',10,'Color',[0.8 0.1 0.1]);
grid on;

sgtitle('Photonic SNN Reservoir Computer — Iris Flower Classification (Paper 5)', ...
        'FontSize',14,'FontWeight','bold');

fprintf('\n   === FINAL RESULTS ===\n');
for i=1:length(results)
    fprintf('   Nv=%4d → Accuracy: %.1f%%\n', results(i).Nv, results(i).accuracy);
end

end

%% =========================================================
%  HELPER: Fisher Iris Data (hardcoded — no toolbox needed)
%% =========================================================
function data = load_iris_data()
% [SepalLen, SepalWid, PetalLen, PetalWid, Class]
% Class: 1=Setosa, 2=Versicolor, 3=Virginica
raw = [
5.1 3.5 1.4 0.2 1; 4.9 3.0 1.4 0.2 1; 4.7 3.2 1.3 0.2 1;
4.6 3.1 1.5 0.2 1; 5.0 3.6 1.4 0.2 1; 5.4 3.9 1.7 0.4 1;
4.6 3.4 1.4 0.3 1; 5.0 3.4 1.5 0.2 1; 4.4 2.9 1.4 0.2 1;
4.9 3.1 1.5 0.1 1; 5.4 3.7 1.5 0.2 1; 4.8 3.4 1.6 0.2 1;
4.8 3.0 1.4 0.1 1; 4.3 3.0 1.1 0.1 1; 5.8 4.0 1.2 0.2 1;
5.7 4.4 1.5 0.4 1; 5.4 3.9 1.3 0.4 1; 5.1 3.5 1.4 0.3 1;
5.7 3.8 1.7 0.3 1; 5.1 3.8 1.5 0.3 1; 5.4 3.4 1.7 0.2 1;
5.1 3.7 1.5 0.4 1; 4.6 3.6 1.0 0.2 1; 5.1 3.3 1.7 0.5 1;
4.8 3.4 1.9 0.2 1; 5.0 3.0 1.6 0.2 1; 5.0 3.4 1.6 0.4 1;
5.2 3.5 1.5 0.2 1; 5.2 3.4 1.4 0.2 1; 4.7 3.2 1.6 0.2 1;
4.8 3.1 1.6 0.2 1; 5.4 3.4 1.5 0.4 1; 5.2 4.1 1.5 0.1 1;
5.5 4.2 1.4 0.2 1; 4.9 3.1 1.5 0.2 1; 5.0 3.2 1.2 0.2 1;
5.5 3.5 1.3 0.2 1; 4.9 3.6 1.4 0.1 1; 4.4 3.0 1.3 0.2 1;
5.1 3.4 1.5 0.2 1; 5.0 3.5 1.3 0.3 1; 4.5 2.3 1.3 0.3 1;
4.4 3.2 1.3 0.2 1; 5.0 3.5 1.6 0.6 1; 5.1 3.8 1.9 0.4 1;
4.8 3.0 1.4 0.3 1; 5.1 3.8 1.6 0.2 1; 4.6 3.2 1.4 0.2 1;
5.3 3.7 1.5 0.2 1; 5.0 3.3 1.4 0.2 1;
7.0 3.2 4.7 1.4 2; 6.4 3.2 4.5 1.5 2; 6.9 3.1 4.9 1.5 2;
5.5 2.3 4.0 1.3 2; 6.5 2.8 4.6 1.5 2; 5.7 2.8 4.5 1.3 2;
6.3 3.3 4.7 1.6 2; 4.9 2.4 3.3 1.0 2; 6.6 2.9 4.6 1.3 2;
5.2 2.7 3.9 1.4 2; 5.0 2.0 3.5 1.0 2; 5.9 3.0 4.2 1.5 2;
6.0 2.2 4.0 1.0 2; 6.1 2.9 4.7 1.4 2; 5.6 2.9 3.6 1.3 2;
6.7 3.1 4.4 1.4 2; 5.6 3.0 4.5 1.5 2; 5.8 2.7 4.1 1.0 2;
6.2 2.2 4.5 1.5 2; 5.6 2.5 3.9 1.1 2; 5.9 3.2 4.8 1.8 2;
6.1 2.8 4.0 1.3 2; 6.3 2.5 4.9 1.5 2; 6.1 2.8 4.7 1.2 2;
6.4 2.9 4.3 1.3 2; 6.6 3.0 4.4 1.4 2; 6.8 2.8 4.8 1.4 2;
6.7 3.0 5.0 1.7 2; 6.0 2.9 4.5 1.5 2; 5.7 2.6 3.5 1.0 2;
5.5 2.4 3.8 1.1 2; 5.5 2.4 3.7 1.0 2; 5.8 2.7 3.9 1.2 2;
6.0 2.7 5.1 1.6 2; 5.4 3.0 4.5 1.5 2; 6.0 3.4 4.5 1.6 2;
6.7 3.1 4.7 1.5 2; 6.3 2.3 4.4 1.3 2; 5.6 3.0 4.1 1.3 2;
5.5 2.5 4.0 1.3 2; 5.5 2.6 4.4 1.2 2; 6.1 3.0 4.6 1.4 2;
5.8 2.6 4.0 1.2 2; 5.0 2.3 3.3 1.0 2; 5.6 2.7 4.2 1.3 2;
5.7 3.0 4.2 1.2 2; 5.7 2.9 4.2 1.3 2; 6.2 2.9 4.3 1.3 2;
5.1 2.5 3.0 1.1 2; 5.7 2.8 4.1 1.3 2;
6.3 3.3 6.0 2.5 3; 5.8 2.7 5.1 1.9 3; 7.1 3.0 5.9 2.1 3;
6.3 2.9 5.6 1.8 3; 6.5 3.0 5.8 2.2 3; 7.6 3.0 6.6 2.1 3;
4.9 2.5 4.5 1.7 3; 7.3 2.9 6.3 1.8 3; 6.7 2.5 5.8 1.8 3;
7.2 3.6 6.1 2.5 3; 6.5 3.2 5.1 2.0 3; 6.4 2.7 5.3 1.9 3;
6.8 3.0 5.5 2.1 3; 5.7 2.5 5.0 2.0 3; 5.8 2.8 5.1 2.4 3;
6.4 3.2 5.3 2.3 3; 6.5 3.0 5.5 1.8 3; 7.7 3.8 6.7 2.2 3;
7.7 2.6 6.9 2.3 3; 6.0 2.2 5.0 1.5 3; 6.9 3.2 5.7 2.3 3;
5.6 2.8 4.9 2.0 3; 7.7 2.8 6.7 2.0 3; 6.3 2.7 4.9 1.8 3;
6.7 3.3 5.7 2.1 3; 7.2 3.2 6.0 1.8 3; 6.2 2.8 4.8 1.8 3;
6.1 3.0 4.9 1.8 3; 6.4 2.8 5.6 2.1 3; 7.2 3.0 5.8 1.6 3;
7.4 2.8 6.1 1.9 3; 7.9 3.8 6.4 2.0 3; 6.4 2.8 5.6 2.2 3;
6.3 2.8 5.1 1.5 3; 6.1 2.6 5.6 1.4 3; 7.7 3.0 6.1 2.3 3;
6.3 3.4 5.6 2.4 3; 6.4 3.1 5.5 1.8 3; 6.0 3.0 4.8 1.8 3;
6.9 3.1 5.4 2.1 3; 6.7 3.1 5.6 2.4 3; 6.9 3.1 5.1 2.3 3;
5.8 2.7 5.1 1.9 3; 6.8 3.2 5.9 2.3 3; 6.7 3.3 5.7 2.5 3;
6.7 3.0 5.2 2.3 3; 6.3 2.5 5.0 1.9 3; 6.5 3.0 5.2 2.0 3;
6.2 3.4 5.4 2.3 3; 5.9 3.0 5.1 1.8 3];

data = raw;
end
