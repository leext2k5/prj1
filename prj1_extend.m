clc;
clearvars;

%% ================== THAM SỐ HỆ THỐNG ==================
Kmax = 10;                 % số sub-channel tối đa
N = 120;                   % số mẫu
sigmaW2 = 1;               % nhiễu Willie
tau = 1.05 * sigmaW2;      % ngưỡng phát hiện

%% ================== CÔNG SUẤT PHÁT ==================
Pa = [ ...
    0.05, ...
    0.09, ...
    0.04, ...
    0.08, ...
    0.03, ...
    0.07, ...
    0.025, ...
    0.06, ...
    0.015, ...
    0.1 ];

%% ================== HÀM TOÁN ==================
gamma_inc = @(a,x) gammainc(x, a, 'lower');

%% ================== MẢNG KẾT QUẢ ==================
P_fa = zeros(Kmax, Kmax);   % P_fa(M,k)
P_md = zeros(Kmax, Kmax);   % P_md(M,k)
P_E  = zeros(Kmax, Kmax);   % P_E(M,k)

%% ================== TÍNH TOÁN (THEO PAPER) ==================
for k = 1:Kmax
    Pa_sorted = sort(Pa(1:k), 'descend');   % chọn kênh mạnh nhất
    % ---- Eq.(13): False Alarm ----
    P_fa(:,k) = 1 - (gamma_inc(N, N*tau/sigmaW2))^k;

    for M = 1:k
        % ---- Eq.(14): Missed Detection (RSCT) ----
        term_signal = prod( ...
            gamma_inc(N, N*tau ./ (Pa_sorted(1:M) + sigmaW2)) );

        term_noise = gamma_inc(N, N*tau/sigmaW2)^(k - M);

        P_md(M,k) = term_signal * term_noise;

        % ---- Eq.(16): Detection Error ----
        P_E(M,k) = P_fa(M,k) + P_md(M,k);
    end
end

K = 1:Kmax;
colors = lines(Kmax);
%% ================== VẼ P_FA ==================
figure; hold on;
for M = 1:5
    plot(M:Kmax, P_fa(M, M:Kmax), '-^', ...
        'LineWidth', 1.8, ...
        'Color', colors(M,:));
end
xlabel('Number of sub-channels K');
ylabel('P_{FA}');
legend(arrayfun(@(x) sprintf('M = %d', x), 1:5, 'UniformOutput', false), ...
       'Location', 'best');
grid on;

%% ================== VẼ P_MD ==================
figure; hold on;
for M = 1:5
    plot(M:Kmax, P_md(M, M:Kmax), '-o', ...
        'LineWidth', 1.8, ...
        'Color', colors(M,:));
end
xlabel('Number of sub-channels K');
ylabel('P_{MD}');
legend(arrayfun(@(x) sprintf('M = %d', x), 1:5, 'UniformOutput', false), ...
       'Location', 'best');
grid on;

%% ================== VẼ P_E ==================
figure; hold on;
for M = 1:5
    plot(M:Kmax, P_E(M, M:Kmax), '-s', ...
        'LineWidth', 1.8, ...
        'Color', colors(M,:));
end
xlabel('Number of sub-channels K');
ylabel('P_{E}');
legend(arrayfun(@(x) sprintf('M = %d', x), 1:5, 'UniformOutput', false), ...
       'Location', 'best');
grid on;
