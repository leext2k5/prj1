clc;
clearvars;

global Pmax sigmaW2 Kmax D epsilon sigmaB2
Pmax = 10^17;
Pa = 10^-4.4;

sigmaW2 = 10^-3.4; 
sigmaB2 = 10^-3.4; 

Kmax = 10;
K = 1:Kmax;
N_list = [120, 140, 160];
D = 12;
epsilon = 0.003;

%Hàm incomplete gamma
gamma_inc = @(a,x) gammainc(x,a,'lower') * gamma(a);

%Hàm Q(x)
Q = @(x) 0.5*erfc(x./sqrt(2));

%Mảng kết quả
tau_star = zeros(length(N_list), Kmax);
PE_star = zeros(length(N_list), Kmax);
Pa_star = zeros(length(N_list), Kmax);
eta = zeros(length(N_list), Kmax);



for n = 1:length(N_list)
    N = N_list(n); 
    R = D/N;

  for k = 1:Kmax
    f = @(tau) tau ...
        - ((Pa + sigmaW2)*sigmaW2 / Pa) * log((Pa + sigmaW2)/sigmaW2) ...
        - ((Pa + sigmaW2)*sigmaW2 / (N * Pa)) * ...
        log(k - (k-1)*gamma_inc(N, N*tau/(Pa + sigmaW2)) / gamma_inc(N, N*tau/sigmaW2));

    %tìm tau*
    tau0 = 0.1;
    tau_star(n,k) = fsolve(f, tau0);
    %disp([tau_star(n,k)]);
    
    %Tính P_E*
    A = gamma_inc(N, N*tau_star(n,k)/sigmaW2)/gamma(N);
    B = gamma_inc(N, N*tau_star(n,k)/(Pa+sigmaW2))/gamma(N);
    PE_star(n, k) = 1 - A^k*(1 - B/A);




    %tìm P_epsilon
    g = @(P_epsilon) 1 - (gamma_inc(N, N*tau_star(n,k)/sigmaW2)/gamma(N))^k...
        *(1 - gamma_inc(N, N*tau_star(n,k)/(P_epsilon+sigmaW2))/gamma_inc(N, N*tau_star(n,k)/sigmaW2))...
        -(1 - epsilon);

    P0 = 0.00000001;
    P_epsilon0 = fsolve(g, P0);
    disp(P_epsilon0);
  %Tìm Pa_star
    if(P_epsilon0 > Pmax)
        Pa_star(n, k) = Pmax;
    else
        Pa_star(n, k) = P_epsilon0;
    end
   %Tính eta
     SNR = Pa_star(n,k) / sigmaB2;
     DELTA = Q((sqrt(N)*(SNR + 1)*log(SNR+1)-R)/sqrt(SNR*(SNR+2)));
     eta(n,k) = 3000*R*(1-DELTA)/k;

   end
end

%Vẽ P_E*
figure;
plot(K, PE_star(1,:), 'b--o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(K, PE_star(2,:), 'r-o', 'LineWidth', 2, 'MarkerSize', 8); 
plot(K, PE_star(3,:), 'k-*', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('The number of sub-channels');
ylabel('Minimum detection error probability');
legend('N = 120', 'N = 140', 'N = 160');
grid on; 


%Vẽ eta
figure;
plot(K, eta(1,:), 'b--o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(K, eta(2,:), 'r-o', 'LineWidth', 2, 'MarkerSize', 8);  
plot(K, eta(3,:), 'k-*', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('The number of sub-channels K');
ylabel('Average effective covertness rate η');
legend('N = 120', 'N = 140', 'N = 160');
grid on;