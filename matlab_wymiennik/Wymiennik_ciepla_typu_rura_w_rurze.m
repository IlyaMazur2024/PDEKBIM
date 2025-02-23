
clear all 
clc

% Wartości parametrów modelu
Asp = 0.2; Apt = 0.226; Atc = 0.327; Aco = 0.352; % m^2
cs = 4185; ct = 4185; cp = 470; cc = 470; % J/kg/K
ms = 0.8; mp = 0.54;mt = 1.1; mc = 0.85; % kg
alpha_sp = 1000; alpha_pt = 1000; alpha_tc = 1000; alpha_co = 10; % W/m^2/K

Tsi = 90; Tti = 30; To = 20; % °C
Ts0 = 90; Tp0 = 70; Tt0 = 50; Tc0 = 20; % Początkowe temperatury °C

% Definiowanie przepływów masowych i ich wartości
m_dot_s = @(t) 0.05; % przepływ masowy m_s(t)
m_dot_t = @(t) 0.01; % przepływ masowy m_t(t)

% Definicja układu równań różniczkowych
ode_system = @(t, T) [
    (1 / (ms * cs)) * (m_dot_s(t) * cs * (Tsi - T(1)) - alpha_sp * Asp * (T(1) - T(2))); % dTs/dt
    (1 / (mp * cp)) * (alpha_sp * Asp * (T(1) - T(2)) - alpha_pt * Apt * (T(2) - T(3))); % dTp/dt
    (1 / (mt * ct)) * (alpha_pt * Apt * (T(2) - T(3)) - alpha_tc * Atc * (T(3) - T(4)) - m_dot_t(t) * ct * (T(3) - Tti)); % dTt/dt
    (1 / (mc * cc)) * (alpha_tc * Atc * (T(3) - T(4)) - alpha_co * Aco * (T(4) - To))  % dTc/dt
];

% Warunki początkowe
T0 = [Ts0; Tp0; Tt0; Tc0];
% Rozwiązywanie układu równań
t_span = [0, 1000]; 
h = 0.1;           
[t, T] = my_ode45(ode_system, t_span, T0, h);


% Rysowanie wyników
figure;
plot(t, T);
grid on;
legend('Ts(t) - temperatura czynnika ogrewającego', 'Tp(t) - temperatura ścianki wewnętrznej', 'Tt(t) - temperatura czynnika ogrzewanego', 'Tc(t) - temperatura ścianki zewnętrznej', 'Location', 'Best');
xlabel('Czas (s)');
ylabel('Temperatura (°C)');
title('Rozwiązanie układu równań różniczkowych');

%%
mdot_s0 = 0.05;  
mdot_t0 = 0.01;  

% macierze
A = [(-alpha_sp*Asp - mdot_s0*cs)/(ms*cs), alpha_sp*Asp/(ms*cs), 0, 0;
     alpha_sp*Asp/(mp*cp), (-alpha_sp*Asp - alpha_pt*Apt)/(mp*cp), alpha_pt*Apt/(mp*cp), 0;
     0, alpha_pt*Apt/(mt*ct), (-alpha_pt*Apt - alpha_tc*Atc - mdot_t0*ct)/(mt*ct), alpha_tc*Atc/(mt*ct);
     0, 0, alpha_tc*Atc/(mc*cc), (-alpha_tc*Atc - alpha_co*Aco)/(mc*cc)];

B = [(Tsi/ms), 0;
     0, 0;
     0, (Tti/mt);
     0, (alpha_co*Aco)/(mc*cc)*To];

C = eye(4);
D = zeros(4,2);

% utworzenie modelu
sys = ss(A, B, C, D);

% czas symulacji
t = 0:0.1:500;

% wejścia
mdot_s = 0.05*ones(size(t)); % 
mdot_t = 0.01*ones(size(t));  % 
u = [mdot_s; mdot_t]'; 


% symulacja systemu
[y, t, x] = lsim(sys, u, t, T0);

% wykresy
figure

plot(t, y(:,1), 'r', t, y(:,2), 'g', t, y(:,3), 'b', t, y(:,4), 'm');
grid on
legend('Ts(t) - temperatura czynnika ogrewającego', 'Tp(t) - temperatura ścianki wewnętrznej', 'Tt(t) - temperatura czynnika ogrzewanego', 'Tc(t) - temperatura ścianki zewnętrznej', 'Location', 'Best');
xlabel('Czas (s)');
ylabel('Temperatura (°C)');
title('Rozwiązanie układu równań stanu');

% figure
% plot(t, u(:,1), 'k--', t, u(:,2),'c--');
% xlabel('Czas (s)');
% ylabel('Przepływ masowy Kg/s');
% title('Sygnały wejściowe przepływów masowych');
% legend('mdot_s(t) - przepływ masowy czynnika ogrzewającego', 'mdot_t - przepływ masowy czynnika ogrzewanego');
% grid on;

disp('Matrix A:'); disp(A);
disp('Matrix B:'); disp(B);
disp('Matrix C:'); disp(C);
disp('Matrix D:'); disp(D);
%%

% ilosc we/wy
n_inputs = size(B, 2);
n_outputs = size(C, 1);

% konwestja ssr do tf
for j = 1:n_inputs 
    [num, den] = ss2tf(A, B, C, D, j); 
    fprintf('TF od wejścia %d:\n', j);
    for i = 1:n_outputs 
        tf_ij = tf(num(i,:), den);
        fprintf('  wyjście %d: ', i);
        disp(tf_ij)
    end
    fprintf('\n') 
end

s = tf('s'); 
disp('TF symbolicznie:')
H = C * inv(s*eye(size(A)) - A) * B + D


%% Charekterystyki częstotliwościowe
[mag,phase,wout] = bode(sys)
[re,im,wout] = nyquist(sys)
figure
bode(sys)
figure
nyquist(sys)
%% Transmitancje

s = tf('s');

% Transmitancja dla T_s(s)
num_T_s = (cs * Tsi) / (ms * cs) * m_dot_s(s) + (alpha_sp * Asp) / (ms * cs) * T_p(s);
den_T_s = s + (cs / (ms * cs)) * m_dot_s(s) + (alpha_sp * Asp) / (ms * cs);
T_s(s) = num_T_s / den_T_s;

% Transmitancja dla T_p(s)
num_T_p = (alpha_sp * Asp) / (mp * cp) * T_s(s) + (alpha_pt * Apt) / (mp * cp) * T_t(s);
den_T_p = s + (alpha_sp * Asp) / (mp * cp) + (alpha_pt * Apt) / (mp * cp);
T_p(s) = num_T_p / den_T_p;

% Transmitancja dla T_t(s)
num_T_t = (alpha_pt * Apt) / (mt * ct) * T_p(s) + (alpha_tc * Atc) / (mt * ct) * T_c(s) + (ct * Tti) / (mt * ct) * m_dot_t(s);
den_T_t = s + (alpha_pt * Apt) / (mt * ct) + (alpha_tc * Atc) / (mt * ct) + (ct / (mt * ct)) * m_dot_t(s);
T_t(s) = num_T_t / den_T_t;

% Transmitancja dla T_c(s)
num_T_c = (alpha_tc * Atc) / (mc * cc) * T_t(s) + (alpha_co * Aco) / (mc * cc) * To;
den_T_c = s + (alpha_tc * Atc) / (mc * cc) + (alpha_co * Aco) / (mc * cc);
T_c(s) = num_T_c / den_T_c;

%% Charakterystyki statyczne
% Wartości przepływu m_dot_t do analizy
m_dot_t_values = [0.005, 0.015]; 

% Zakres wartości m_dot_s
m_dot_s_values = linspace(0, 0.1, 25); 

% Przechowywanie wyników
Tt_steady = zeros(length(m_dot_s_values), length(m_dot_t_values));

for j = 1:length(m_dot_t_values)
    m_dot_t = @(t) m_dot_t_values(j); % Ustawienie wartości m_dot_t

    for i = 1:length(m_dot_s_values)
        m_dot_s = @(t) m_dot_s_values(i); % Ustawienie wartości m_dot_s
        
        % Rozwiązywanie układu równań różniczkowych dla danej kombinacji
        [t, T] = my_ode45(ode_system, t_span, T0, h);
        
        % Pobranie wartości ustalonej Tt (ostatnia wartość z rozwiązania)
        Tt_steady(i, j) = T(end, 3);
    end
end

% Rysowanie charakterystyk statycznych
figure; hold on;
plot(m_dot_s_values, Tt_steady(:,1), 'ro-', 'DisplayName', 'mdot_s = 0.005 kg/s');
plot(m_dot_s_values, Tt_steady(:,2), 'go-', 'DisplayName', 'mdot_s = 0.015 kg/s');

xlabel('Przepływ masowy [kg/s]');
ylabel('Temperatura w stanie ustalonym [°C]');
title('Charakterystyki statyczne T_t w zależności od mdot_s');
legend;
grid on;
hold off;


%% Pomiary charakterystyk częstotliwościowych
omegas = linspace(0.001, 1, 15);  % lub można użyć logspace(np. logspace(-3,0,15))

for k = 1:length(omegas)
    omega = omegas(k);
    m_dot_s = @(t) 0.05 + 0.01*sin(omega*t);
    m_dot_t = @(t) 0.01;
    

    ode_system = @(t, T) [
        (1/(ms*cs)) * ( m_dot_s(t)*cs*(Tsi - T(1)) - alpha_sp*Asp*(T(1) - T(2)) );
        (1/(mp*cp)) * ( alpha_sp*Asp*(T(1) - T(2)) - alpha_pt*Apt*(T(2) - T(3)) );
        (1/(mt*ct)) * ( alpha_pt*Apt*(T(2) - T(3)) - alpha_tc*Atc*(T(3) - T(4)) - m_dot_t(t)*ct*(T(3) - Tti) );
        (1/(mc*cc)) * ( alpha_tc*Atc*(T(3) - T(4)) - alpha_co*Aco*(T(4) - To) )
    ];
    
    t_span = [0, 500];
    h = 0.1;
    
    [t_sim, T_sim] = my_ode45(ode_system, t_span, T0, h);
    
    m_dot_s_signal = arrayfun(m_dot_s, t_sim);
    
    % Rysowanie przebiegów sygnałów
    figure;
    % Sygnał wejściowy: m_dot_s(t)
    subplot(2,1,1);
    plot(t_sim, m_dot_s_signal, 'b-', 'LineWidth', 1.5);
    xlabel('Czas [s]');
    ylabel('Przepływ masowy [kg/s]');
    title(['Sygnał wejściowy \dot{m}_s(t) dla \omega = ', num2str(omega), ' rad/s']);
    grid on;
    
    % Odpowiedź wyjściowa: T_t(t) (trzeci element wektora stanu)
    subplot(2,1,2);
    plot(t_sim, T_sim(:,3), 'r-', 'LineWidth', 1.5);
    xlabel('Czas [s]');
    ylabel('Temperatura T_t [°C]');
    title(['Odpowiedź wyjściowa T_t(t) dla \omega = ', num2str(omega), ' rad/s']);
    grid on;
    pause
end
%%