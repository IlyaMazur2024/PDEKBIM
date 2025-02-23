%% Czyszczenie środowiska
clear; close all; clc;

%% Parametry modelu
Asp = 0.2;    Apt = 0.226;   Atc = 0.327;   Aco = 0.352;    % [m^2]
cs  = 4185;   ct  = 4185;    cp  = 470;     cc  = 470;      % [J/kg/K]
ms  = 0.8;    mp  = 0.54;     mt  = 1.1;     mc  = 0.85;     % [kg]
alpha_sp = 1000;  alpha_pt = 1000;  alpha_tc = 1000;  alpha_co = 10; % [W/m^2/K]

Tsi = 90;     Tti = 30;      To = 20;      % [°C]
Ts0 = 90;     Tp0 = 70;      Tt0 = 50;     Tc0 = 20;     % Warunki początkowe [°C]
T0 = [Ts0; Tp0; Tt0; Tc0];

%% Parametry symulacji
dt = 0.1;                % krok czasowy [s]
T_end = 4000;            % czas symulacji [s]
time = 0:dt:T_end;
N = length(time);

%% Parametry regulatora PID
% Dostosuj współczynniki PID do dynamiki układu
Kp = 0.01;
Ki = 0.0001;
Kd = 0.005;



% Inicjalizacja zmiennych PID
error_int = 0;                   % całka błędu
error_prev = 50 - Tt0;            % początkowy błąd (dla setpointu 50°C)

%% Inicjalizacja wektorów wynikowych
T_history = zeros(N, 4);         % stany: [Ts, Tp, Tt, Tc]
T_history(1,:) = T0';
u_history = zeros(N, 1);         % sygnał sterujący: m_dot_s
Tt_ref_history = zeros(N, 1);      % historia setpointu dla T_t
m_dot_t_history = zeros(N, 1);     % historia zakłócenia m_dot_t

%% Parametr dla zakłócenia sinusoidalnego
omega_d = 2*pi/1000;  % pulsacja zapewniająca okres oscylacji 1000 s

%% Symulacja układu z regulatorem PID
for i = 1:N-1
    t_current = time(i);
    
    % Ustalanie setpointu T_t (wartość zadana)
    if t_current < 1000
        Tt_ref = 50;
    elseif t_current < 2000
        Tt_ref = 60;
    else
        Tt_ref = 70;
    end
    Tt_ref_history(i) = Tt_ref;
    
    % Ustalanie zakłócenia: sinusoidalnie zmieniające się m_dot_t
    if t_current < 1000
         m_dot_t_val = 0.005 + 0.005 * sin(omega_d * t_current);
    elseif t_current < 2000
        m_dot_t_val = 0.01 + 0.005 * sin(omega_d * t_current);
    else
        m_dot_t_val = 0.001 + 0.005 * sin(omega_d * t_current);
    end
   
    m_dot_t_history(i) = m_dot_t_val;
    
    % Aktualny stan układu
    T_current = T_history(i,:)';
    
    % Mierzona wartość: T_t (trzeci stan)
    Tt_meas = T_current(3);
    
    % Obliczenie błędu regulacji
    error = Tt_ref - Tt_meas;
    
    % Składowe regulatora PID
    error_int = error_int + error * dt;           % składnik całkujący
    error_der = (error - error_prev) / dt;          % składnik różniczkujący
    
    % Obliczenie sygnału sterującego (PID)
    u = Kp * error + Ki * error_int + Kd * error_der;
    
    % Ograniczenie sygnału sterującego do przedziału [0, 0.1] kg/s
    u = max(0, min(u, 0.1));
    u_history(i) = u;
    
    % Obliczenie pochodnych stanu według modelu
    % Stany: T = [Ts; Tp; Tt; Tc]
    Ts = T_current(1);  
    Tp = T_current(2);  
    Tt_val = T_current(3);  
    Tc = T_current(4);
    
    dTs_dt = (1/(ms*cs)) * ( u*cs*(Tsi - Ts) - alpha_sp*Asp*(Ts - Tp) );
    dTp_dt = (1/(mp*cp)) * ( alpha_sp*Asp*(Ts - Tp) - alpha_pt*Apt*(Tp - Tt_val) );
    dTt_dt = (1/(mt*ct)) * ( alpha_pt*Apt*(Tp - Tt_val) - alpha_tc*Atc*(Tt_val - Tc) - m_dot_t_val*ct*(Tt_val - Tti) );
    dTc_dt = (1/(mc*cc)) * ( alpha_tc*Atc*(Tt_val - Tc) - alpha_co*Aco*(Tc - To) );
    
    % Aktualizacja stanu (metoda Eulera)
    T_next = T_current + dt * [dTs_dt; dTp_dt; dTt_dt; dTc_dt];
    T_history(i+1,:) = T_next';
    
    % Aktualizacja poprzedniego błędu
    error_prev = error;
end
% Uzupełnienie ostatnich wartości historii setpointu i zakłócenia
Tt_ref_history(end) = Tt_ref_history(end-1);
m_dot_t_history(end) = m_dot_t_history(end-1);

%% Wizualizacja wyników
figure;

% Przebieg T_t oraz zmieniający się setpoint
subplot(3,1,1);
plot(time, T_history(:,3), 'b', 'LineWidth', 1.5);
hold on;
plot(time, Tt_ref_history, 'r--', 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('T_t [°C]');
title('Odpowiedź wyjściowa T_t(t) i wartość zadana');
legend('T_t','Setpoint');
grid on;

% Sygnał sterujący: m_dot_s
subplot(3,1,2);
plot(time, u_history, 'g', 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('\dot{m}_s [kg/s]');
title('Sygnał sterujący \dot{m}_s(t)');
grid on;

% Zakłócenie: m_dot_t
subplot(3,1,3);
plot(time, m_dot_t_history, 'm', 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('\dot{m}_t [kg/s]');
title('Zakłócenie \dot{m}_t(t) - sinusoidalne');
grid on;
