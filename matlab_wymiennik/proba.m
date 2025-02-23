
clc; clear; close all;

% Parametry symulacji
L = 4;         % Długość rury [m]
Nx = 50;        % Liczba podziałów wzdłuż rury
dx = L / Nx;    % Krok przestrzenny
dt = 0.01;      % Krok czasowy
Tmax = 500;      % Czas symulacji [s]
Nt = Tmax / dt; % Liczba kroków czasowych

% Prędkości przepływu
vt = 1.0;  % Prędkość płynu wewnętrznego [m/s]
vs = 0.5;  % Prędkość płynu zewnętrznego [m/s]

% Parameters for heat exchange 

% Geometry
d1 = 16e-3; % Inner tube diameter (m)
d2 = 18e-3; % Outer wall diameter (m)
d3 = 26e-3; % Shell internal diameter (m)
d4 = 28e-3; % Shell external diameter (m)

% Fluid and material properties
rho_t = 1000; % Inner fluid density (kg/m^3) 
rho_s = rho_t;% Shell-side fluid density (kg/m^3)
rho_w = 7800; % Wall material density (kg/m^3) 
rho_z = rho_w;% Ambient density (kg/m^3)

c_t = 4200; % Inner fluid specific heat (J/kg*K)
c_s = c_t;  % Shell-side fluid specific heat (J/kg*K)
c_w = 500;  % Wall material specific heat (J/kg*K)
c_z = c_w;  % Ambient specific heat (J/kg*K)

h_t = 2000; % Heat transfer coefficient for inner fluid (W/m^2*K)
h_s = h_t;  % Heat transfer coefficient for shell-side fluid (W/m^2*K)
h_z = 5;    % Heat transfer coefficient for ambient (W/m^2*K)

% Calculate coefficients
k1 = (4 * h_t) / (rho_t * c_t * d1);
k2 = (4 * d1 * h_t) / (rho_w * c_w * (d2^2 - d1^2));
k3 = (4 * d2 * h_s) / (rho_w * c_w * (d2^2 - d1^2));
k4 = (4 * d2 * h_s) / (rho_s * c_s * (d3^2 - d2^2));
k5 = (4 * d3 * h_s) / (rho_s * c_s * (d3^2 - d2^2));
k6 = (4 * d3 * h_s) / (rho_z * c_z * (d4^2 - d3^2));
k7 = (4 * h_z) / (rho_z * c_z * d4);

% Warunki początkowe
Tt = 80 * ones(Nx,1); % Temperatura płynu wewnętrznego
Tw = 50 * ones(Nx,1); % Temperatura ścianki rury
Ts = 30 * ones(Nx,1); % Temperatura płynu zewnętrznego
Tz = 20 * ones(Nx,1); % Temperatura otoczenia
Ta = 20;              % Temperatura zewnętrzna (stała)

% Macierze wyników
Tt_hist = zeros(Nx, Nt);
Tw_hist = zeros(Nx, Nt);
Ts_hist = zeros(Nx, Nt);
Tz_hist = zeros(Nx, Nt);

% Symulacja numeryczna
for n = 1:Nt
    % Zachowanie poprzednich wartości
    Tt_old = Tt;
    Tw_old = Tw;
    Ts_old = Ts;
    Tz_old = Tz;

    % Aktualizacja temperatury płynu wewnętrznego 
    for i = 2:Nx
        Tt(i) = Tt_old(i) - dt * vt * (Tt_old(i) - Tt_old(i-1)) / dx ...
                + dt * k1 * (Tw_old(i) - Tt_old(i));
    end

    % Aktualizacja temperatury ścianki rury 
    Tw = Tw_old + dt * (k2 * (Tt_old - Tw_old) + k3 * (Ts_old - Tw_old));

    % Aktualizacja temperatury płynu zewnętrznego 
    for i = 2:Nx
        Ts(i) = Ts_old(i) - dt * vs * (Ts_old(i) - Ts_old(i-1)) / dx ...
                + dt * (k4 * (Tw_old(i) - Ts_old(i)) + k5 * (Tz_old(i) - Ts_old(i)));
    end

    % Aktualizacja temperatury otoczenia 
    Tz = Tz_old + dt * (k6 * (Ts_old - Tz_old) + k7 * (Ta - Tz_old));

    % Zapis wyników do historii
    Tt_hist(:, n) = Tt;
    Tw_hist(:, n) = Tw;
    Ts_hist(:, n) = Ts;
    Tz_hist(:, n) = Tz;
end

% Wizualizacja wyników
x = linspace(0, L, Nx);
t = linspace(0, Tmax, Nt);

figure;
subplot(2,2,1);
imagesc(t, x, Tt_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura T_t');

subplot(2,2,2);
imagesc(t, x, Tw_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura T_w');

subplot(2,2,3);
imagesc(t, x, Ts_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura T_s');

subplot(2,2,4);
imagesc(t, x, Tz_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura T_z');
%% przeciwzbiezny

clc; clear; close all;

% Parametry symulacji
L = 4;         % Długość rury [m]
Nx = 50;        % Liczba podziałów wzdłuż rury
dx = L / Nx;    % Krok przestrzenny
dt = 0.01;      % Krok czasowy
Tmax = 500;      % Czas symulacji [s]
Nt = Tmax / dt; % Liczba kroków czasowych

% Prędkości przepływu
vt = 1.0;  % Prędkość płynu wewnętrznego [m/s]
vs = -0.5;  % Prędkość płynu zewnętrznego [m/s]

% Parameters for heat exchange 

% Given data for inner and outer pipes
% Geometry
d1 = 16e-3; % Inner tube diameter (m)
d2 = 18e-3; % Outer wall diameter (m)
d3 = 26e-3; % Shell internal diameter (m)
d4 = 28e-3; % Shell external diameter (m)

% Fluid and material properties
rho_t = 1000; % Inner fluid density (kg/m^3) 
rho_s = rho_t;% Shell-side fluid density (kg/m^3)
rho_w = 7800; % Wall material density (kg/m^3) 
rho_z = rho_w;% Ambient density (kg/m^3)

c_t = 4200; % Inner fluid specific heat (J/kg*K)
c_s = c_t;  % Shell-side fluid specific heat (J/kg*K)
c_w = 500;  % Wall material specific heat (J/kg*K)
c_z = c_w;  % Ambient specific heat (J/kg*K)

h_t = 2000; % Heat transfer coefficient for inner fluid (W/m^2*K)
h_s = h_t;  % Heat transfer coefficient for shell-side fluid (W/m^2*K)
h_z = 5;    % Heat transfer coefficient for ambient (W/m^2*K)

% Obliczenie współczynników k
k1 = (4 * h_t) / (rho_t * c_t * d1);
k2 = (4 * d1 * h_t) / (rho_w * c_w * (d2^2 - d1^2));
k3 = (4 * d2 * h_s) / (rho_w * c_w * (d2^2 - d1^2));
k4 = (4 * d2 * h_s) / (rho_s * c_s * (d3^2 - d2^2));
k5 = (4 * d3 * h_s) / (rho_s * c_s * (d3^2 - d2^2));
k6 = (4 * d3 * h_s) / (rho_z * c_z * (d4^2 - d3^2));
k7 = (4 * h_z) / (rho_z * c_z * d4);


% Warunki początkowe
Tt = 80 * ones(Nx,1); % Temperatura płynu wewnętrznego (wchodzi na lewym końcu)
Tw = 50 * ones(Nx,1); % Temperatura ścianki rury
Ts = 30 * ones(Nx,1); % Temperatura płynu zewnętrznego (wchodzi na prawym końcu)
Tz = 20 * ones(Nx,1); % Temperatura otoczenia
Ta = 20;              % Temperatura zewnętrzna (stała)

% **Warunki brzegowe**
Tt_in = 80; % Temperatura wejściowa Tt (wejście na lewej krawędzi)
Ts_in = 30; % Temperatura wejściowa Ts (wejście na prawej krawędzi)

% **Macierze wyników**
Tt_hist = zeros(Nx, Nt);
Tw_hist = zeros(Nx, Nt);
Ts_hist = zeros(Nx, Nt);
Tz_hist = zeros(Nx, Nt);

% **Symulacja numeryczna**
for n = 1:Nt
  
    Tt_old = Tt;
    Tw_old = Tw;
    Ts_old = Ts;
    Tz_old = Tz;

    % **Warunek brzegowy dla Tt**
    Tt(1) = Tt_in;

    %  temperatury płynu wewnętrznego (konwekcja + wymiana ciepła)
    for i = 2:Nx
        Tt(i) = Tt_old(i) - (dt / dx) * vt * (Tt_old(i) - Tt_old(i-1)) ...
                + dt * k1 * (Tw_old(i) - Tt_old(i));
    end

    %temperatury ścianki rury (wymiana ciepła)**
    Tw = Tw_old + dt * (k2 * (Tt_old - Tw_old) + k3 * (Ts_old - Tw_old));

    % **Warunek brzegowy dla Ts (wejście na końcu rury!)**
    Ts(Nx) = Ts_in;

    %temperatury płynu zewnętrznego (przeciwprąd)
    for i = Nx-1:-1:1  % Przepływ przeciwny → zaczynamy od prawej
        Ts(i) = Ts_old(i) - (dt / dx) * vs * (Ts_old(i+1) - Ts_old(i)) ...
                + dt * (k4 * (Tw_old(i) - Ts_old(i)) + k5 * (Tz_old(i) - Ts_old(i)));
    end

    % **Aktualizacja temperatury otoczenia**
    Tz = Tz_old + dt * (k6 * (Ts_old - Tz_old) + k7 * (Ta - Tz_old));

    % **Zapis wyników**
    Tt_hist(:, n) = Tt;
    Tw_hist(:, n) = Tw;
    Ts_hist(:, n) = Ts;
    Tz_hist(:, n) = Tz;
end

% **Wizualizacja wyników**
x = linspace(0, L, Nx);
t = linspace(0, Tmax, Nt);

figure;
subplot(2,2,1);
imagesc(t, x, Tt_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura Tt');

subplot(2,2,2);
imagesc(t, x, Tw_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura Tw');

subplot(2,2,3);
imagesc(t, x, Ts_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura Ts (przeciwprąd)');

subplot(2,2,4);
imagesc(t, x, Tz_hist);
colorbar; xlabel('Czas [s]'); ylabel('Pozycja [m]'); title('Temperatura Tz');
