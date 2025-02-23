function [t, T] = my_ode45(odefun, tspan, T0, h)
    % my_ode45 - Rozwiązywanie układu ODE metodą Rungego-Kutty 4 rzędu.
    %
    % Wejście:
    % odefun - funkcja definiująca układ ODE (t, T) -> dT/dt
    % tspan  - zakres czasu [t_start, t_end]
    % T0     - warunki początkowe (wektor)
    % h      - krok czasowy
    %
    % Wyjście:
    % t      - wektor czasu
    % T      - macierz wyników (każda kolumna odpowiada stanowi w danym czasie)
    
    % Inicjalizacja
    t_start = tspan(1);
    t_end = tspan(2);
    t = t_start:h:t_end; % Dyskretyzacja czasu
    n_steps = length(t);
    n_vars = length(T0); % Liczba zmiennych w układzie ODE
    T = zeros(n_steps, n_vars);
    
    % Ustawienie warunków początkowych
    T(1, :) = T0;
    
    % Pętla RK4
    for i = 1:n_steps-1
        ti = t(i);        % Aktualny czas
        Ti = T(i, :)';    % Aktualny stan (jako wektor kolumnowy)
        
        % RK4 kroki
        k1 = odefun(ti, Ti);
        k2 = odefun(ti + h/2, Ti + h/2 * k1);
        k3 = odefun(ti + h/2, Ti + h/2 * k2);
        k4 = odefun(ti + h, Ti + h * k3);
        
        % Aktualizacja stanu
        T(i+1, :) = Ti' + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)';
    end
end