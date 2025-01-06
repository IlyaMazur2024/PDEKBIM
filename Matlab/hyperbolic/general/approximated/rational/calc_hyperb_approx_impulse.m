
function [g11Ni,g12Ni,g21Ni,g22Ni] = calc_hyperb_approx_impulse(GN,t)

% funkcja na podstawie macierzy transmitancji operatorowych sekcyjnego modelu 
% aproksymacyjnego uk�adu hiperbolicznego 2x2 z�o�onego z N sekcji
% oblicza odpowiedzi impulsowe poszczeg�lnych tor�w transmitancji

g11Ni = impulse(GN(1,1),t);
g12Ni = impulse(GN(1,2),t);
g21Ni = impulse(GN(2,1),t);
g22Ni = impulse(GN(2,2),t);


