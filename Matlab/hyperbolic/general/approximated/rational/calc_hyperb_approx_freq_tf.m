
function [G11Nw,G12Nw,G21Nw,G22Nw] = calc_hyper_approx_freq_tf(GN,omega)

% funkcja na podstawie macierzy transmitancji operatorowych sekcyjnego modelu 
% aproksymacyjnego uk³adu hiperbolicznego 2x2 z³o¿onego z N sekcji
% oblicza odpowiedzi czêstotliwoœciowe poszczególnych torów transmitancji

F = freqresp(GN(1,1),omega);
G11Nw = reshape(F,[length(omega) 1]);
F = freqresp(GN(1,2),omega);
G12Nw = reshape(F,[length(omega) 1]);
F = freqresp(GN(2,1),omega);
G21Nw = reshape(F,[length(omega) 1]);
F = freqresp(GN(2,2),omega);
G22Nw = reshape(F,[length(omega) 1]);




