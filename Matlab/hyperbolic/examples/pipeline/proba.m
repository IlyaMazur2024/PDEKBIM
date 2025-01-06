
load odpowiedzi_uszk

Fs=10;
t=[0:0.1:60];
L=length(t);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
omegas = Fs*2*pi;


a=zeros(5,NFFT/2+1) % frequency responses matrix


% dla ró¿nych impedancji wzgl. wycieku (0.01 0.1 0.5 1)
figure(1)
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (rad)')
ylabel('|Y(f)|')
col = ['k' 'r' 'c' 'g' 'b'];

for i=1:5
    y=hG1(i,:);
    Y = fft(y,NFFT)/L;
    omega = omegas/2*linspace(0,1,NFFT/2+1);
    a(i,:) = 2*abs(Y(1:NFFT/2+1));
    plot(omega,a(i,:),col(i)) 
    hold on
end    


% dla ró¿nych lokalizacji wzgl. wycieku (0.25 0.5 0.75 0.99)
figure(2)
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (rad)')
ylabel('|Y(f)|')
col = ['k' 'r' 'c' 'g' 'b'];

for i=1:5
    y=hG2(i,:);
    Y = fft(y,NFFT)/L;
    omega = omegas/2*linspace(0,1,NFFT/2+1);
    a(i,:) = 2*abs(Y(1:NFFT/2+1));
    plot(omega,a(i,:),col(i)) 
    hold on
end    


