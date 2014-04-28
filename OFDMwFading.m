%% OFDM with a Fading Channel
%% Constant Assignment
clear all
close all
clc
m=2;    %m = log2(M) is the number bits per symbol or M = 2^m
M=4;    %number of symbols in M-ary modulation
L=16;   %size of IFFT/FFT
N=200;  %number of OFDM symbols

%% (a) Generate Input Bit Sequence
% Generate a random sequence bk of 1?s and 0?s with equal probability of 
%length m×L×N. Suggested values are: m = 2; L = 16 (64 and 256 are 
%optional); and N = 200 or more.
bits_bk =randi([0 1],[N*L*M 1]);    %Input Bit Sequence

%% (b) Convert Bit Sequence to Bit Group Sequence
bitsTo_M_Vectors_mk = reshape(bits_bk,length(bits_bk)/M,M);     %Bit Group
DecimalVector = bi2de(bitsTo_M_Vectors_mk);
Symbols_Xk = qammod(DecimalVector,L,0,'gray');  %Constellation Symbol Group
SymbolsIFFT = ifft(Symbols_Xk);    
sigmaValues = [0.0, 0.02, 0.08];    %Standard Deviation 
randQ=randn(length(Symbols_Xk),1)*i;    %Creates random variations of Q
randI=randn(length(Symbols_Xk),1);      %Creates random variations of I
Noise=(randQ+randI)*sigmaValues;        %Imposes Noise to the signal vector

%h is the fading simulation
h=sqrt(abs(Noise));

%SingleCarrier Modulation,fading Channel
SymbolsOut_rk_SingleCarrier(:,1)=Symbols_Xk+h(:,1)+h(:,1)*i;
SymbolsOut_rk_SingleCarrier(:,2)=Symbols_Xk+h(:,2)+h(:,2)*i;
SymbolsOut_rk_SingleCarrier(:,3)=Symbols_Xk+h(:,3)+h(:,3)*i;

%Multicarrier Modulation,fading channel
SymbolsOut_rk_woNoise(:,1)=SymbolsIFFT+h(:,1)+h(:,1)*i;
SymbolsOut_rk_woNoise(:,2)=SymbolsIFFT+h(:,2)+h(:,2)*i;
SymbolsOut_rk_woNoise(:,3)=SymbolsIFFT+h(:,3)+h(:,3)*i;

%MultiCarrier Modulation, AWGN channel
SymbolsOut_rk(:,1)=SymbolsIFFT+Noise(:,1);
SymbolsOut_rk(:,2)=SymbolsIFFT+Noise(:,2);
SymbolsOut_rk(:,3)=SymbolsIFFT+Noise(:,3);

%MultiCarrier Modulation, Raleigh fading
SymbolsInFFT_Rayleigh(:,1) = fft(SymbolsOut_rk_woNoise(:,1));
SymbolsInFFT_Rayleigh(:,2) = fft(SymbolsOut_rk_woNoise(:,2));
SymbolsInFFT_Rayleigh(:,3) = fft(SymbolsOut_rk_woNoise(:,3));

%Input Symbols with AWGN 
SymbolsInFFT_AWGM(:,1) = fft(SymbolsOut_rk(:,1));
SymbolsInFFT_AWGM(:,2) = fft(SymbolsOut_rk(:,2));
SymbolsInFFT_AWGM(:,3) = fft(SymbolsOut_rk(:,3));

%Demodulation

SymbolsIn_dk_Single(:,1) = qamdemod(SymbolsOut_rk_SingleCarrier(:,1),L,0,'gray');
SymbolsIn_dk_Single(:,2) = qamdemod(SymbolsOut_rk_SingleCarrier(:,2),L,0,'gray');
SymbolsIn_dk_Single(:,3) = qamdemod(SymbolsOut_rk_SingleCarrier(:,3),L,0,'gray');

%Rayleigh simulated fading
SymbolsIn_dk_Rayleigh(:,1) = qamdemod(SymbolsInFFT_Rayleigh(:,1),L,0,'gray');
SymbolsIn_dk_Rayleigh(:,2) = qamdemod(SymbolsInFFT_Rayleigh(:,2),L,0,'gray');
SymbolsIn_dk_Rayleigh(:,3) = qamdemod(SymbolsInFFT_Rayleigh(:,3),L,0,'gray');

%AWGM-MultiCarrier
SymbolsIn_dk(:,1) = qamdemod(SymbolsInFFT_AWGM(:,1),L,0,'gray');
SymbolsIn_dk(:,2) = qamdemod(SymbolsInFFT_AWGM(:,2),L,0,'gray');
SymbolsIn_dk(:,3) = qamdemod(SymbolsInFFT_AWGM(:,3),L,0,'gray');

%Graphs
splotfig1=scatterplot(SymbolsOut_rk_SingleCarrier(:,1),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig1);
title(strcat('Single Carrier: L = ', {' '},num2str(L),' sigma = 0.0'))
axis([-M M -M M])

splotfig2=scatterplot(SymbolsOut_rk_SingleCarrier(:,1),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig2);
title(strcat('Single Carrier: L = ', {' '},num2str(L),' sigma = 0.02'))
axis([-M M -M M])

splotfig3=scatterplot(SymbolsOut_rk_SingleCarrier(:,1),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig3);
title(strcat('Single Carrier: L = ', {' '},num2str(L),' sigma = 0.08'))
axis([-M M -M M])

splotfig4=scatterplot(SymbolsInFFT_Rayleigh(:,1),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig4);
title(strcat('Multiple Carrier: L = ', {' '},num2str(L),' sigma = 0.0'))
axis([-M M -M M])

splotfig5=scatterplot(SymbolsInFFT_Rayleigh(:,2),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig5);
title(strcat('Multiple Carrier: L = ', {' '},num2str(L),' sigma = 0.02'))
axis([-M M -M M])

splotfig6=scatterplot(SymbolsInFFT_Rayleigh(:,3),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig6);
title(strcat('Multiple Carrier: L = ', {' '},num2str(L),' sigma = 0.08'))
axis([-M M -M M])

splotfig7=scatterplot(SymbolsInFFT_AWGM(:,1),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig7);
title(strcat('Multiple Carrier AWGN: L = ', {' '},num2str(L),' sigma = 0.0'))
axis([-M M -M M])

splotfig8=scatterplot(SymbolsInFFT_AWGM(:,2),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig8);
title(strcat('Multiple Carrier AWGN: L = ', {' '},num2str(L),' sigma = 0.02'))
axis([-M M -M M])

splotfig9=scatterplot(SymbolsInFFT_AWGM(:,3),1,0,'g.');
hold on
scatterplot(Symbols_Xk,1,0,'k*',splotfig9);
title(strcat('Multiple Carrier AWGN: L =  ', {' '},num2str(L),', sigma = 0.08'))
axis([-M M -M M])

%Convert to binary

%Rayleigh-Single Carrier
DecimalVector_To_Binary_ck_Single(:,1:M) = de2bi(SymbolsIn_dk_Single(:,1),M);
DecimalVector_To_Binary_ck_Single(:,M+1:2*M) = de2bi(SymbolsIn_dk_Single(:,2),M);
DecimalVector_To_Binary_ck_Single(:,2*M+1:3*M) = de2bi(SymbolsIn_dk_Single(:,3),M);

%Rayleigh-MultiCarrier
DecimalVector_To_Binary_ck_Rayleigh(:,1:M) = de2bi(SymbolsIn_dk_Rayleigh(:,1),M);
DecimalVector_To_Binary_ck_Rayleigh(:,M+1:2*M) = de2bi(SymbolsIn_dk_Rayleigh(:,2),M);
DecimalVector_To_Binary_ck_Rayleigh(:,2*M+1:3*M) = de2bi(SymbolsIn_dk_Rayleigh(:,3),M);

%AWGM-MultiCarrier
DecimalVector_To_Binary_ck(:,1:M) = de2bi(SymbolsIn_dk(:,1),M);
DecimalVector_To_Binary_ck(:,M+1:2*M) = de2bi(SymbolsIn_dk(:,2),M);
DecimalVector_To_Binary_ck(:,2*M+1:3*M) = de2bi(SymbolsIn_dk(:,3),M);

%bits received

bitsIn_ck_Single(:,1)=reshape(DecimalVector_To_Binary_ck_Single(:,1:M),[N*L*M 1]);
bitsIn_ck_Single(:,2)=reshape(DecimalVector_To_Binary_ck_Single(:,M+1:2*M),[N*L*M 1]);
bitsIn_ck_Single(:,3)=reshape(DecimalVector_To_Binary_ck_Single(:,2*M+1:3*M),[N*L*M 1]);

bitsIn_ck_Rayleigh(:,1)=reshape(DecimalVector_To_Binary_ck_Rayleigh(:,1:M),[N*L*M 1]);
bitsIn_ck_Rayleigh(:,2)=reshape(DecimalVector_To_Binary_ck_Rayleigh(:,M+1:2*M),[N*L*M 1]);
bitsIn_ck_Rayleigh(:,3)=reshape(DecimalVector_To_Binary_ck_Rayleigh(:,2*M+1:3*M),[N*L*M 1]);

bitsIn_ck(:,1)=reshape(DecimalVector_To_Binary_ck(:,1:M),[N*L*M 1]);
bitsIn_ck(:,2)=reshape(DecimalVector_To_Binary_ck(:,M+1:2*M),[N*L*M 1]);
bitsIn_ck(:,3)=reshape(DecimalVector_To_Binary_ck(:,2*M+1:3*M),[N*L*M 1]);




fprintf('\nL = %d     N = %d\n',L,N)

fprintf('\n Sigma = 0.0\n')

[numErrors,ber]=biterr(bits_bk,bitsIn_ck_Single(:,1));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for single carrier fading.\n', ber,numErrors)

[numErrors,ber]=biterr(bits_bk,bitsIn_ck_Rayleigh(:,1));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for multiple carrier fading.\n', ber,numErrors)

[numErrors,ber]=biterr(bits_bk,bitsIn_ck(:,1));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for AWGN.\n', ber,numErrors)

fprintf('\n Sigma = 0.02\n')

[numErrors,ber]=biterr(bits_bk,bitsIn_ck_Single(:,2));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for single carrier fading.\n', ber,numErrors)

[numErrors,ber]=biterr(bits_bk,bitsIn_ck_Rayleigh(:,2));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for multiple carrier fading.\n', ber,numErrors)

[numErrors,ber]=biterr(bits_bk,bitsIn_ck(:,2));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for AWGN.\n', ber,numErrors)

fprintf('\n Sigma = 0.08\n')

[numErrors,ber]=biterr(bits_bk,bitsIn_ck_Single(:,3));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for single carrier fading.\n', ber,numErrors)

[numErrors,ber]=biterr(bits_bk,bitsIn_ck_Rayleigh(:,3));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for multiple carrier fading.\n', ber,numErrors)

[numErrors,ber]=biterr(bits_bk,bitsIn_ck(:,3));

fprintf('\nThe bit error rate = %5.2e, based on %d errors for AWGN.\n', ber,numErrors)

