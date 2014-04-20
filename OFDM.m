clc;
clear all;
%M: number of symbols in M-ary modulation 
%m = log2(M) is the number bits per symbol or M = 2^m
%L: size of IFFT/FFT 
%N:  number of OFDM symbols 
%k:  index 
%bsubk   input bit sequence 
%msubk   sequence of groups of m bits  
%Xsubk   sequence of constellation symbols to be input to IFFT 
%ssubk   IFFT output sequence 
%rsubk received sequence to be input to FFT 
%Ysubk    FFT output values 
%dsubk   sequence of recovered groups of m bits 
%csubk   received bit sequence 
%sigma standard deviation of AWGN
m=4;    %bits per symbol
L=16;  %size of IFFT/FFT  0
L64=64;  %size of IFFT/FFT 1
L254=256;  %size of IFFT/FFT 2
N=10;  %Number of OFDM symboL254 Must be N>=100

%% Generate Random Sequence
%(a) Generate a random sequence bkof 1’s and 0’s with equal probability of 
%   length mxLxN. Suggested values are: m = 4; L = 16, 64, and 256; and 
%   N = 100 or more.  
for Nu = 1:N
    bsubk(:,:,Nu)=randi([0 1],m,L);  %is a random sequence of ones and 
                                            %zeros with equal probability
end

%% Generate Bit Group Sequence msubk
%(b) Let us consider QAM16, i.e., m = 4 or M = 16. A higher value is an 
%   option. Convert groups of b bits mk into a sequence of unsigned decimal 
%   values (e.g., 0, 1, 2, …, 14, 15). 
for Nu = 1:N
    msubk(:,Nu)=sum(bsubk(:,:,Nu));
end
bsubk(:,:,1)
msubk(:,1)'

%% Make constellation Symbol Sequence Xsubk
%(c) Use the constellation diagram below to map mk into constellation 
%   symbol sequence Xk. These are complex values (I-Q data).
%   This block is modified code from:
%   http://www.mathworks.com/help/comm/ug/constellation-visualization.html#bs2bgcd
% Create 8-QAM Gray encoded modulator
hMod = modem.qammod('M',16,'SymbolOrder','User-defined')
hMod.SymbolMapping = [6 4 12 14 7 5 13 15 3 1 9 11 2 0 8 10]
% Create a scatter plot
scatterPlot = commscope.ScatterPlot('SamplesPerSymbol',4,...
    'Constellation',hMod.Constellation);
% Show constellation
scatterPlot.PlotSettings.Constellation = 'on';
scatterPlot.PlotSettings.ConstellationStyle = '.';
% Add symbol labels
hold on;
k=log2(hMod.M);
for jj=1:hMod.M
        text(real(hMod.Constellation(jj))+0.15,...,
        imag(hMod.Constellation(jj)),...
        dec2base(hMod.SymbolMapping(jj),2,k));
end
hold off;


%% IFFT
%(d) Take a block of size L constellation symbols and apply the IFFT 
%   algorithm. Repeat this N times. These are complex valued too. 

%% White Noise
%(e) Add (I-Q) white Gaussian noise with zero mean; suggested standard 
%   deviation values are sigma=0 , 0.02, 0.08


%% FFT
%(f)  Input blocks of L received symbols from the channel to FFT. The 
%   FFT output values are the recovered constellation symbols. 

%% Demap
%(g) Implement demapping to convert the received constellation symbols 
%   into a bit stream. 

%% Calculate error ratio
%(h) Compare the received bit stream with the original and determine the 
%   error ratio (number of errors divided by the total number of bits). 

%% L=64
%(i)  Repeat (a)-(h) for L = 16, 64, and 256.

%% L=256