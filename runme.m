% Clean the previous runs
close all;
clear;
clc;

%%  Read Image
img = imread('./Graph_Source/lena512noisy.bmp');
% img = imread('./Graph_Source/Example1.png');
% show the read image
img_FFT = fft2(img);

[LoD, HiD, LoR, HiR] = wfilters('db1'); 
% expansion for avoiding edge disortion
%   'sym' or 'symh': Symmetric extension (half point): boundary value
%   symmetric replication for even-length filters

%% Decompose noisy image into subbands using
%   16-bancd dyadic (LoD - LPF, HiD - HPF)

[cA, cH, cV, cD] = dwt2(img,LoD,HiD,'mode','symh');  % cA=LL, cH=LH, cV=HL, cD=HH
figure;
subplot(2,2,1);imagesc(cA);
title("Level 1 Dec LL")
subplot(2,2,3);imagesc(cH);
title("Level 1 Dec LH")
subplot(2,2,2);imagesc(cV);
title("Level 1 Dec HL")
subplot(2,2,4);imagesc(cD);
title("Level 1 Dec HH")
colormap(gray)

[y1,y2,y3,y4] = dwt2(cA,LoD,HiD,'mode','symh');
[y5,y6,y7,y8] = dwt2(cH,LoD,HiD,'mode','symh');
[y9,y10,y11,y12] = dwt2(cV,LoD,HiD,'mode','symh');
[y13,y14,y15,y16] = dwt2(cD,LoD,HiD,'mode','symh');

% [y01, y02, y03, y04] = dwt2(cA,2,LoR,HiR);

%   22-band modified pyramid
%   initially decomposed into 16 equal-sized subbands
[z1,z2,z3,z4] = dwt2(cA,LoD,HiD,'mode','symh');
[z5,z6,z7,z8] = dwt2(cH,LoD,HiD,'mode','symh');
[z9,z10,z11,z12] = dwt2(cV,LoD,HiD,'mode','symh');
[z13,z14,z15,z16] = dwt2(cD,LoD,HiD,'mode','symh');
%   then add two additional decomposition to the lowest subband.
[z1_1,z1_2,z1_3,z1_4] = dwt2(z1,LoD,HiD,'mode','symh');
[z1_1_1, z1_1_2, z1_1_3, z1_1_4] = dwt2(z1_1,LoD,HiD,'mode','symh');
%   Use dwt filter with eliminate edge effects

%%  Remove high-frequency noise signal from frequency domain subbands
%   Dyadic case
%       highest-frequency = 0 (size = 131*131)
%       3 highest-frequency = 0 (sizes = 131*131)
%       6 highest-frequency = 0 (sizes = 131*131)
%   Modified pyramid case
%       3 highest-frequency = 0 (sizes = 131*131)
%       10 highest-frequency = 0    (sizes = 131*131)
%       15 highest-frequency = 0    (sizes = 131*131)
zero_arr = zeros(length(y16)); % will apply this zero_arr back




%%  Form reconstructed image in each case by Inverse dwt
%   Dyadic case
%       highest-frequency = 0
%   y16
cAR = idwt2(y1,y2,y3,y4,LoR,HiR);
cHR = idwt2(y5,y6,y7,y8,LoR,HiR);
cVR = idwt2(y9,y10,y11,y12,LoR,HiR);
cDR = idwt2(y13,y14,y15,zero_arr,LoR,HiR);
img_D1 = idwt2(cAR,cHR,cVR,cDR,LoR,HiR);
img_FFT_D1 = fft2(img_D1);
% double to uint8 for output format (0~255); or double (0~1)->img./256 
img_D1 = uint8(img_D1);

%       3 highest-frequency = 0
%   y14, y15, y16
cAR = idwt2(y1,y2,y3,y4,LoR,HiR);
cHR = idwt2(y5,y6,y7,y8,LoR,HiR);
cVR = idwt2(y9,y10,y11,y12,LoR,HiR);
cDR = idwt2(y13,zero_arr,zero_arr,zero_arr,LoR,HiR);
img_D2 = idwt2(cAR,cHR,cVR,cDR,LoR,HiR);
img_FFT_D2 = fft2(img_D2);
% double to uint8 for output format (0~255); or double (0~1)->img./256 
img_D2 = uint8(img_D2);

%       6 highest-frequency = 0
%   y8, y12, y13, y14, yh15, y16
cAR = idwt2(y1,y2,y3,y4,LoR,HiR);
cHR = idwt2(y5,y6,y7,zero_arr,LoR,HiR);
cVR = idwt2(y9,y10,y11,zero_arr,LoR,HiR);
cDR = idwt2(zero_arr,zero_arr,zero_arr,zero_arr,LoR,HiR);
img_D3 = idwt2(cAR,cHR,cVR,cDR,LoR,HiR);
img_FFT_D3 = fft2(img_D3);
% double to uint8 for output format (0~255); or double (0~1)->img./256 
img_D3 = uint8(img_D3);

%   Modified pyramid case
%       3 highest-frequency = 0
%   z14, z15, z16
z1_1 = idwt2(z1_1_1,z1_1_2,z1_1_3,z1_1_4,LoR,HiR);
z1 = idwt2(z1_1,z1_2,z1_3,z1_4,LoR,HiR);
cAR = idwt2(z1,z2,z3,z4,LoR,HiR);
cHR = idwt2(z5,z6,z7,z8,LoR,HiR);
cVR = idwt2(z9,z10,z11,z12,LoR,HiR);
cDR = idwt2(z13,zero_arr,zero_arr,zero_arr,LoR,HiR);
img_MP1 = idwt2(cAR,cHR,cVR,cDR,LoR,HiR);
img_FFT_MP1 = fft2(img_MP1);
% double to uint8 for output format (0~255); or double (0~1)->img./256 
img_MP1 = uint8(img_MP1);

%       10 highest-frequency = 0
%   z6, z7, z8, z10, z11, z12, z13, z14, z15, z16
z1_1 = idwt2(z1_1_1,z1_1_2,z1_1_3,z1_1_4,LoR,HiR);
z1 = idwt2(z1_1,z1_2,z1_3,z1_4,LoR,HiR);
cAR = idwt2(z1,z2,z3,z4,LoR,HiR);
cHR = idwt2(z5,zero_arr,zero_arr,zero_arr,LoR,HiR);
cVR = idwt2(z9,zero_arr,zero_arr,zero_arr,LoR,HiR);
cDR = idwt2(zero_arr,zero_arr,zero_arr,zero_arr,LoR,HiR);
img_MP2 = idwt2(cAR,cHR,cVR,cDR,LoR,HiR);
img_FFT_MP2 = fft2(img_MP2);
% double to uint8 for output format (0~255); or double (0~1)->img./256 
img_MP2 = uint8(img_MP2);

%       15 highest-frequency = 0
%   z2:Z16
z1_1 = idwt2(z1_1_1,z1_1_2,z1_1_3,z1_1_4,LoR,HiR);
z1 = idwt2(z1_1,z1_2,z1_3,z1_4,LoR,HiR);
cAR = idwt2(z1,zero_arr,zero_arr,zero_arr,LoR,HiR);
cHR = idwt2(zero_arr,zero_arr,zero_arr,zero_arr,LoR,HiR);
cVR = idwt2(zero_arr,zero_arr,zero_arr,zero_arr,LoR,HiR);
cDR = idwt2(zero_arr,zero_arr,zero_arr,zero_arr,LoR,HiR);
img_MP3 = idwt2(cAR,cHR,cVR,cDR,LoR,HiR);
img_FFT_MP3 = fft2(img_MP3);
% double to uint8 for output format (0~255); or double (0~1)->img./256 
img_MP3 = uint8(img_MP3);

%%  Figure Plotting
%   Original noisy plot
figure;
imagesc(img)
title("Noisy Figure")
colormap(gray)

figure
plot(10*log10(abs(img_FFT)))
xlim([1 length(img_FFT)])
title("Noisy Figure FFT (Mag)")

%   Dyadic case
%       highest-frequency = 0
figure
imagesc(img_D1)
title("Dyadic, One HF = 0")
colormap(gray)

imwrite(img_D1,'img_Dyadic_one_HF.bmp')

figure
plot(10*log10(abs(img_FFT_D1)))
xlim([1 length(img_FFT_D1)])
title("Dyadic, one HF FFT (Mag)")

%       3 highest-frequency = 0
figure
imagesc(img_D2)
title("Dyadic, Three HF = 0")
colormap(gray)

imwrite(img_D2,'img_Dyadic_three_HF.bmp')

figure
plot(10*log10(abs(img_FFT_D2)))
xlim([1 length(img_FFT_D2)])
title("Dyadic, three HF FFT (Mag)")

%       6 highest-frequency = 0
figure
imagesc(img_D3)
title("Dyadic, Six HF = 0")
colormap(gray)

imwrite(img_D3,'img_Dyadic_six_HF.bmp')

figure
plot(10*log10(abs(img_FFT_D3)))
xlim([1 length(img_FFT_D3)])
title("Dyadic, six HF FFT (Mag)")

%   Modified pyramid case
%       3 highest-frequency = 0
figure
imagesc(img_MP1)
title("Modified Pyramid, Three HF = 0")
colormap(gray)

imwrite(img_MP1,'img_MP_three_HF.bmp')

figure
plot(10*log10(abs(img_FFT_MP1)))
xlim([1 length(img_FFT_MP1)])
title("MP, three HF FFT (Mag)")

%       10 highest-frequency = 0
figure
imagesc(img_MP2)
title("Modified Pyramid, Ten HF = 0")
colormap(gray)
imwrite(img_MP2,'img_MP_ten_HF.bmp')

figure
plot(10*log10(abs(img_FFT_MP2)))
xlim([1 length(img_FFT_MP2)])
title("MP, ten HF FFT (Mag)")

%       15 highest-frequency = 0
figure
imagesc(img_MP3)
title("Modified Pyramid, 15 HF = 0")
colormap(gray)

imwrite(img_MP3,'img_MP_15_HF.bmp')

figure
plot(10*log10(abs(img_FFT_MP3)))
xlim([1 length(img_FFT_MP3)])
title("MP, 15 HF FFT (Mag)")
