clear all;
close all;
clc;

N=512;
O=phantom('Modified Shepp-Logan', N);
figure(1),imshow(O);title('Original Image');


%cutoff frequency of the filters
cutoff= [0.5, 0.6, 0.7, 0.8, 0.9]; %cutoff frequencies in percentage (%) of Nyquist Frequency:
%cutoff in freq= [0.25, 0.3, 0.35, 0.4, 0.45] (mm-1);

%order of the filters
order=floor(linspace(0, 20, 5));

%Assumed Pixel size
pixel_size = 1     % [mm]
 
fs = 1/pixel_size              % sampling frequency [1/mm] = 1 mm-1
delta_f = fs/N          % frequency resolution

freq_vector2=[-fs/2:delta_f:fs/2 - delta_f];

%% PROJECTIONS
% proj= number of projections

proj=512;
angle= (180/proj);
sino_O= zeros(N,proj);
for i=0:proj-1
   
    %rotation
    O_rot=imrotate(O,i*angle,'crop');
    %projection
    i=i+1;
    for j=1:N
        sino_O(j,i)= sum(O_rot(:,j));
    end
end
% figure
% imagesc(sino_O);title('Non-filtered Sinogram');

%% 0 Ramp filtering and reconstruction

% Creation of the ramp filter
freqs=linspace(-1, 1, N).';
filterline = abs( freqs );
figure
plot(freq_vector2,filterline), xlim([-0.5 0.5]), title('Ramp Filter'), xlabel('Frequency (mm-1)'), ylabel('Amplitude');
Filter = repmat(filterline, [1 proj]);

% FT domain filtering
ft_sino_O = fftshift(fft(sino_O,[],1),1);
figure,plot(freq_vector2,abs(ft_sino_O(:,1))), hold on, plot(freq_vector2, filterline), legend('FFT of projection', 'Ramp Filter');
filteredProj_0 = ft_sino_O .* Filter;
filteredProj_0 = ifftshift(filteredProj_0,1);
ift_sino_O = real(ifft(filteredProj_0,[],1));
figure
subplot(2,1,1),imagesc(abs(ift_sino_O)), title('Filtered Sinogram');
subplot(2,1,2),plot(abs(ift_sino_O(:,1))), title('Ramp Filtered Projection');

% Sinogram in the frequency domain
% imagesc(abs(ft_sino)); title('FFT of the sinogram');
% Filtered sinogram
% imagesc(ift_sino);
% title('filtered sinogram');


% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_O(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Reconstructed image, ramp filter');
R0= Rsum;

%% 1 Ramp * Hamming filtering
R1_cutoff = [];
filteredProj_1 = [];
for i = 1:length(cutoff)
    Filter1=[];
    filterline1=[];
    % Creation of the window*ramp filter
    freqs=linspace(-1, 1, N).';
    m= fir1(order(3)-1, cutoff(i), 'low', hamming(order(3)));
    [whamming,w]=freqz(m,1,N/2);
    % figure
    % plot(abs(whamming))
    %Mirroring the Hamming window and centered it at the middle as the
    %fftshifted sinogram is done to put the zero frequencies in the middle
    whamming=abs((whamming)');
    windowHamming=zeros(1,N);
    windowHamming(N/2+1:N)=whamming(1:N/2);
    count=1;
    i=0;
    for i=N/2:-1:1
        windowHamming(i)=whamming(count);
        count=count+1;
    end
   
    filterline1 = abs( freqs ).* windowHamming';
    
    Filter1 = repmat(filterline1, [1 proj]);
    % FT domain filtering
    ft_sino_O = fftshift(fft(sino_O,[],1),1);
    filteredProj = ft_sino_O .* Filter1;
    filteredProj = ifftshift(filteredProj,1);
    filteredProj_1 = [filteredProj_1, filteredProj];
    ift_sino_O = real(ifft(filteredProj,[],1));
    
 
    % Retroprojection
    R=zeros(N,N);
    Rsum=zeros(N,N);
    for i=1:proj
        for j= 1:N
            R(:,j)= ift_sino_O(j,i);
        end
        R= imrotate(R,-i*angle,'crop');
        Rsum=Rsum+R;
    end
    
    % figure
    % imshow(Rsum,[]);title('Recontructed image, ramp * Hamming window filter');
    R1= Rsum;
    R1_cutoff = [R1_cutoff,R1];
end
nor_freq=linspace(0,1,N);
figure 
plot(nor_freq, windowHamming), xlabel('Normalized Frequency'), ylabel('Amplitude'), title('Hamming Window Filter for cutoff=0.45');
    
R1_1=R1_cutoff(:,1:N);
R1_2=R1_cutoff(:,N+1:N+N);
R1_3=R1_cutoff(:,N+N+1:3*N);
R1_4=R1_cutoff(:,3*N+1:4*N);
R1_5=R1_cutoff(:,4*N+1:5*N);


mse_1 = [immse(R1_1,O), immse(R1_2,O) ,immse(R1_3,O), immse(R1_4,O), immse(R1_5,O)];
% mse_1_mean = mean(mse_1);

figure
subplot(2,3,1), imshow(R1_1,[]);title('Ramp * Hamming , cutoff = 0.25');
subplot(2,3,2), imshow(R1_2,[]);title('Ramp * Hamming ,cutoff = 0.30');
subplot(2,3,3),imshow(R1_3,[]);title('Ramp * Hamming ,cutoff = 0.35');
subplot(2,3,4),imshow(R1_4,[]);title('Ramp * Hamming ,cutoff = 0.40');
subplot(2,3,5),imshow(R1_5,[]);title('Ramp * Hamming ,cutoff = 0.45');

R1_1nor=abs(R1_1);
max_11= max(max(R1_1nor));
R1_1nor=R1_1./max_11;
% figure, imshow(R1_nor);

R1_2nor=abs(R1_2);
max_12= max(max(R1_2nor));
R1_2nor=R1_2./max_12;
% figure, imshow(R2_nor);

R1_3nor=abs(R1_3);
max_13= max(max(R1_3nor));
R1_3nor=R1_3./max_13;
%figure, imshow(R3_nor);


R1_4nor=abs(R1_4);
max_14= max(max(R1_4nor));
R1_4nor=R1_4./max_14;
%figure, imshow(R4_nor);

R1_5nor=abs(R1_5);
max_15= max(max(R1_5nor));
R1_5nor=R1_5./max_15;
%figure, imshow(R5_nor);


%% Metrics for choosing the cutoff freq
% Peak SNR
peaksnr=zeros(5,1);
peaksnr(1) = psnr(R1_1nor,O);
peaksnr(2) = psnr(R1_2nor,O);
peaksnr(3) = psnr(R1_3nor,O);
peaksnr(4) = psnr(R1_4nor,O);
peaksnr(5) = psnr(R1_5nor,O);
% figure;
% bar(peaksnr); title('peak SNR for different cutoff freq');

% Mean square ratio
mse=zeros(5,1);
mse(1)=immse(R1_1nor,O);
mse(2)=immse(R1_2nor,O);
mse(3)=immse(R1_3nor,O);
mse(4)=immse(R1_4nor,O);
mse(5)=immse(R1_5nor,O);
% figure;
% bar(mse); title('MSE for different cutoff freq');

% Structural similarity
ssimval=zeros(5,1);
ssimval(1)= ssim(R1_1nor,O);
ssimval(2)= ssim(R1_2nor,O);
ssimval(3)= ssim(R1_3nor,O);
ssimval(4)= ssim(R1_4nor,O);
ssimval(5)= ssim(R1_5nor,O);
% figure;
% bar(ssimval); title('SSIM for different cutoff freq');

%Retrive max/min
min_p=min(peaksnr);
index_min_peaksnr=find(peaksnr==min(peaksnr));
index_min_mse=find(mse==min(mse));
index_max_ssim=find(ssimval==min(ssimval))


%% 2 Ramp * Gaussian filtering
R2_cutoff = [];
for i = 1:length(cutoff)
Filter2=[];
filterline2=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
M= fir1(order(3)-1, cutoff(i), 'low', gausswin(order(3)));
[wgaussian,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wgaussian=abs((wgaussian)');
windowGaussian=zeros(1,N);
windowGaussian(N/2+1:N)=wgaussian(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowGaussian(i)=wgaussian(count);
    count=count+1;
end

filterline2 = abs( freqs ).* windowGaussian';


Filter2 = repmat(filterline2, [1 proj]);

% FT domain filtering
ft_sino_O= fftshift(fft(sino_O,[],1),1);
filteredProj = ft_sino_O .* Filter2;
filteredProj = ifftshift(filteredProj,1);
ift_sino_O = real(ifft(filteredProj,[],1));



% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_O(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, ramp * Gaussian window filter');
R2=Rsum;
R2_cutoff = [R2_cutoff,R2];
end
R2_1=R2_cutoff(:,1:N);
R2_2=R2_cutoff(:,N+1:N+N);
R2_3=R2_cutoff(:,N+N+1:N*3);
R2_4=R2_cutoff(:,3*N+1:4*N);
R2_5=R2_cutoff(:,4*N+1:5*N);

mse_2 = [immse(R2_1,O), immse(R2_2,O) ,immse(R2_3,O),immse(R2_4,O), immse(R2_5,O)];
mse_2_mean = mean(mse_2);
% figure
% subplot(3,5,1), imshow(R2_1,[]);title('Ramp * Gaussian , cutoff = 0.25');
% subplot(3,5,2), imshow(R2_2,[]);title('Ramp * Gaussian ,cutoff = 0.30');
% subplot(3,5,3),imshow(R2_3,[]);title('Ramp * Gaussian ,cutoff = 0.35');
% subplot(3,5,4),imshow(R2_4,[]);title('Ramp * Gaussian ,cutoff = 0.40');
% subplot(3,5,5),imshow(R2_5,[]);title('Ramp * Gaussian ,cutoff = 0.45');
%% 3 Ramp * Hanning filtering
R3_cutoff=[];
for i = 1:length(cutoff)
Filter3=[];
filterline3=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', hann(order(3)));
[whann,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
whann=abs((whann)');
windowHanning=zeros(1,N);
windowHanning(N/2+1:N)=whamming(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowHanning(i)=whann(count);
    count=count+1;
end
%plot(windowHanning);title('windowhanning');
filterline3 = abs( freqs ).* windowHanning';


Filter3 = repmat(filterline3, [1 proj]);

% FT domain filtering
ft_sino_O = fftshift(fft(sino_O,[],1),1);
filteredProj = ft_sino_O.* Filter3;
filteredProj = ifftshift(filteredProj,1);
ift_sino_O = real(ifft(filteredProj,[],1));



% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_O(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
    
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Hann window filter');
R3=Rsum;
R3_cutoff = [R3_cutoff,R3];
end
R3_1=R3_cutoff(:,1:N);
R3_2=R3_cutoff(:,N+1:N+N);
R3_3=R3_cutoff(:,N+N+1:3*N);
R3_4=R3_cutoff(:,3*N+1:4*N);
R3_5=R3_cutoff(:,4*N+1:5*N);


mse_3 = [immse(R3_1,O), immse(R3_2,O) ,immse(R3_3,O), immse(R3_4,O), immse(R3_5,O)];
mse_3_mean = mean(mse_3);

% subplot(3,5,6), imshow(R3_1,[]);title('Ramp * Hanning , cutoff = 0.25');
% subplot(3,5,7), imshow(R3_2,[]);title('Ramp * Hanning ,cutoff = 0.3');
% subplot(3,5,8), imshow(R3_3,[]);title('Ramp * Hanning ,cutoff = 0.35');
% subplot(2,5,9),imshow(R3_4,[]);title('Ramp * Hanning ,cutoff = 0.40');
% subplot(2,5,10),imshow(R3_5,[]);title('Ramp * Hanning ,cutoff = 0.45');
%% 4 Ramp * Parzen filtering
R4_cutoff=[];
for i = 1:length(cutoff)
Filter4=[];
filterline4=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', parzenwin(order(3)));
[wparzen,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wparzen=abs((wparzen)');
windowParzen=zeros(1,N);
windowParzen(N/2+1:N)=wparzen(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowParzen(i)=wparzen(count)
    count=count+1;
end
filterline4 = abs( freqs ).* windowParzen';

Filter4 = repmat(filterline4, [1 proj]);

% FT domain filtering
ft_sino_O = fftshift(fft(sino_O,[],1),1);
filteredProj = ft_sino_O .* Filter4;
filteredProj = ifftshift(filteredProj,1);
ift_sino_O = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_O(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Parzen window filter');
R4=Rsum;
R4_cutoff = [R4_cutoff,R4];
end
R4_1=R4_cutoff(:,1:N);
R4_2=R4_cutoff(:,N+1:N+N);
R4_3=R4_cutoff(:,N+N+1:3*N);
R4_4=R4_cutoff(:,3*N+1:4*N);
R4_5=R4_cutoff(:,4*N+1:5*N);

mse_4 = [immse(R4_1,O), immse(R4_2,O) ,immse(R4_3,O), immse(R4_4,O), immse(R4_5,O)];
mse_4_mean = mean(mse_4);

% subplot(3,5,11), imshow(R4_1,[]);title('Ramp * Parzen , cutoff = 0.25');
% subplot(3,5,12), imshow(R4_2,[]);title('Ramp * Parzen ,cutoff = 0.3');
% subplot(3,5,13), imshow(R4_3,[]);title('Ramp * Parzen ,cutoff = 0.35');
% subplot(2,5,14),imshow(R4_4,[]);title('Ramp * Parzen ,cutoff = 0.40');
% subplot(2,5,15),imshow(R4_5,[]);title('Ramp * Parzen ,cutoff = 0.45');


%% 5 Ramp * Bartlett filtering
R5_cutoff=[];
for i = 1:length(cutoff)
Filter5=[];
filterline5=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', bartlett(order(3)));
[wbartlett,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wbartlett=abs((wbartlett)');
windowBartlett=zeros(1,N);
windowBartlett(N/2+1:N)=wbartlett(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowBartlett(i)=wbartlett(count)
    count=count+1;
end
% plot(windowBartlett);title('windowbartlett');
filterline5 = abs( freqs ).* windowBartlett';


Filter5 = repmat(filterline5, [1 proj]);

% FT domain filtering
ft_sino_O = fftshift(fft(sino_O,[],1),1);
filteredProj = ft_sino_O .* Filter5;
filteredProj = ifftshift(filteredProj,1);
ift_sino_O = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_O(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Bartlett window filter');
R5=Rsum;
R5_cutoff = [R5_cutoff,R5];
end
R5_1=R5_cutoff(:,1:N);
R5_2=R5_cutoff(:,N+1:N+N);
R5_3=R5_cutoff(:,N+N+1:3*N);
R5_4=R5_cutoff(:,3*N+1:4*N);
R5_5=R5_cutoff(:,4*N+1:5*N);

mse_5 = [immse(R5_1,O), immse(R5_2,O) ,immse(R5_3,O), immse(R5_4,O), immse(R5_5,O)];
mse_5_mean = mean(mse_5);
% figure
% subplot(3,5,1), imshow(R5_1,[]);title('Ramp * Barlett , cutoff = 0.25');
% subplot(3,5,2), imshow(R5_2,[]);title('Ramp * Barlett ,cutoff = 0.3');
% subplot(3,5,3), imshow(R5_3,[]);title('Ramp * Barlett ,cutoff = 0.35');
% subplot(2,5,4),imshow(R5_4,[]);title('Ramp * Bartlett ,cutoff = 0.40');
% subplot(2,5,5),imshow(R5_5,[]);title('Ramp * Bartlett ,cutoff = 0.45');

%% 6 Ramp * Butterworth filtering
R6_cutoff=[];
for i = 1:length(cutoff)
Filter6=[];
filterline6=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
[m,l]= butter(order(2)-2, cutoff(i), 'low');
[h, w]= freqz(m,l,N/2);
windowButterworth = abs(h)';
windowButterworth = [flip(windowButterworth) windowButterworth];
%plot(windowButterworth);title('windowbutterworth');


filterline6 = abs( freqs ).* (windowButterworth)';


Filter6 = repmat(filterline6, [1 proj]);

% FT domain filtering
ft_sino_O = fftshift(fft(sino_O,[],1),1);
filteredProj = ft_sino_O .* Filter6;
filteredProj = ifftshift(filteredProj,1);
ift_sino_O = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_O(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
%      
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Butterworth window filter');
R6=Rsum;
R6_cutoff = [R6_cutoff,R6];
end
R6_1=R6_cutoff(:,1:N);
R6_2=R6_cutoff(:,N+1:N+N);
R6_3=R6_cutoff(:,N+N+1:3*N);
R6_4=R6_cutoff(:,3*N+1:4*N);
R6_5=R6_cutoff(:,4*N+1:5*N);


mse_6 = [immse(R6_1,O), immse(R6_2,O) ,immse(R6_3,O),immse(R6_4,O), immse(R6_5,O)];
mse_6_mean = mean(mse_6);

% figure
% subplot(2,3,1), imshow(R6_1,[]);title('Ramp * Butter , cutoff = 0.25');
% subplot(2,3,2), imshow(R6_2,[]);title('Ramp * Butter ,cutoff = 0.30');
% subplot(2,3,3), imshow(R6_3,[]);title('Ramp * Butter ,cutoff = 0.35');
% subplot(2,3,4),imshow(R6_4,[]);title('Ramp * Butter ,cutoff = 0.40');
% subplot(2,3,5),imshow(R6_5,[]);title('Ramp * Butter ,cutoff = 0.45');

%%
figure
plot(freq_vector2,filterline), xlabel('Frequency (1/mm)'), ylabel('Amplitude');
hold on
plot(freq_vector2,filterline1);
hold on
plot(freq_vector2,filterline2);
hold on
plot(freq_vector2,filterline3);
hold on
plot(freq_vector2, filterline4);
hold on
plot(freq_vector2,filterline5);
hold on
plot(freq_vector2, filterline6);
legend('Ramp Filter','Ramp*Hamming Filter', 'Ramp*Hanning Filter', 'Ramp*Gaussian Filter', 'Ramp*Parzen Filter', 'Ramp*Bartlett Filter', 'Ramp*ButterWorth');

figure
subplot(2,1,1),imshow(R0,[]);title('Recontructed image, Ramp window filter');
subplot(2,1,2), imagesc(sino_O);title('Non-filtered Sinogram');

figure
subplot(4,2,1), imshow(R1_3,[]);title('Recontructed image, Ramp * Hamming window filter');
subplot(4,2,2), imshow(R2_3, []);title('Recontructed image, Ramp * Hanning window filter');
subplot(4,2,3),imshow(R3_3,[]);title('Recontructed image, Ramp * Gaussian window filter');
subplot(4,2,4), imshow(R4_3,[]);title('Recontructed image, Ramp * Parzen window filter');
subplot(4,2,5),imshow(R5_3,[]);title('Recontructed image, Ramp * Bartlett window filter');
subplot(4,2,6), imshow(R6_3,[]);title('Recontructed image, Ramp * Butterworth window filter');
subplot(4,2,7), imshow(R0,[]);title('Recontructed image, Ramp window filter');

%Normalizing the reconstructed images for the gray scale range [0-1]
R1_nor=abs(R1_3);
max_1= max(max(R1_nor));
R1_nor=R1_3./max_1;
% figure, imshow(R1_nor);

R2_nor=abs(R2_3);
max_2= max(max(R2));
R2_nor=R2./max_2;
% figure, imshow(R2_nor);

R3_nor=abs(R3_3);
max_3= max(max(R3_nor));
R3_nor=R3./max_3;
%figure, imshow(R3_nor);


R4_nor=abs(R4_3);
max_4= max(max(R4_nor));
R4_nor=R4./max_4;
%figure, imshow(R4_nor);

R5_nor=abs(R5_3);
max_5= max(max(R5));
R5_nor=R5./max_5;
%figure, imshow(R5_nor);

R6_nor=abs(R6_3);
max_6= max(max(R6));
R6_nor=R6./max_6;
figure, imshow(R6_nor);

R0_nor=abs(R0);
max_0= max(max(R0));
R0_nor=R0./max_0;
%% Metrics for Original Image
% %teste for optimal Results: comparing original with original
% peaksnr_t=psnr(O,O);
% figure;
% bar(peaksnr_t), title('peasksnr')
% mse_t=immse(O,O);
% figure;
% bar(mse_t), title('mse');
% ssim_t=ssim(O,O);
% figure;
% bar(ssim_t), title('ssim');
% Peak SNR
peaksnr=zeros(7);
peaksnr(1) = psnr(R0_nor,O);
peaksnr(2) = psnr(R1_nor,O);
peaksnr(3) = psnr(R2_nor,O);
peaksnr(4) = psnr(R3_nor,O);
peaksnr(5) = psnr(R4_nor,O);
peaksnr(6) = psnr(R5_nor,O);
peaksnr(7) = psnr(R6_nor,O);
figure;
bar(peaksnr); title('peak SNR for Reconstructed Image from Original');

% Mean square ratio
mse=zeros(7);
mse(1)=immse(R0_nor,O);
mse(2)=immse(R1_nor,O);
mse(3)=immse(R2_nor,O);
mse(4)=immse(R3_nor,O);
mse(5)=immse(R4_nor,O);
mse(6)=immse(R5_nor,O);
mse(7)=immse(R6_nor,O);
figure;
bar(mse); title('MSE for  for Reconstructed Image from Original');

% Structural similarity
ssimval=zeros(7);
ssimval(1)= ssim(R0_nor,O);
ssimval(2)= ssim(R1_nor,O);
ssimval(3)= ssim(R2_nor,O);
ssimval(4)= ssim(R3_nor,O);
ssimval(5)= ssim(R4_nor,O);
ssimval(6)= ssim(R5_nor,O);
ssimval(7)= ssim(R6_nor,O);
figure;
bar(ssimval); title('SSIM for  for Reconstructed Image from Original');

%Comparing Original, Ramp and Rampp* Hamming
figure
subplot(2,2,1), imshow(O),title ('Original');
subplot(2,2,2), imshow(R0, []), title('Reconstructed using Ramp Filter');
subplot(2,2,3), imshow(R1_3, []), title('Reconstructed using Ramp*Hamming Filter');

%% Gaussian Noise Imposed

% adding noise
G = imnoise(O,'gaussian');
figure
imshow(G);title('Added Noise');


%% PROJECTIONS
% proj= number of projections
proj=512;
angle= (180/proj);
sino_G= zeros(N,proj);
for i=0:proj-1
   
    %rotation
    G_rot=imrotate(G,i*angle,'crop');
    %projection
    i=i+1;
    for j=1:N
        sino_G(j,i)= sum(G_rot(:,j));
    end
end
% figure
% imagesc(sino_G);title('Non-filtered Sinogram');

%% 0 Ramp filtering and reconstruction

% Creation of the ramp filter
freqs=linspace(-1, 1, N).';
filterline = abs( freqs );
Filter = repmat(filterline, [1 proj]);

% FT domain filtering
ft_sino_G = fftshift(fft(sino_G,[],1),1);
filteredProj_G_0 = ft_sino_G .* Filter;
filteredProj_G_0 = ifftshift(filteredProj_G_0,1);
ift_sino_G = real(ifft(filteredProj_G_0,[],1));

%Noise spectrum 
G_ft = abs(filteredProj_G_0) - abs(filteredProj_0);

% %Showing the filtering results
% pixel = 1:N;
% figure, plot(freqs,abs(ft_sino_G(:,100))),hold on, plot(freqs,abs(ft_sino_O(:,100))),title('FT of  projections Origina and Gaussian for projection Nb 100');
% %figure,plot(freqs,abs(filteredProj_0(:,100))),hold on, plot(freqs,abs(filteredProj_G_0(:,100))),title('FT of filtered projections of Original and Gaussian for projection Nb 100');
% 
figure,plot(freqs,abs(filteredProj_0(:,100)),'g'),hold on, plot(freqs,abs(filteredProj_G_0(:,100)),'r'),hold on,plot(freqs,abs(G_ft(:,100)),'b'), xlabel('Frequency'), ylabel('Amplitude'), title('Noise Spectrum');
legend('Original Filtered Spectrum Projection', 'Gaussian Filter Spectrum Projection', 'Gaussian Noise Spectrum');
% figure, imshow(ift_sino_O),title('Filtered Original sinogram');
% figure, imshow(ift_sino_G),title('Filtered Gaussian sinogram');

% Calculate SNR 
Power_Sino_G_0 = trapz(freqs,abs(filteredProj_G_0)); %Power of the sinogram
Power_Noise_G_0 = trapz(freqs,abs(G_ft));%Power of the Gaussian Noise
SNR_G_0 = sum(Power_Noise_G_0)/sum(Power_Sino_G_0);


% Sinogram in the frequency domain
% imagesc(abs(ft_sino_G)); title('FFT of the sinogram');
% Filtered sinogram
% imagesc(ift_sino_G);
% title('filtered sinogram');

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_G(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Reconstructed image, ramp filter');
R0_G= Rsum;

%% 1 Ramp * Hamming filtering
R1_G_cutoff=[];
SNR_G_1 = [];
for i =1:length(cutoff)
    Filter1=[];
    filterline1=[];
    % Creation of the window*ramp filter
    freqs=linspace(-1, 1, N).';
    m= fir1(order(3)-1, cutoff(i), 'low', hamming(order(3)));
    [whamming,w]=freqz(m,1,N/2);
    % figure
    % plot(abs(whamming))
    whamming=abs((whamming)');
    windowHamming=zeros(1,N);
    windowHamming(N/2+1:N)=whamming(1:N/2);
    count=1;
    i=0;
    for i=N/2:-1:1
        windowHamming(i)=whamming(count);
        count=count+1;
    end
    filterline1 = abs( freqs ).* windowHamming';
    
    
    Filter1 = repmat(filterline1, [1 proj]);
    
    % FT domain filtering
    ft_sino_G = fftshift(fft(sino_G,[],1),1);
    filteredProj = ft_sino_G .* Filter1;
    filteredProj = ifftshift(filteredProj,1);
    ift_sino_G = real(ifft(filteredProj,[],1));
    
    
    %Noise spectrum 
    G_ft = abs(filteredProj) - abs(filteredProj_1(:,(i-1)*proj+1:i*proj));
    % Calculate SNR 
    figure, plot(freqs,abs(G_ft),'b'),title('Noise Spectrum for Filtered Projections');
    hold on, plot(freqs,abs(filteredProj),'r'), xlabel('Frequency'), ylabel('Amplitude');

    Power_Sino_G_1 = trapz(freqs,abs(filteredProj)); %Power of the sinogram
    Power_Noise_G_1 = trapz(freqs,abs(G_ft));%Power of the Gaussian Noise
    SNR_G_1 = [SNR_G_1, sum(Power_Noise_G_1)/sum(Power_Sino_G_1)];
 
    
    % Retroprojection
    R=zeros(N,N);
    Rsum=zeros(N,N);
    for i=1:proj
        for j= 1:N
            R(:,j)= ift_sino_G(j,i);
        end
        R= imrotate(R,-i*angle,'crop');
        Rsum=Rsum+R;
    end
    
    % figure
    % imshow(Rsum,[]);title('Recontructed image, ramp * Hamming window filter');
    R1_G= Rsum;
    R1_G_cutoff = [R1_G_cutoff,R1_G];
end
R1_G_1=R1_G_cutoff(:,1:N);
R1_G_2=R1_G_cutoff(:,N+1:2*N);
R1_G_3=R1_G_cutoff(:,2*N+1:3*N);
R1_G_4=R1_G_cutoff(:,3*N+1:4*N);
R1_G_5=R1_G_cutoff(:,4*N+1:5*N);

mse_G = immse(G,O);
mse_0 = immse(R0_G,O);
mse_1_G = [immse(R1_G_1,O), immse(R1_G_2,O) ,immse(R1_G_3,O), immse(R1_G_4,O), immse(R1_G_5,O)];
mse_1_G_mean = mean(mse_1_G);

% figure
% subplot(3,6,1), imshow(R0_G,[]);title('Ramp');
% subplot(3,6,2), imshow(R1_G_1,[]);title('Ramp * Hamming ,cutoff = 0.25');
% subplot(3,6,3), imshow(R1_G_2,[]);title('Ramp * Hamming ,cutoff = 0.30');
% subplot(3,6,4), imshow(R1_G_3,[]);title('Ramp * Hamming ,cutoff = 0.35');
% subplot(3,6,5), imshow(R1_G_4,[]);title('Ramp * Hamming ,cutoff = 0.40'); 
% subplot(3,6,6), imshow(R1_G_5,[]);title('Ramp * Hamming ,cutoff = 0.45'); 
%% 2 Ramp * Gaussian filtering
R2_G_cutoff =[];
for i =1:length(cutoff)
Filter2=[];
filterline2=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
M= fir1(order(3)-1, cutoff(i), 'low', gausswin(order(3)));
[wgaussian,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wgaussian=abs((wgaussian)');
windowGaussian=zeros(1,N);
windowGaussian(N/2+1:N)=wgaussian(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowGaussian(i)=wgaussian(count);
    count=count+1;
end

filterline2 = abs( freqs ).* windowGaussian';


Filter2 = repmat(filterline2, [1 proj]);

% FT domain filtering
ft_sino_G= fftshift(fft(sino_G,[],1),1);
filteredProj = ft_sino_G .* Filter2;
filteredProj = ifftshift(filteredProj,1);
ift_sino_G = real(ifft(filteredProj,[],1));



% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_G(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, ramp * Gaussian window filter');
R2_G=Rsum;
R2_G_cutoff = [R2_G_cutoff,R2_G];
end
R2_G_1=R2_G_cutoff(:,1:N);
R2_G_2=R2_G_cutoff(:,N+1:2*N);
R2_G_3=R2_G_cutoff(:,2*N+1:3*N);
R2_G_4=R2_G_cutoff(:,3*N+1:4*N);
R2_G_5=R2_G_cutoff(:,4*N+1:5*N);

mse_2_G = [immse(R2_G_1,O), immse(R2_G_2,O) ,immse(R2_G_3,O), immse(R2_G_4,O), immse(R2_G_5,O)];
mse_2_G_mean = mean(mse_2_G);

% subplot(3,6,7), imshow(R0_G,[]);title('Ramp');
% subplot(3,6,8), imshow(R2_G_1,[]);title('Ramp * Gaussian ,cutoff = 0.25');
% subplot(3,6,9),imshow(R2_G_2,[]);title('Ramp * Gaussian ,cutoff = 0.30');
% subplot(3,6,10), imshow(R2_G_3,[]);title('Ramp * Gaussian ,cutoff = 0.35');
% subplot(3,6,11), imshow(R2_G_4,[]);title('Ramp * Gaussian ,cutoff = 0.40'); 
% subplot(3,6,12), imshow(R2_G_5,[]);title('Ramp * Gaussian ,cutoff = 0.45'); 

%% 3 Ramp * Hanning filtering
R3_G_cutoff=[],
for i =1:length(cutoff)
Filter3=[];
filterline3=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', hann(order(3)));
[whann,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
whann=abs((whann)');
windowHanning=zeros(1,N);
windowHanning(N/2+1:N)=whamming(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowHanning(i)=whann(count);
    count=count+1;
end
%plot(windowHanning);title('windowhanning');
filterline3 = abs( freqs ).* windowHanning';


Filter3 = repmat(filterline3, [1 proj]);

% FT domain filtering
ft_sino_G = fftshift(fft(sino_G,[],1),1);
filteredProj = ft_sino_G.* Filter3;
filteredProj = ifftshift(filteredProj,1);
ift_sino_G = real(ifft(filteredProj,[],1));



% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_G(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Hann window filter');
R3_G=Rsum;
R3_G_cutoff = [R3_G_cutoff,R3_G];
end
R3_G_1=R3_G_cutoff(:,1:N);
R3_G_2=R3_G_cutoff(:,N+1:2*N);
R3_G_3=R3_G_cutoff(:,2*N+1:3*N);
R3_G_4=R3_G_cutoff(:,3*N+1:4*N);
R3_G_5=R3_G_cutoff(:,4*N+1:5*N);

mse_3_G = [immse(R3_G_1,O), immse(R3_G_2,O) ,immse(R3_G_3,O),immse(R3_G_4,O), immse(R3_G_5,O)];
mse_3_G_mean = mean(mse_3_G);

% subplot(3,6,13), imshow(R0_G,[]);title('Ramp');
% subplot(3,6,14), imshow(R3_G_1,[]);title('Ramp * Hanning ,cutoff = 0.25');
% subplot(3,6,15),imshow(R3_G_2,[]);title('Ramp * Hanning ,cutoff = 0.30');
% subplot(3,6,16), imshow(R3_G_3,[]);title('Ramp * Hanning ,cutoff = 0.35');
% subplot(3,6,17), imshow(R3_G_4,[]);title('Ramp * Hanning ,cutoff = 0.40'); 
% subplot(3,6,18), imshow(R3_G_5,[]);title('Ramp * Hanning ,cutoff = 0.45'); 

%% 4 Ramp * Parzen filtering
R4_G_cutoff=[];
for i = 1:length(cutoff)
Filter4=[];
filterline4=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', parzenwin(order(3)));
[wparzen,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wparzen=abs((wparzen)');
windowParzen=zeros(1,N);
windowParzen(N/2+1:N)=wparzen(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowParzen(i)=wparzen(count);
    count=count+1;
end
filterline4 = abs( freqs ).* windowParzen';

Filter4 = repmat(filterline4, [1 proj]);

% FT domain filtering
ft_sino_G = fftshift(fft(sino_G,[],1),1);
filteredProj = ft_sino_G .* Filter4;
filteredProj = ifftshift(filteredProj,1);
ift_sino_G = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_G(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Parzen window filter');
R4_G=Rsum;
R4_G_cutoff = [R4_G_cutoff,R4_G];
end
R4_G_1=R4_G_cutoff(:,1:N);
R4_G_2=R4_G_cutoff(:,N+1:2*N);
R4_G_3=R4_G_cutoff(:,2*N+1:3*N);
R4_G_4=R4_G_cutoff(:,3*N+1:4*N);
R4_G_5=R4_G_cutoff(:,4*N+1:5*N);

mse_4_G = [immse(R4_G_1,O), immse(R4_G_2,O) ,immse(R4_G_3,O), immse(R4_G_4,O), immse(R4_G_5,O)];
mse_4_G_mean = mean(mse_4_G);

% figure
% subplot(3,6,1), imshow(R0_G,[]);title('Ramp');
% subplot(3,6,2), imshow(R4_G_1,[]);title('Ramp * Parzen ,cutoff = 0.25');
% subplot(3,6,3),imshow(R4_G_2,[]);title('Ramp * Parzen ,cutoff = 0.30');
% subplot(3,6,4), imshow(R4_G_3,[]);title('Ramp * Parzen ,cutoff = 0.35');
% subplot(3,6,5), imshow(R4_G_4,[]);title('Ramp * Parzen ,cutoff = 0.40'); 
% subplot(3,6,6), imshow(R4_G_5,[]);title('Ramp * Parzen ,cutoff = 0.45'); 

%% 5 Ramp * Bartlett filtering
R5_G_cutoff=[];
for i =1:length(cutoff)
Filter5=[];
filterline5=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', bartlett(order(3)));
[wbartlett,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wbartlett=abs((wbartlett)');
windowBartlett=zeros(1,N);
windowBartlett(N/2+1:N)=wbartlett(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowBartlett(i)=wbartlett(count);
    count=count+1;
end
% plot(windowBartlett);title('windowbartlett');
filterline5 = abs( freqs ).* windowBartlett';


Filter5 = repmat(filterline5, [1 proj]);

% FT domain filtering
ft_sino_G = fftshift(fft(sino_G,[],1),1);
filteredProj = ft_sino_G .* Filter5;
filteredProj = ifftshift(filteredProj,1);
ift_sino_G = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_G(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Bartlett window filter');
R5_G=Rsum;
R5_G_cutoff = [R5_G_cutoff,R5_G];
end
R5_G_1=R5_G_cutoff(:,1:N);
R5_G_2=R5_G_cutoff(:,N+1:2*N);
R5_G_3=R5_G_cutoff(:,2*N+1:3*N);
R5_G_4=R5_G_cutoff(:,3*N+1:4*N);
R5_G_5=R5_G_cutoff(:,4*N+1:5*N);

mse_5_G = [immse(R5_G_1,O), immse(R5_G_2,O) ,immse(R5_G_3,O), immse(R5_G_4,O), immse(R5_G_5,O)];
mse_5_G_mean = mean(mse_5_G);

% subplot(3,6,7), imshow(R0_G,[]);title('Ramp');
% subplot(3,6,8), imshow(R5_G_1,[]);title('Ramp * Barlett ,cutoff = 0.25');
% subplot(3,6,9),imshow(R5_G_2,[]);title('Ramp * Barlett ,cutoff = 0.30');
% subplot(3,6,10), imshow(R5_G_3,[]);title('Ramp * Barlett ,cutoff = 0.35');
% subplot(3,6,11), imshow(R5_G_4,[]);title('Ramp * Barlett ,cutoff = 0.40'); 
% subplot(3,6,12), imshow(R5_G_5,[]);title('Ramp * Barlett ,cutoff = 0.45'); 

%% 6 Ramp * Butterworth filtering
R6_G_cutoff=[];
for i =1:length(cutoff)
Filter6=[];
filterline6=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
[m,l]= butter(order(3)-2, cutoff(i), 'low');
[h, w]= freqz(m,l,N/2);
windowButterworth = abs(h)';
windowButterworth = [flip(windowButterworth) windowButterworth];
%plot(windowButterworth);title('windowbutterworth');


filterline6 = abs( freqs ).* (windowButterworth)';


Filter6 = repmat(filterline6, [1 proj]);

% FT domain filtering
ft_sino_G = fftshift(fft(sino_G,[],1),1);
filteredProj = ft_sino_G .* Filter6;
filteredProj = ifftshift(filteredProj,1);
ift_sino_G = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_G(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
%      
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Butterworth window filter');
R6_G_3=Rsum;
R6_G_cutoff = [R6_G_cutoff,R6_G_3];
end
R6_G_1=R6_G_cutoff(:,1:N);
R6_G_2=R6_G_cutoff(:,N+1:2*N);
R6_G_3=R6_G_cutoff(:,2*N+1:3*N);
R6_G_4=R6_G_cutoff(:,3*N+1:4*N);
R6_G_5=R6_G_cutoff(:,4*N+1:5*N);


mse_6_G = [immse(R6_G_1,O), immse(R6_G_2,O) ,immse(R6_G_3,O), immse(R6_G_4,O), immse(R6_G_5,O)];
mse_6_G_mean = mean(mse_6_G);
% subplot(3,6,13), imshow(R0_G,[]);title('Ramp');
% subplot(3,6,14), imshow(R6_G_1,[]);title('Ramp * Butter ,cutoff = 0.25');
% subplot(3,6,15),imshow(R6_G_2,[]);title('Ramp * Butter ,cutoff = 0.30');
% subplot(3,6,16), imshow(R6_G_3,[]);title('Ramp * Butter ,cutoff = 0.35');
% subplot(3,6,17), imshow(R6_G_4,[]);title('Ramp * Butter ,cutoff = 0.40'); 
% subplot(3,6,18), imshow(R6_G_5,[]);title('Ramp * Butter ,cutoff = 0.45'); 

%% Post Processing : Average Filtering for gaussian noise with cut off = 0.35

%Average filetring on cut off = 0.35
h_average=fspecial('average',2);
R0_G_average =imfilter(R0_G,h_average);
R1_G_3_average =imfilter(R1_G_3,h_average);
R2_G_3_average =imfilter(R2_G_3,h_average);
R3_G_3_average =imfilter(R3_G_3,h_average);
R4_G_3_average =imfilter(R4_G_3,h_average);
R5_G_3_average =imfilter(R5_G_3,h_average);
R6_G_3_average =imfilter(R6_G_3,h_average);
mse_0_average = immse(R0_G_average,O);

mse_average_3 = [immse(R1_G_3_average,O),immse(R2_G_3_average,O),immse(R3_G_3_average,O),immse(R4_G_3_average,O)]; %immse(R5_G_3_average,O),immse(R6_G_3_average,O)];

% %Average filetring on cut off = 0.4
% R1_G_4_average =imfilter(R1_G_4,h_average);
% R2_G_4_average =imfilter(R2_G_4,h_average);
% R3_G_4_average =imfilter(R3_G_4,h_average);
% R4_G_4_average =imfilter(R4_G_4,h_average);
% % R5_G_4_average =imfilter(R5_G_4,h_average);
% % R6_G_4_average =imfilter(R6_G_4,h_average);

% mse_average_4 = [immse(R1_G_4_average,O),immse(R2_G_4_average,O),immse(R3_G_4_average,O), immse(R4_G_4_average,O)]; %immse(R5_G_4_average,O),immse(R6_G_4_average,O)];

%figure, imshow(R0_G_average,[]);
% subplot(3,7,8),imshow(R1_G_4_average,[]), title('Average filtering Hamming 0.4');
% subplot(3,7,9),imshow(R2_G_4_average,[]), title('Average filtering Gaussian 0.4');
% subplot(3,7,10),imshow(R3_G_4_average,[]), title('Average filtering Hanning 0.4');
% subplot(3,7,11),imshow(R4_G_4_average,[]), title('Average filtering Parzen 0.4');
% % subplot(3,7,12),imshow(R5_G_4_average,[]), title('Average filtering Barlett 0.4');
% % subplot(3,7,13),imshow(R6_G_4_average,[]), title('Average filtering Butter 0.4');
% subplot(3,7,14),imshow(R0_G_average,[]), title('Average filtering Ramp');

figure,
subplot(4,2,1),imshow(R1_G_3_average,[]), title('Average filtering for Gaussian Noise: Hamming 0.35');
subplot(4,2,2),imshow(R2_G_3_average,[]), title('Average filtering Gaussian 0.35');
subplot(4,2,3),imshow(R3_G_3_average,[]), title('Average filtering Hanning 0.35');
subplot(4,2,4),imshow(R4_G_3_average,[]), title('Average filtering Parzen 0.35');
subplot(4,2,5),imshow(R5_G_3_average,[]), title('Average filtering Barlett 0.35');
subplot(4,2,6),imshow(R6_G_3_average,[]), title('Average filtering Butter 0.35');
subplot(4,2,7),imshow(R0_G_average,[]), title('Average filtering Ramp');

figure,
subplot(2,4,1),imshow(O),title('Original Image');
subplot(2,4,2),imshow(G),title('Original Image + Gaussian Noise');
subplot(2,4,3),imshow(R0_G,[]), title('filtering Ramp');
subplot(2,4,4),imshow(R0_G_average,[]), title('filtering Ramp + Average');
subplot(2,4,5),imshow(O),title('Original Image');
subplot(2,4,6),imshow(G),title('Original Image + Gaussian Noise');
subplot(2,4,7),imshow(R1_G_3,[]), title('Ramp + Hamming 0.35');
subplot(2,4,8),imshow(R1_G_3_average,[]), title('Ramp + Hammin 0.35 + Average');

%% Final Results for Gaussian Noise Reconstruction

% figure
% plot(freq_vector2,filterline);
% hold on
% plot(freq_vector2,filterline1);
% hold on
% plot(filterline2);
% hold on
% plot(filterline3);
% hold on
% plot(filterline4);
% hold on
% plot(filterline5);
% hold on
% plot(filterline6);
% legend('Ramp*Hamming Filter', 'Ramp*Gaussian Filter', 'Ramp*Hanning Filter', 'Ramp*Parzen Filter', 'Ramp*Bartlett Filter', 'Ramp*ButterWorth');

%Using cutoff frequency of 0.35
figure
subplot(2,1,1),imshow(R0_G,[]);title('Recontructed image for Gaussian Noise, Ramp window filter');
subplot(2,1,2), imagesc(sino_G);title('Non-filtered Sinogram');

figure
subplot(3,3,1), imshow(R1_G_3,[]);title('Recontructed images for Gaussian Noise, Ramp * Hamming window filter');
subplot(3,3,2), imshow(R2_G_3,[]);title('Ramp * Gaussian window filter');
subplot(3,3,3),imshow(R3_G_3,[]);title(' Ramp * Hanning window filter');
subplot(3,3,4), imshow(R4_G_3,[]);title(' Ramp * Parzen window filter');
subplot(3,3,5),imshow(R5_G_3,[]);title(' Ramp * Bartlett window filter');
subplot(3,3,6), imshow(R6_G_3,[]);title(' Ramp * Butterworth window filter');
subplot(3,3,7),imshow(R0_G,[]);title(' Ramp window filter');
 
%Normalization of Image Range - To Gray Scale
R1_G_3=abs(R1_G_3);
max_1G= max(max(R1_G_3));
R1_G_3=R1_G_3./max_1G;
% figure, imshow(R1_G);

R2_G_3=abs(R2_G_3);
max_2G= max(max(R2_G_3));
R2_G_3=R2_G_3./max_2G;
% figure, imshow(R2_G);

R3_G_3=abs(R3_G_3);
max_3G= max(max(R3_G_3));
R3_G_3=R3_G_3./max_3G;
% figure, imshow(R3_G);


R4_G_3=abs(R4_G_3);
max_4G= max(max(R4_G_3));
R4_G_3=R4_G_3./max_4G;
% figure, imshow(R4_G);

R5_G_3=abs(R5_G_3);
max_5G= max(max(R5_G_3));
R5_G_3=R5_G_3./max_5G;
% figure, imshow(R5_G);

R6_G_3=abs(R6_G_3);
max_6G= max(max(R6_G_3));
R6_G_3=R6_G_3./max_6G;
% figure, imshow(R6_G);

R0_G=abs(R0_G);
max_0G= max(max(R0_G));
R0_G=R0_G./max_0G;

%% Metrics

% Peak SNR
peaksnr=zeros(7);
peaksnr(1) = psnr(R0_G,O);
peaksnr(2) = psnr(R1_G_3,O);
peaksnr(3) = psnr(R2_G_3,O);
peaksnr(4) = psnr(R3_G_3,O);
peaksnr(5) = psnr(R4_G_3,O);
peaksnr(6) = psnr(R5_G_3,O);
peaksnr(7) = psnr(R6_G_3,O);
figure;
bar(peaksnr); title('peak SNR for Gaussian Noise Image Reconstruction');

% Mean square ratio
mse=zeros(7);
mse(1)=immse(R0_G,O);
mse(2)=immse(R1_G_3,O);
mse(3)=immse(R2_G_3,O);
mse(4)=immse(R3_G_3,O);
mse(5)=immse(R4_G_3,O);
mse(6)=immse(R5_G_3,O);
mse(7)=immse(R6_G_3,O);
figure;
bar(mse); title('MSE for Gaussian Noise Image Reconstruction');

% Structural similarity
ssimval=zeros(7);
ssimval(1)= ssim(R0_G,O);
ssimval(2)= ssim(R1_G_3,O);
ssimval(3)= ssim(R2_G_3,O);
ssimval(4)= ssim(R3_G_3,O);
ssimval(5)= ssim(R4_G_3,O);
ssimval(6)= ssim(R5_G_3,O);
ssimval(7)= ssim(R6_G_3,O);
figure
bar(ssimval); title('SSIM for Gaussian Noise Image Reconstruction');

%Final Results: Comparing Original, Ramp and Ramp*Hamming 
figure
subplot(2,2,1), imshow(O),title ('Original');
subplot(2,2,2), imshow(R0_G, []), title('Reconstructed for Gaussian Noise using Ramp Filter');
subplot(2,2,3), imshow(R1_G_3, []), title('Reconstructed for Gaussian Noise using Ramp*Hamming Filter');

%% Pepper and Salt Noise
P = imnoise(O,'salt & pepper')
figure
imshow(P);title('Pepper and Salt Noise');

% P_filtered = medfilt2(P,[5 5]);
% figure
% imshow(P_filtered);title('Filtered with median Added Noise');

%% PROJECTIONS
% proj= number of projections
proj=1026;
angle= (180/proj);
sino_P= zeros(N,proj);
for i=0:proj-1
   
    %rotation
    P_rot=imrotate(P,i*angle,'crop');
    %projection
    i=i+1;
    for j=1:N
        sino_P(j,i)= sum(P_rot(:,j));
    end
end
% figure
% imagesc(sino_P);title('Non-filtered SinoPram');

%% 0 Ramp filtering and reconstruction

% Creation of the ramp filter
freqs=linspace(-1, 1, N).';
filterline = abs( freqs );
Filter = repmat(filterline, [1 proj]);

% FT domain filtering
ft_sino_P = fftshift(fft(sino_P,[],1),1);
filteredProj = ft_sino_P .* Filter;
filteredProj = ifftshift(filteredProj,1);
ift_sino_P = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_P(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Reconstructed image, ramp filter');
R0_P= Rsum;

%% 1 Ramp * Hamming filtering
R1_P_cutoff=[];

for i =1:length(cutoff)
    Filter1=[];
    filterline1=[];
    % Creation of the window*ramp filter
    freqs=linspace(-1, 1, N).';
    m= fir1(order(3)-1, cutoff(i), 'low', hamming(order(3)));
    [whamming,w]=freqz(m,1,N/2);
    % figure
    % plot(abs(whamming))
    whamming=abs((whamming)');
    windowHamming=zeros(1,N);
    windowHamming(N/2+1:N)=whamming(1:N/2);
    count=1;
    i=0;
    for i=N/2:-1:1
        windowHamming(i)=whamming(count);
        count=count+1;
    end
    filterline1 = abs( freqs ).* windowHamming';
    
    
    Filter1 = repmat(filterline1, [1 proj]);
    
    % FT domain filtering
    ft_sino_P = fftshift(fft(sino_P,[],1),1);
    filteredProj = ft_sino_P .* Filter1;
    filteredProj = ifftshift(filteredProj,1);
    ift_sino_P = real(ifft(filteredProj,[],1));
    
    
    % Retroprojection
    R=zeros(N,N);
    Rsum=zeros(N,N);
    for i=1:proj
        for j= 1:N
            R(:,j)= ift_sino_P(j,i);
        end
        R= imrotate(R,-i*angle,'crop');
        Rsum=Rsum+R;
    end
    
    % figure
    % imshow(Rsum,[]);title('Recontructed image, ramp * Hamming window filter');
    R1_P= Rsum;
    R1_P_cutoff = [R1_P_cutoff,R1_P];
end
R1_P_1=R1_P_cutoff(:,1:N);
R1_P_2=R1_P_cutoff(:,N+1:2*N);
R1_P_3=R1_P_cutoff(:,2*N+1:3*N);
R1_P_4=R1_P_cutoff(:,3*N+1:4*N);
R1_P_5=R1_P_cutoff(:,4*N+1:5*N);

mse_P = immse(P,O);
mse_0 = immse(R0_P,O);
mse_1_P = [immse(R1_P_1,O), immse(R1_P_2,O) ,immse(R1_P_3,O),immse(R1_P_4,O), immse(R1_P_5,O)];
mse_1_P_mean = mean(mse_1_P);

% figure
% subplot(3,6,1), imshow(R0_P,[]);title('Ramp');
% subplot(3,6,2), imshow(R1_P_1,[]);title('Ramp * Hamming ,cutoff = 0.25');
% subplot(3,6,3), imshow(R1_P_2,[]);title('Ramp * Hamming ,cutoff = 0.3');
% subplot(3,6,4), imshow(R1_P_3,[]);title('Ramp * Hamming ,cutoff = 0.35');
% % subplot(3,6,5), imshow(R1_P_4,[]);title('Ramp * Hamming ,cutoff = 0.40'); 
% subplot(3,6,6), imshow(R1_P_5,[]);title('Ramp * Hamming ,cutoff = 0.45'); 
%% 2 Ramp * Gaussian filtering
R2_P_cutoff =[];
for i =1:length(cutoff)
Filter2=[];
filterline2=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
M= fir1(order(3)-1, cutoff(i), 'low', gausswin(order(3)));
[wgaussian,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wgaussian=abs((wgaussian)');
windowGaussian=zeros(1,N);
windowGaussian(N/2+1:N)=wgaussian(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowGaussian(i)=wgaussian(count);
    count=count+1;
end

filterline2 = abs( freqs ).* windowGaussian';


Filter2 = repmat(filterline2, [1 proj]);

% FT domain filtering
ft_sino_P= fftshift(fft(sino_P,[],1),1);
filteredProj = ft_sino_P .* Filter2;
filteredProj = ifftshift(filteredProj,1);
ift_sino_P = real(ifft(filteredProj,[],1));



% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_P(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, ramp * Paussian window filter');
R2_P=Rsum;
R2_P_cutoff = [R2_P_cutoff,R2_P];
end
R2_P_1=R2_P_cutoff(:,1:N);
R2_P_2=R2_P_cutoff(:,N+1:2*N);
R2_P_3=R2_P_cutoff(:,2*N+1:3*N);
R2_P_4=R2_P_cutoff(:,3*N+1:4*N);
R2_P_5=R2_P_cutoff(:,4*N+1:5*N);

mse_2_P = [immse(R2_P_1,O), immse(R2_P_2,O) ,immse(R2_P_3,O), immse(R2_P_4,O), immse(R2_P_5,O)];
mse_2_P_mean = mean(mse_2_P);

% subplot(3,6,7), imshow(R0_P,[]);title('Ramp');
% subplot(3,6,8), imshow(R2_P_1,[]);title('Ramp * Paussian ,cutoff = 0.25');
% subplot(3,6,9),imshow(R2_P_2,[]);title('Ramp * Gaussian ,cutoff = 0.30');
% subplot(3,6,10), imshow(R2_P_3,[]);title('Ramp * Gaussian ,cutoff = 0.35');
% subplot(3,6,11), imshow(R2_P_4,[]);title('Ramp * Gaussian ,cutoff = 0.4'); 
% subplot(3,6,12), imshow(R2_P_5,[]);title('Ramp * Gaussian ,cutoff = 0.45'); 

%% 3 Ramp * Hanning filtering
R3_P_cutoff=[],
for i =1:length(cutoff)
Filter3=[];
filterline3=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', hann(order(3)));
[whann,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
whann=abs((whann)');
windowHanning=zeros(1,N);
windowHanning(N/2+1:N)=whamming(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowHanning(i)=whann(count);
    count=count+1;
end
%plot(windowHanning);title('windowhanning');
filterline3 = abs( freqs ).* windowHanning';


Filter3 = repmat(filterline3, [1 proj]);

% FT domain filtering
ft_sino_P = fftshift(fft(sino_P,[],1),1);
filteredProj = ft_sino_P.* Filter3;
filteredProj = ifftshift(filteredProj,1);
ift_sino_P = real(ifft(filteredProj,[],1));



% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_P(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Hann window filter');
R3_P=Rsum;
R3_P_cutoff = [R3_P_cutoff,R3_P];
end
R3_P_1=R3_P_cutoff(:,1:N);
R3_P_2=R3_P_cutoff(:,N+1:2*N);
R3_P_3=R3_P_cutoff(:,2*N+1:3*N);
R3_P_4=R3_P_cutoff(:,3*N+1:4*N);
R3_P_5=R3_P_cutoff(:,4*N+1:5*N);

mse_3_P = [immse(R3_P_1,O), immse(R3_P_2,O) ,immse(R3_P_3,O),immse(R3_P_4,O), immse(R3_P_5,O)];
mse_3_P_mean = mean(mse_3_P);

% subplot(3,6,13), imshow(R0_P,[]);title('Ramp');
% subplot(3,6,14), imshow(R3_P_1,[]);title('Ramp * Hanning ,cutoff = 0.25');
% subplot(3,6,15),imshow(R3_P_2,[]);title('Ramp * Hanning ,cutoff = 0.30');
% subplot(3,6,16), imshow(R3_P_3,[]);title('Ramp * Hanning ,cutoff = 0.35');
% % subplot(3,6,17), imshow(R3_P_4,[]);title('Ramp * Hanning ,cutoff = 0.40'); 
% subplot(3,6,18), imshow(R3_P_5,[]);title('Ramp * Hanning ,cutoff = 0.45'); 

%% 4 Ramp * Parzen filtering
R4_P_cutoff=[];
for i = 1:length(cutoff)
Filter4=[];
filterline4=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', parzenwin(order(3)));
[wparzen,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wparzen=abs((wparzen)');
windowParzen=zeros(1,N);
windowParzen(N/2+1:N)=wparzen(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowParzen(i)=wparzen(count);
    count=count+1;
end
filterline4 = abs( freqs ).* windowParzen';

Filter4 = repmat(filterline4, [1 proj]);

% FT domain filtering
ft_sino_P = fftshift(fft(sino_P,[],1),1);
filteredProj = ft_sino_P .* Filter4;
filteredProj = ifftshift(filteredProj,1);
ift_sino_P = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_P(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Parzen window filter');
R4_P=Rsum;
R4_P_cutoff = [R4_P_cutoff,R4_P];
end
R4_P_1=R4_P_cutoff(:,1:N);
R4_P_2=R4_P_cutoff(:,N+1:2*N);
R4_P_3=R4_P_cutoff(:,2*N+1:3*N);
R4_P_4=R4_P_cutoff(:,3*N+1:4*N);
R4_P_5=R4_P_cutoff(:,4*N+1:5*N);

mse_4_P = [immse(R4_P_1,O), immse(R4_P_2,O) ,immse(R4_P_3,O), immse(R4_P_4,O), immse(R4_P_5,O)];
mse_4_P_mean = mean(mse_4_P);

% figure
% subplot(3,6,1), imshow(R0_P,[]);title('Ramp');
% subplot(3,6,2), imshow(R4_P_1,[]);title('Ramp * Parzen ,cutoff = 0.25');
% subplot(3,6,3),imshow(R4_P_2,[]);title('Ramp * Parzen ,cutoff = 0.30');
% subplot(3,6,4), imshow(R4_P_3,[]);title('Ramp * Parzen ,cutoff = 0.35');
% subplot(3,6,5), imshow(R4_P_4,[]);title('Ramp * Parzen ,cutoff = 0.40'); 
% subplot(3,6,6), imshow(R4_P_5,[]);title('Ramp * Parzen ,cutoff = 0.45'); 

%% 5 Ramp * Bartlett filtering
R5_P_cutoff=[];
for i =1:length(cutoff)
Filter5=[];
filterline5=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
m= fir1(order(3)-1, cutoff(i), 'low', bartlett(order(3)));
[wbartlett,w]=freqz(m,1,N/2);
% figure
% plot(abs(whamming))
wbartlett=abs((wbartlett)');
windowBartlett=zeros(1,N);
windowBartlett(N/2+1:N)=wbartlett(1:N/2);
count=1;
i=0;
for i=N/2:-1:1
    windowBartlett(i)=wbartlett(count);
    count=count+1;
end
% plot(windowBartlett);title('windowbartlett');
filterline5 = abs( freqs ).* windowBartlett';


Filter5 = repmat(filterline5, [1 proj]);

% FT domain filtering
ft_sino_P = fftshift(fft(sino_P,[],1),1);
filteredProj = ft_sino_P .* Filter5;
filteredProj = ifftshift(filteredProj,1);
ift_sino_P = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_P(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
     
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Bartlett window filter');
R5_P=Rsum;
R5_P_cutoff = [R5_P_cutoff,R5_P];
end
R5_P_1=R5_P_cutoff(:,1:N);
R5_P_2=R5_P_cutoff(:,N+1:2*N);
R5_P_3=R5_P_cutoff(:,2*N+1:3*N);
R5_P_4=R5_P_cutoff(:,3*N+1:4*N);
R5_P_5=R5_P_cutoff(:,4*N+1:5*N);

mse_5_P = [immse(R5_P_1,O), immse(R5_P_2,O) ,immse(R5_P_3,O), immse(R5_P_4,O), immse(R5_P_5,O)];
mse_5_P_mean = mean(mse_5_P);

% subplot(3,6,7), imshow(R0_P,[]);title('Ramp');
% subplot(3,6,8), imshow(R5_P_1,[]);title('Ramp * Barlett ,cutoff = 0.25');
% subplot(3,6,9),imshow(R5_P_2,[]);title('Ramp * Barlett ,cutoff = 0.30');
% subplot(3,6,10), imshow(R5_P_3,[]);title('Ramp * Barlett ,cutoff = 0.35');
% subplot(3,6,11), imshow(R5_P_4,[]);title('Ramp * Barlett ,cutoff = 0.40'); 
% subplot(3,6,12), imshow(R5_P_5,[]);title('Ramp * Barlett ,cutoff = 0.45'); 

%% 6 Ramp * Butterworth filtering
R6_P_cutoff=[];
for i =1:length(cutoff)
Filter6=[];
filterline6=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
[m,l]= butter(order(3)-2, cutoff(i), 'low');
[h, w]= freqz(m,l,N/2);
windowButterworth = abs(h)';
windowButterworth = [flip(windowButterworth) windowButterworth];
%plot(windowButterworth);title('windowbutterworth');


filterline6 = abs( freqs ).* (windowButterworth)';


Filter6 = repmat(filterline6, [1 proj]);

% FT domain filtering
ft_sino_P = fftshift(fft(sino_P,[],1),1);
filteredProj = ft_sino_P .* Filter6;
filteredProj = ifftshift(filteredProj,1);
ift_sino_P = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_P(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
%      
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Butterworth window filter');
R6_P=Rsum;
R6_P_cutoff = [R6_P_cutoff,R6_P];
end
R6_P_1=R6_P_cutoff(:,1:N);
R6_P_2=R6_P_cutoff(:,N+1:2*N);
R6_P_3=R6_P_cutoff(:,2*N+1:3*N);
R6_P_4=R6_P_cutoff(:,3*N+1:4*N);
R6_P_5=R6_P_cutoff(:,4*N+1:5*N);


mse_6_P = [immse(R6_P_1,O), immse(R6_P_2,O) ,immse(R6_P_3,O), immse(R6_P_4,O), immse(R6_P_5,O)];
mse_6_P_mean = mean(mse_6_P);
% subplot(3,6,13), imshow(R0_P,[]);title('Ramp');
% subplot(3,6,14), imshow(R6_P_1,[]);title('Ramp * Butter ,cutoff = 0.25');
% subplot(3,6,15),imshow(R6_P_2,[]);title('Ramp * Butter ,cutoff = 0.30');
% subplot(3,6,16), imshow(R6_P_3,[]);title('Ramp * Butter ,cutoff = 0.35');
% subplot(3,6,17), imshow(R6_P_4,[]);title('Ramp * Butter ,cutoff = 0.40'); 
% subplot(3,6,18), imshow(R6_P_5,[]);title('Ramp * Butter ,cutoff = 0.45'); 
%%
    
% figure
% plot(filterline);
% hold on
% plot(filterline1);
% hold on
% plot(filterline2);
% hold on
% plot(filterline3);
% hold on
% plot(filterline4);
% hold on
% plot(filterline5);
% hold on
% plot(filterline6);
% legend('Ramp', 'Ramp*Hamming Filter', 'Ramp*Gaussian Filter', 'Ramp*Hanning Filter', 'Ramp*Parzen Filter', 'Ramp*Bartlett Filter', 'Ramp*ButterWorth');

%Values taken at cutofffreq= 0.35

R1_P_3=abs(R1_P_3);
max_1P= max(max(R1_P_3));
R1_P_3=R1_P_3./max_1P;
% figure, imshow(R1_P);

R2_P_3=abs(R2_P_3);
max_2P= max(max(R2_P_3));
R2_P_3=R2_P_3./max_2P;
% figure, imshow(R2_P);

R3_P_3=abs(R3_P_3);
max_3P= max(max(R3_P_3));
R3_P_3=R3_P_3./max_3P;
% figure, imshow(R3_P);


R4_P_3=abs(R4_P_3);
max_4P= max(max(R4_P_3));
R4_P_3=R4_P_3./max_4P;
% figure, imshow(R4_P);

R5_P_3=abs(R5_P_3);
max_5P= max(max(R5_P_3));
R5_P_3=R5_P_3./max_5P;
% figure, imshow(R5_P);

R6_P_3=abs(R6_P_3);
max_6P= max(max(R6_P_3));
R6_P_3=R6_P_3./max_6P;
% figure, imshow(R6_P);

R0_P=abs(R0_P);
max_0P= max(max(R0_P));
R0_P=R0_P./max_0P;

figure
subplot(2,1,1),imshow(R0_P, []);title('Recontructed PS Noisy image, Ramp window filter');
subplot(2,1,2), imagesc(sino_P); xlabel('projection number'), ylabel('Row Number'),,title('Non-filtered Sinogram');

% angle_vect=(0:180/1026:180);
figure
subplot(4,2,1),imshow(R0_P,[]);title('Recontructed PS Noisy image, Ramp window filter');
subplot(4,2,2), imshow(R1_P_3,[]);title('Recontructed image, Ramp * Hamming window filter');
subplot(4,2,3), imshow(R2_P_3,[]);title('Recontructed image, Ramp * Gaussian window filter');
subplot(4,2,4),imshow(R3_P_3,[]);title('Recontructed image, Ramp * Hanning window filter');
subplot(4,2,5), imshow(R4_P_3,[]);title('Recontructed image, Ramp * Parzen window filter');
subplot(4,2,6),imshow(R5_P_3,[]);title('Recontructed image, Ramp * Bartlett window filter');
subplot(4,2,7), imshow(R6_P_3,[]);title('Recontructed image, Ramp * Butterworth window filter');

P_filtered_R0 = medfilt2(R0_P,[3 3], 'symmetric')

P_filtered_R1 = medfilt2(R1_P_3,[3 3], 'symmetric')

P_filtered_R2 = medfilt2(R2_P_3,[3 3], 'symmetric')

P_filtered_R3 = medfilt2(R3_P_3,[3 3], 'symmetric')

P_filtered_R4 = medfilt2(R4_P_3, [3 3], 'symmetric')

P_filtered_R5 = medfilt2(R5_P_3,[3 3], 'symmetric')

P_filtered_R6 = medfilt2(R6_P_3,[3 3])

%2nd iteration
P_filtered_R0 = medfilt2(P_filtered_R0,[3 3])

P_filtered_R1 = medfilt2(P_filtered_R1,[3 3])

P_filtered_R2 = medfilt2(P_filtered_R2,[3 3])

P_filtered_R3 = medfilt2(P_filtered_R3,[3 3])

P_filtered_R4 = medfilt2(P_filtered_R4 ,[3 3])

P_filtered_R5 = medfilt2(P_filtered_R5,[3 3])

P_filtered_R6 = medfilt2(P_filtered_R6,[3 3])

figure
subplot(3,3,1), imshow(P_filtered_R1,[]);title('Recontructed images after Median Filter, Ramp * Hamming window filter');
subplot(3,3,2), imshow(P_filtered_R2,[]);title(' Ramp * Hanning window filter');
subplot(3,3,3),imshow(P_filtered_R3,[]);title('Ramp * Gaussian window filter');
subplot(3,3,4), imshow(P_filtered_R4,[]);title('Ramp * Parzen window filter');
subplot(3,3,5),imshow(P_filtered_R5,[]);title('Ramp * Bartlett window filter');
subplot(3,3,6), imshow(P_filtered_R6,[]);title('Ramp * Butterworth window filter');
subplot(3,3,7), imshow(P_filtered_R0,[]);title(' Ramp filter');


% %% Enhance Edges and superimpose 
% %Go get the edges of the original image and superimpose to the
% %reconstructed image
% hy = [0 1 0; 0 0 0; 0 -1 0]; % y direction
% hx = hy'; % x direction
% 
% %derivative computation in the two directions
% Gx = imfilter(O, hx, 'conv');
% Gy = imfilter(O, hy, 'conv');
% 
% figure
% subplot(121), imshow(Gx, []), title('Gx'), colorbar
% subplot(122), imshow(Gy, []), title('Gy'), colorbar
% 
% %sum of absolute partial derivative values
% O_gradient = abs(Gx)+abs(Gy);
% 
% threshold = linspace(0, max(O(:))/2, 10);
% 
% O_edge = zeros(size(O,1), size(O,2), 10);
% figure
% for i=1:length(threshold)
%     O_edge( :,:,i ) = imbinarize(O_gradient, threshold(i));
%     subplot(2,5,i), imshow(O_edge(:,:,i), [])
%     title( strcat('thr = ', num2str(threshold(i))) )
%     drawnow;
% end
% 
% thr = threshold(2)              % optimal threshold 
% O_thr = O_edge(:, :, 2);
% figure
% imshow(O_thr)
% Superimpose the edgest to the Reconstructed Image
% for i=1:512
%     for j=1:512
%         if(O_thr(i,j) == 1)
%             R1_P_3(i,j)=O_thr(i,j);
%         end
%     end
% end

% figure
% imshow(R1_P_3,'InitialMagnification', 'fit')
% drawnow;

%% Metrics for Peper and Salt Noise

%Reconstructed Image 
% Peak SNR
peaksnr=zeros(7);
peaksnr(1) = psnr(R0_P,O);
peaksnr(2) = psnr(R1_P_3,O);
peaksnr(3) = psnr(R2_P_3,O);
peaksnr(4) = psnr(R3_P_3,O);
peaksnr(5) = psnr(R4_P_3,O);
peaksnr(6) = psnr(R5_P_3,O);
peaksnr(7) = psnr(R6_P_3,O);
figure
bar(peaksnr); title('peak SNR for Peper and Salt Noisy Image Reconstruction ');

% Mean square ratio
mse=zeros(7);
mse(1)=immse(R0_P,O);
mse(2)=immse(R1_P_3,O);
mse(3)=immse(R2_P_3,O);
mse(4)=immse(R3_P_3,O);
mse(5)=immse(R4_P_3,O);
mse(6)=immse(R5_P_3,O);
mse(7)=immse(R6_P_3,O);
figure
bar(mse); title('MSE for Peper and Salt Noisy Image Reconstruction ');

% Structural similarity
ssimval=zeros(7);
ssimval(1)= ssim(R0_P,O);
ssimval(2)= ssim(R1_P_3,O);
ssimval(3)= ssim(R2_P_3,O);
ssimval(4)= ssim(R3_P_3,O);
ssimval(5)= ssim(R4_P_3,O);
ssimval(6)= ssim(R5_P_3,O);
ssimval(7)= ssim(R6_P_3,O);
figure
bar(ssimval); title('SSIM for Peper and Salt Noisy Image Reconstruction ');

%After Post-processing with median filter 
% Peak SNR
peaksnr=zeros(7);
peaksnr(1) = psnr(P_filtered_R0,O);
peaksnr(2) = psnr(P_filtered_R1,O);
peaksnr(3) = psnr(P_filtered_R2,O);
peaksnr(4) = psnr(P_filtered_R3,O);
peaksnr(5) = psnr(P_filtered_R4,O);
peaksnr(6) = psnr(P_filtered_R5,O);
peaksnr(7) = psnr(P_filtered_R6,O);
figure
bar(peaksnr); title('peak SNR for Peper and Salt Noisy Image Reconstruction and Applied Median Filter');

% Mean square ratio
mse=zeros(7);
mse(1)=immse(P_filtered_R0,O);
mse(2)=immse(P_filtered_R1,O);
mse(3)=immse(P_filtered_R2,O);
mse(4)=immse(P_filtered_R3,O);
mse(5)=immse(P_filtered_R4,O);
mse(6)=immse(P_filtered_R5,O);
mse(7)=immse(P_filtered_R6,O);
figure
bar(mse); title('MSE for Peper and Salt Noisy Image Reconstruction and Applied Median Filter');

% Structural similarity
ssimval=zeros(7);
ssimval(1)= ssim(P_filtered_R0,O);
ssimval(2)= ssim(P_filtered_R1,O);
ssimval(3)= ssim(P_filtered_R2,O);
ssimval(4)= ssim(P_filtered_R3,O);
ssimval(5)= ssim(P_filtered_R4,O);
ssimval(6)= ssim(P_filtered_R5,O);
ssimval(7)= ssim(P_filtered_R6,O);
figure
bar(ssimval); title('SSIM for Peper and Salt Noisy Image Reconstruction and Applied Median Filter');

%Final Results: Comparing Original, Noisy, Ramp, Ramp*Hamming and after
%Median Filter, Ramp and Ramp*Hamming

figure
subplot(2,3,1), imshow(O), title('Original');
subplot(2,3,2), imshow(P), title('Peper and Salt Noise');
subplot(2,3,3), imshow(R0_P), title('Ramp Filter');
subplot(2,3,4), imshow(R1_P_3),title('Ramp*Hamming Filter');
subplot(2,3,5), imshow(P_filtered_R0), title('Ramp Filter + Median Filter');
subplot(2,3,6), imshow(P_filtered_R1), title('Ramp* Hamming Filter + Median Filter');

%% CT Chest Image 
load ('imChest.mat');
ct = im;
figure
imshow(ct);
N_ct = length(ct(:,1));
%% PROJECTIONS
% proj= number of projections

proj=512;
angle= (180/proj);
sino_ct= zeros(N_ct,proj);
for i=0:proj-1
    %rotation
    ct_rot=imrotate(ct,i*angle,'crop');
    %projection
    i=i+1;
    for j=1:N_ct
        sino_ct(j,i)= sum(ct_rot(:,j));
    end
end
% figure
% imagesc(sino_ct), colormap gray ,title('Non-filtered Sinogram');

%% 0 Ramp filtering and reconstruction of CT

% Creation of the ramp filter
freqs=linspace(-1, 1, N_ct).';
filterline = abs( freqs );
Filter = repmat(filterline, [1 proj]);

% FT domain filtering
ft_sino_ct = fftshift(fft(sino_ct,[],1),1);
% figure,plot(freqs,abs(ft_sino_ct));
filteredProj_ct = ft_sino_ct .* Filter;
filteredProj_ct = ifftshift(filteredProj_ct,1);
ift_sino_ct = real(ifft(filteredProj_ct,[],1));

% Sinogram in the frequency domain
% imagesc(abs(ft_sino)); title('FFT of the sinogram');
% Filtered sinogram
% imagesc(ift_sino);
% title('filtered sinogram');


% Retroprojection
R=zeros(N_ct,N_ct);
Rsum=zeros(N_ct,N_ct);
for i=1:proj
    for j= 1:N_ct
        R(:,j)= ift_sino_ct(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end   
% figure
% imshow(Rsum,[]);title('Reconstructed image, ramp filter');
R0_ct= Rsum;

%% 1 Ramp * Hamming filtering

Filter1=[];
filterline1=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N_ct).';
m= fir1(order(3)-1, 0.6, 'low', hamming(order(3)));
[whamming,w]=freqz(m,1,N_ct/2);
% figure
% plot(abs(whamming))
whamming=abs((whamming)');
windowHamming=zeros(1,N_ct);
windowHamming(N_ct/2+1:N_ct)=whamming(1:N_ct/2);
count=1;
i=0;
for i=N_ct/2:-1:1
    windowHamming(i)=whamming(count);
    count=count+1;
end
filterline1 = abs( freqs ).* windowHamming';


Filter1 = repmat(filterline1, [1 proj]);

% FT domain filtering
ft_sino_ct = fftshift(fft(sino_ct,[],1),1);
filteredProj = ft_sino_ct .* Filter1;
filteredProj = ifftshift(filteredProj,1);
ift_sino_ct = real(ifft(filteredProj,[],1));


% Retroprojection
R=zeros(N_ct,N_ct);
Rsum=zeros(N_ct,N_ct);
for i=1:proj
    for j= 1:N_ct
        R(:,j)= ift_sino_ct(j,i);
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end

% figure
% imshow(Rsum,[]);title('Recontructed image, ramp * Hamming window filter');
R1_ct= Rsum;


%% 2 Ramp * Gaussian filtering
Filter2=[];
filterline2=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N_ct).';
M= fir1(order(3)-1, cutoff(2), 'low', gausswin(order(3)));
[wgaussian,w]=freqz(m,1,N_ct/2);
% figure
% plot(abs(wgaussian))
wgaussian=abs((wgaussian)');
windowGaussian=zeros(1,N_ct);
windowGaussian(N_ct/2+1:N_ct)=wgaussian(1:N_ct/2);
count=1;
i=0;
for i=N_ct/2:-1:1
    windowGaussian(i)=wgaussian(count);
    count=count+1;
end

filterline2 = abs( freqs ).* windowGaussian';


Filter2 = repmat(filterline2, [1 proj]);

% FT domain filtering
ft_sino_ct= fftshift(fft(sino_ct,[],1),1);
filteredProj = ft_sino_ct .* Filter2;
filteredProj = ifftshift(filteredProj,1);
ift_sino_ct = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N_ct,N_ct);
Rsum=zeros(N_ct,N_ct);
for i=1:proj
    for j= 1:N_ct
        R(:,j)= ift_sino_ct(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
% figure
% imshow(Rsum,[]);title('Recontructed image, ramp * Gaussian window filter');
R2_ct=Rsum;

%% 6 Ramp * Butterworth filtering

Filter6=[];
filterline6=[];
% Creation of the window*ramp filter
freqs=linspace(-1, 1, N).';
[m,l]= butter(order(2)-2, cutoff(2), 'low');
[h, w]= freqz(m,l,N/2);
windowButterworth = abs(h)';
windowButterworth = [flip(windowButterworth) windowButterworth];
%plot(windowButterworth);title('windowbutterworth');

filterline6 = abs( freqs ).* (windowButterworth)';

Filter6 = repmat(filterline6, [1 proj]);

% FT domain filtering
ft_sino_ct = fftshift(fft(sino_ct,[],1),1);
filteredProj = ft_sino_ct .* Filter6;
filteredProj = ifftshift(filteredProj,1);
ift_sino_ct = real(ifft(filteredProj,[],1));

% Retroprojection
R=zeros(N,N);
Rsum=zeros(N,N);
for i=1:proj
    for j= 1:N
        R(:,j)= ift_sino_ct(j,i);  
    end
    R= imrotate(R,-i*angle,'crop');
    Rsum=Rsum+R;
end
%      
% figure
% imshow(Rsum,[]);title('Recontructed image, Ramp * Butterworth window filter');
R6_ct=Rsum;

%%
    
figure
plot(freq_vector2, filterline);
hold on
plot(freq_vector2, filterline1);
hold on
plot(freq_vector2, filterline2);
hold on
plot(freq_vector2, filterline6);
legend('Ramp', 'Ramp*Hamming Filter', 'Ramp*Gaussian Filter', 'Ramp*ButterWorth');

%Values taken at cutofffreq= 0.5

R1_ct=abs(R1_ct);
max_1ct= max(max(R1_ct));
R1_ct=R1_ct./max_1ct;
% figure, imshow(R1_P);

R2_ct=abs(R2_ct);
max_2ct= max(max(R2_ct));
R2_ct=R2_ct./max_2ct;
% figure, imshow(R2_P);

R6_ct=abs(R6_ct);
max_6ct= max(max(R6_ct));
R6_ct=R6_ct./max_6ct;
% figure, imshow(R6_P);

R0_ct=abs(R0_ct);
max_0ct= max(max(R0_ct));
R0_ct=R0_ct./max_0ct;

figure
subplot(2,1,1),imshow(R0_ct, []);title('Recontructed  CT image, Ramp window filter');
subplot(2,1,2), imagesc(sino_ct);title('Non-filtered Sinogram for CT image');


figure
subplot(2,2,1),imshow(R0_ct,[]);title('Recontructed CT image, Ramp window filter');
subplot(2,2,2), imshow(R1_ct,[]);title('Recontructed CT image, Ramp * Hamming window filter');
subplot(2,2,3), imshow(R2_ct,[]);title('Recontructed CT image, Ramp * Gaussian window filter');
subplot(2,2,4), imshow(R6_ct,[]);title('Recontructed CT image, Ramp * Butterworth window filter');

%% Performance Metrics
%Reconstructed Image 
% Peak SNR
peaksnr=zeros(4);
peaksnr(1) = psnr(R0_ct,ct);
peaksnr(2) = psnr(R1_ct,ct);
peaksnr(3) = psnr(R2_ct,ct);
peaksnr(4) = psnr(R6_ct,ct);
figure
bar(peaksnr); title('peak SNR for CT Image');

% Mean square ratio
mse=zeros(4);
mse(1)=immse(R0_ct,ct);
mse(2)=immse(R1_ct,ct);
mse(3)=immse(R2_ct,ct);
mse(4)=immse(R6_ct,ct);
figure
bar(mse); title('MSE for CT image');

% Structural similarity
ssimval=zeros(4);
ssimval(1)= ssim(R0_ct,ct);
ssimval(2)= ssim(R1_ct,ct);
ssimval(3)= ssim(R2_ct,ct);
ssimval(4)= ssim(R6_ct,ct);
figure
bar(ssimval); title('SSIM for CT image');





