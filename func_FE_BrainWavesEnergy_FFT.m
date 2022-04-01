function [delta,theta,alpha,beta,gamma,power]=func_FE_BrainWavesEnergy_FFT(data,Fs,opt)

%%---------------------------------------------------------------------------------------------------------
%%% This function extracts five brain wave signals in time domain and also 
%%% the power of each band usning FFT algorithm

%%% input: data (channels,samples)
%%%        Fs  frequecy rate (Hz)
%%%        opt  1-> show the results   0->do not show the results

%%% output: delta -> delta band in time domain (channels,samples)
%%%         theta -> theta band in time domain (channels,samples)
%%%         alpha -> alpha band in time domain (channels,samples)
%%%         beta -> beta band in time domain (channels,samples)
%%%         gamma -> gamma band in time domain (channels,samples)
%%%         power -> brain waves power (channels,power of each band [delta, theta, alpha, beta, gamma])

%%% writen by Golnaz Baghdadi 6/23/2021
%%---------------------------------------------------------------------------------------------------------


[channel,L]=size(data);
time=linspace(0,fix(L/Fs),L);
lenTime=length(time);

    for i=1:channel
        x=data(i,:);       % Your signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        X = fft(x,NFFT)/L;
        f = Fs/2*linspace(0,1,NFFT/2+1);
        fLen=min(length(X),lenTime);

        %%% ---------- frequecy bands extraction -----------%
        % delta
        Xf=zeros(1,NFFT/2+1);
        Xf(find(f>=0.5 & f<4))=X(find(f>=0.5 & f<4))*L;
        delta(i,:)=ifft(Xf,NFFT,'symmetric');

        % theta
        Xf=zeros(1,NFFT/2+1);
        Xf(find(f>=4 & f<8))=X(find(f>=4 & f<8))*L;
        theta(i,:)=ifft(Xf,NFFT,'symmetric');

        % alpha
        Xf=zeros(1,NFFT/2+1);
        Xf(find(f>=8 & f<13))=X(find(f>=8 & f<13))*L;
        alpha(i,:)=ifft(Xf,NFFT,'symmetric');

        % beta
        Xf=zeros(1,NFFT/2+1);
        Xf(find(f>=13 & f<30))=X(find(f>=13 & f<30))*L;
        beta(i,:)=ifft(Xf,NFFT,'symmetric');


        % gamma
        Xf=zeros(1,NFFT/2+1);
        Xf(find(f>=30 & f<80))=X(find(f>=30 & f<80))*L;
        gamma(i,:)=ifft(Xf,NFFT,'symmetric');

    end
    
    if opt==1
        chan=input('channel number: ');
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(6,1,1)
        plot(time,data(chan,1:fLen));title('Signal');ylabel('amplitude');%ylim([-30 30])
        subplot(6,1,2)
        plot(time,delta(chan,1:fLen));title('Delta');ylabel('amplitude');%ylim([-30 30])
        subplot(6,1,3)
        plot(time,theta(chan,1:fLen));title('Theta');ylabel('amplitude');%ylim([-30 30])
        subplot(6,1,4)
        plot(time,alpha(chan,1:fLen));title('Alpha');ylabel('amplitude');%ylim([-30 30])
        subplot(6,1,5)
        plot(time,beta(chan,1:fLen));title('Beta');ylabel('amplitude');%ylim([-30 30])
        subplot(6,1,6)
        plot(time,gamma(chan,1:fLen));title('Gamma');ylabel('amplitude');%ylim([-30 30])
    end
%%% ---------------------------- power calculation --------------------------------
    power=[sum(delta'.^2)'/length(X) sum(theta'.^2)'/length(X) sum(alpha'.^2)'/length(X) sum(beta'.^2)'/length(X) sum(gamma'.^2)'/length(X)];

    if opt==1
        figure;
        image(power)
        xticks([1:5]);
        yticks([1:channel]);
        ylabel('Channels')
        xticklabels({'Delta','Theta','Alpha','Beta','Gamma'})
        title('Power')
        colorbar
    end
