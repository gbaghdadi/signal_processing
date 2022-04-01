function [delta,theta,alpha,beta,gamma,power]=func_FE_BrainWavesEnergy_WT(data,Fs,Wname,opt)

%%% This function extracts five brain wave signals in time domain and also 
%%% the power of each band usning Wavelet transform algorithm

%%% input: data (channels,samples)
%%%        Fs  frequecy rate (Hz)
%%%        Wname   'amor' 'bump' 'Morse'
%%%        opt   option for plotting the results 0=no plot 

%%% output: delta -> delta band in time domain (channels,samples)
%%%         theta -> theta band in time domain (channels,samples)
%%%         alpha -> alpha band in time domain (channels,samples)
%%%         beta -> beta band in time domain (channels,samples)
%%%         gamma -> gamma band in time domain (channels,samples)
%%%         power -> brain waves power (channels,power of each band [delta, theta, alpha, beta, gamma])

%%% writen by Golnaz Baghdadi 7/25/2021
%%---------------------------------------------------------------------------------------------------------

[channel,N]=size(data);
time=linspace(0,fix(N/Fs),N);
lenTime=length(time);

for ch=1:channel
    
    s=data(ch,:);
    [wt,f] = cwt(s, Fs,Wname); 

    gammaC = icwt(wt,f,[30 max(f)]);
    betaC = icwt(wt,f,[13 30]);
    alphaC = icwt(wt,f,[8 13]);
    thetaC = icwt(wt,f,[4 8]);
    deltaC = icwt(wt,f,[min(f) 4]);
    
    delta(ch,:)=deltaC;
    theta(ch,:)=thetaC;
    alpha(ch,:)=alphaC;
    beta(ch,:)=betaC;
    gamma(ch,:)=gammaC;


end
power=[sum(delta'.^2)'/N sum(theta'.^2)'/N sum(alpha'.^2)'/N sum(beta'.^2)'/N sum(gamma'.^2)'/N];


if opt==1
    I=func_PL_fftSpect(gamma,Fs,opt);
    fprintf('Gamma:Maximum occurs at %3.2f Hz.\n',I);

    I=func_PL_fftSpect(beta,Fs,opt);
    fprintf('Gamma:Maximum occurs at %3.2f Hz.\n',I);

    I=func_PL_fftSpect(alpha,Fs,opt);
    fprintf('Gamma:Maximum occurs at %3.2f Hz.\n',I);

    I=func_PL_fftSpect(theta,Fs,opt);
    fprintf('Gamma:Maximum occurs at %3.2f Hz.\n',I);

    I=func_PL_fftSpect(delta,Fs,opt);
    fprintf('Gamma:Maximum occurs at %3.2f Hz.\n',I);
    
    figure;
    plot(s,'r')
    hold on
    subplot(5,1,1); plot(1:1:length(gamma),gamma);title('GAMMA');
    subplot(5,1,2); plot(1:1:length(beta), beta); title('BETA');
    subplot(5,1,3); plot(1:1:length(alpha),alpha); title('ALPHA'); 
    subplot(5,1,4); plot(1:1:length(theta),theta);title('THETA');
    subplot(5,1,5);plot(1:1:length(delta),delta);title('DELTA');

    chan=input('channel number: ');
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(6,1,1)
    plot(time,data(chan,1:lenTime));title('Signal');ylabel('amplitude');%ylim([-30 30])
    subplot(6,1,2)
    plot(time,delta(chan,1:lenTime));title('Delta');ylabel('amplitude');%ylim([-30 30])
    subplot(6,1,3)
    plot(time,theta(chan,1:lenTime));title('Theta');ylabel('amplitude');%ylim([-30 30])
    subplot(6,1,4)
    plot(time,alpha(chan,1:lenTime));title('Alpha');ylabel('amplitude');%ylim([-30 30])
    subplot(6,1,5)
    plot(time,beta(chan,1:lenTime));title('Beta');ylabel('amplitude');%ylim([-30 30])
    subplot(6,1,6)
    plot(time,gamma(chan,1:lenTime));title('Gamma');ylabel('amplitude');%ylim([-30 30])

end
