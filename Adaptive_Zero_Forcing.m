close all;
clear all;
clc;
train_bit=100; % number of bits to be transmitted
EN=[0:30]; % E/N values
x=1;
while x~=length(EN)+1
            t=2;
            r=2;
            % Transmitter
            bitgentb = randi([0 1],1,train_bit); % Random generation of 0 and 1
            s = 2*bitgentb-1; % BPSK 0 -> -1  &   1 -> 0
            sModtb = kron(s,ones(r,1)); 
            sModtb = reshape(sModtb,[r,t,train_bit/t]); 
            ngtb = 1/sqrt(2)*[randn(r,train_bit/t) + j*randn(r,train_bit/t)]; % white gaussian noise
            hrtb = 1/sqrt(2)*[randn(r,t,train_bit/t) + j*randn(r,t,train_bit/t)]; % Rayleigh channel
            nrtb = 1/sqrt(2)*[randn(r,train_bit/t) + j*randn(r,train_bit/t)]; % white gaussian noise
            stb=sqrt(1/2); %Non-Centrality Parameter
            sigmatb=1/sqrt(2*2);
            hritb=((sigmatb*randn(r,t,train_bit/t)+stb)+1i*(randn(r,t,train_bit/t)*sigmatb+0)); %Rician Fading
            nritb = 1/sqrt(2)*[randn(r,train_bit/t) + j*randn(r,train_bit/t)]; % white gaussian noise
            % Addition of channel and noise
            ygtb = squeeze(sum(sModtb,2)) + 10^(-EN(x)/20)*ngtb; % Gaussian
            yrtb = squeeze(sum(hrtb.*sModtb,2)) + 10^(-EN(x)/20)*nrtb; % Rayleigh
            yritb = squeeze(sum(hritb.*sModtb,2)) + 10^(-EN(x)/20)*nritb; % Rician
            % Receiver in Gaussian Channel
            yModgtb = kron(ygtb,ones(1,2)); % formatting the received symbol for equalization
            yModgtb = sum(yModgtb,1); % H^H * y 
            yModgtb = kron(reshape(yModgtb,2,train_bit/t),ones(1,2)); % formatting
            
            hCofrtb = zeros(2,2,train_bit/t)  ; 
            hCofrtb(1,1,:) = sum(hrtb(:,2,:).*conj(hrtb(:,2,:)),1);  % d term
            hCofrtb(2,2,:) = sum(hrtb(:,1,:).*conj(hrtb(:,1,:)),1);  % a term
            hCofrtb(2,1,:) = -sum(hrtb(:,2,:).*conj(hrtb(:,1,:)),1); % c term
            hCofrtb(1,2,:) = -sum(hrtb(:,1,:).*conj(hrtb(:,2,:)),1); % b term
            hDenrtb = ((hCofrtb(1,1,:).*hCofrtb(2,2,:)) - (hCofrtb(1,2,:).*hCofrtb(2,1,:))); % ad-bc term
            hDenrtb = reshape(kron(reshape(hDenrtb,1,train_bit/t),ones(2,2)),2,2,train_bit/t);  % formatting for division
            hInvrtb = hCofrtb./hDenrtb; % inv(H^H*H)
            hModrtb =  reshape(conj(hrtb),r,train_bit); % H^H operation
            yModrtb = kron(yrtb,ones(1,2)); % formatting the received symbol for equalization
            yModrtb = sum(hModrtb.*yModrtb,1); % H^H * y 
            yModrtb = kron(reshape(yModrtb,2,train_bit/t),ones(1,2)); % formatting
            yHatrtb = sum(reshape(hInvrtb,2,train_bit).*yModrtb,1); % inv(H^H*H)*H^H*y
            
            
            hCofritb = zeros(2,2,train_bit/t)  ; 
            hCofritb(1,1,:) = sum(hritb(:,2,:).*conj(hritb(:,2,:)),1);  % d term
            hCofritb(2,2,:) = sum(hritb(:,1,:).*conj(hritb(:,1,:)),1);  % a term
            hCofritb(2,1,:) = -sum(hritb(:,2,:).*conj(hritb(:,1,:)),1); % c term
            hCofritb(1,2,:) = -sum(hritb(:,1,:).*conj(hritb(:,2,:)),1); % b term
            hDenritb = ((hCofritb(1,1,:).*hCofritb(2,2,:)) - (hCofritb(1,2,:).*hCofritb(2,1,:))); % ad-bc term
            hDenritb = reshape(kron(reshape(hDenritb,1,train_bit/t),ones(2,2)),2,2,train_bit/t);  % formatting for division
            hInvritb = hCofritb./hDenritb; % inv(H^H*H)
            hModritb =  reshape(conj(hritb),r,train_bit); % H^H operation
            yModritb = kron(yritb,ones(1,2)); % formatting the received symbol for equalization
            yModritb = sum(hModritb.*yModritb,1); % H^H * y 
            yModritb = kron(reshape(yModritb,2,train_bit/t),ones(1,2)); % formatting
            yHatritb = sum(reshape(hInvritb,2,train_bit).*yModritb,1); % inv(H^H*H)*H^H*y
            
            x=x+1;
end

ag=zeros(2,2);
ar=zeros(2,2);
ari=zeros(2,2);
for i=1:50
ar=ar+hInvrtb(:,:,i);
ari=ari+hInvritb(:,:,i);
end
ar=ar./50;
ari=ari./50;


y=1;        
inp=10000;
while y~=length(EN)+1
            t=2;
            r=2;
            % Transmitter
            bitgen = randi([0 1],1,inp); % Random generation of 0 and 1
            s1 = 2*bitgen-1; % BPSK 0 -> -1  &   1 -> 0
            sMod = kron(s1,ones(r,1)); 
            sMod = reshape(sMod,[r,t,inp/t]); 
            
            rc = 1/sqrt(2)*[randn(r,t,inp/t) + j*randn(r,t,inp/t)]; % Rayleigh channel
            hra=zeros(2,2);
            for i=1:5000
                hra(:,:,i)=rc(:,:,i).*ar;
            end
            nr = 1/sqrt(2)*[randn(r,inp/t) + j*randn(r,inp/t)]; % white gaussian noise
            hr = 1/sqrt(2)*[randn(r,t,inp/t) + j*randn(r,t,inp/t)]; % Rayleigh channel
            
            s=sqrt(1/2); %Non-Centrality Parameter
            sigma=1/sqrt(2*2);
            ric=((sigma*randn(r,t,inp/t)+s)+1i*(randn(r,t,inp/t)*sigma+0)); %Rician Fading
            hria=zeros(2,2);
            for i=1:5000
                hria(:,:,i)=ric(:,:,i).*ari;
            end
            nri = 1/sqrt(2)*[randn(r,inp/t) + 1j*randn(r,inp/t)]; % white gaussian noise
            hri=((sigma*randn(r,t,inp/t)+s)+1i*(randn(r,t,inp/t)*sigma+0)); %Rician Fading
            
            % Addition of channel and noise
            yra = squeeze(sum(hra.*sMod,2)) + 10^(-EN(y)/20)*nr; % Rayleigh
            yr = squeeze(sum(hr.*sMod,2)) + 10^(-EN(y)/20)*nr; % Rayleigh
            
            yria = squeeze(sum(hria.*sMod,2)) + 10^(-EN(y)/20)*nri; % Rician
            yri = squeeze(sum(hri.*sMod,2)) + 10^(-EN(y)/20)*nri; % Rician
            % Receiver in Gaussian Channel
            yHatga = awgn(s1,y,'measured'); % inv(H^H*H)*H^H*y
            % receiver - hard decision decoding
            bitrecga = real(yHatga)>0;
            yHatg = awgn(s1,y,'measured'); % inv(H^H*H)*H^H*y
            % receiver - hard decision decoding
            bitrecg = real(yHatg)>0;
            
            
            
            hCofra = zeros(2,2,inp/t)  ; 
            hCofra(1,1,:) = sum(hra(:,2,:).*conj(hra(:,2,:)),1);  % d term
            hCofra(2,2,:) = sum(hra(:,1,:).*conj(hra(:,1,:)),1);  % a term
            hCofra(2,1,:) = -sum(hra(:,2,:).*conj(hra(:,1,:)),1); % c term
            hCofra(1,2,:) = -sum(hra(:,1,:).*conj(hra(:,2,:)),1); % b term
            hDenra = ((hCofra(1,1,:).*hCofra(2,2,:)) - (hCofra(1,2,:).*hCofra(2,1,:))); % ad-bc term
            hDenra = reshape(kron(reshape(hDenra,1,inp/t),ones(2,2)),2,2,inp/t);  % formatting for division
            hInvra = hCofra./hDenra; % inv(H^H*H)
            hModra =  reshape(conj(hra),r,inp); % H^H operation
            yModra = kron(yra,ones(1,2)); % formatting the received symbol for equalization
            yModra = sum(hModra.*yModra,1); % H^H * y 
            yModra = kron(reshape(yModra,2,inp/t),ones(1,2)); % formatting
            yHatra = sum(reshape(hInvra,2,inp).*yModra,1); % inv(H^H*H)*H^H*y
            % receiver - hard decision decoding
            bitrecra = real(yHatra)>0;
            
            hCofr = zeros(2,2,inp/t)  ; 
            hCofr(1,1,:) = sum(hr(:,2,:).*conj(hr(:,2,:)),1);  % d term
            hCofr(2,2,:) = sum(hr(:,1,:).*conj(hr(:,1,:)),1);  % a term
            hCofr(2,1,:) = -sum(hr(:,2,:).*conj(hr(:,1,:)),1); % c term
            hCofr(1,2,:) = -sum(hr(:,1,:).*conj(hr(:,2,:)),1); % b term
            hDenr = ((hCofr(1,1,:).*hCofr(2,2,:)) - (hCofr(1,2,:).*hCofr(2,1,:))); % ad-bc term
            hDenr = reshape(kron(reshape(hDenr,1,inp/t),ones(2,2)),2,2,inp/t);  % formatting for division
            hInvr = hCofr./hDenr; % inv(H^H*H)
            hModr =  reshape(conj(hr),r,inp); % H^H operation
            yModr = kron(yr,ones(1,2)); % formatting the received symbol for equalization
            yModr = sum(hModr.*yModr,1); % H^H * y 
            yModr = kron(reshape(yModr,2,inp/t),ones(1,2)); % formatting
            yHatr = sum(reshape(hInvr,2,inp).*yModr,1); % inv(H^H*H)*H^H*y
            % receiver - hard decision decoding
            bitrecr = real(yHatr)>0;
            
            
            hCofria = zeros(2,2,inp/t)  ; 
            hCofria(1,1,:) = sum(hria(:,2,:).*conj(hria(:,2,:)),1);  % d term
            hCofria(2,2,:) = sum(hria(:,1,:).*conj(hria(:,1,:)),1);  % a term
            hCofria(2,1,:) = -sum(hria(:,2,:).*conj(hria(:,1,:)),1); % c term
            hCofria(1,2,:) = -sum(hria(:,1,:).*conj(hria(:,2,:)),1); % b term
            hDenria = ((hCofria(1,1,:).*hCofria(2,2,:)) - (hCofria(1,2,:).*hCofria(2,1,:))); % ad-bc term
            hDenria = reshape(kron(reshape(hDenria,1,inp/t),ones(2,2)),2,2,inp/t);  % formatting for division
            hInvria = hCofria./hDenria; % inv(H^H*H)
            hModria =  reshape(conj(hria),r,inp); % H^H operation
            yModria = kron(yria,ones(1,2)); % formatting the received symbol for equalization
            yModria = sum(hModria.*yModria,1); % H^H * y 
            yModria = kron(reshape(yModria,2,inp/t),ones(1,2)); % formatting
            yHatria = sum(reshape(hInvria,2,inp).*yModria,1); % inv(H^H*H)*H^H*y
            % Receiver - hard decision decoding
            bitrecria = real(yHatria)>0;
            
            hCofri = zeros(2,2,inp/t)  ; 
            hCofri(1,1,:) = sum(hri(:,2,:).*conj(hri(:,2,:)),1);  % d term
            hCofri(2,2,:) = sum(hri(:,1,:).*conj(hri(:,1,:)),1);  % a term
            hCofri(2,1,:) = -sum(hri(:,2,:).*conj(hri(:,1,:)),1); % c term
            hCofri(1,2,:) = -sum(hri(:,1,:).*conj(hri(:,2,:)),1); % b term
            hDenri = ((hCofri(1,1,:).*hCofri(2,2,:)) - (hCofri(1,2,:).*hCofri(2,1,:))); % ad-bc term
            hDenri = reshape(kron(reshape(hDenri,1,inp/t),ones(2,2)),2,2,inp/t);  % formatting for division
            hInvri = hCofri./hDenri; % inv(H^H*H)
            hModri =  reshape(conj(hri),r,inp); % H^H operation
            yModri = kron(yri,ones(1,2)); % formatting the received symbol for equalization
            yModri = sum(hModri.*yModri,1); % H^H * y 
            yModri = kron(reshape(yModri,2,inp/t),ones(1,2)); % formatting
            yHatri = sum(reshape(hInvri,2,inp).*yModri,1); % inv(H^H*H)*H^H*y
            % Receiver - hard decision decoding
            bitrecri = real(yHatri)>0;
            
            % Counting errors
            nErrga(y) = size(find([bitgen- bitrecga]),2);
            nErrg(y) = size(find([bitgen- bitrecg]),2);
            nErrra(y) = size(find([bitgen- bitrecra]),2);
            nErrr(y) = size(find([bitgen- bitrecr]),2);
            nErrria(y) = size(find([bitgen- bitrecria]),2);
            nErrri(y) = size(find([bitgen- bitrecri]),2);
            y=y+1;
end
simBerga = nErrga/inp; % Simulated ber of Gaussian
simBerg = nErrg/inp; % Simulated ber of Gaussian
simBerra = nErrra/inp; % Simulated ber of Rayleigh
simBerr = nErrr/inp; % Simulated ber of Rayleigh
simBerria = nErrria/inp; % Simulated ber of Rician
simBerri = nErrri/inp; % Simulated ber of Rician

figure
semilogy(EN,simBerg,'b*-','LineWidth',2);
hold on
semilogy(EN,simBerr,'r*-','LineWidth',2);
hold on
semilogy(EN,simBerri,'g*-','LineWidth',2);
axis([0 30 10^-5 0.5])
grid on
legend('Simulated(t=2, r=2, ZFGaussian)', 'Simulated(t=1,r=2, ZFRayleigh)','Simulated(t=1,r=2,ZFRician)');
xlabel('Average E/N(in dB)');
ylabel('Bit Error Rate');
title('BER for BPSK modulation having 2x2 MIMO & Zero Forcing equalizer')

figure
semilogy(EN,simBerga,'b*-','LineWidth',2);
hold on
semilogy(EN,simBerra,'r*-','LineWidth',2);
hold on
semilogy(EN,simBerria,'g*-','LineWidth',2);
axis([0 30 10^-5 0.5])
grid on
legend('Simulated(t=2, r=2, AZFGaussian)', 'Simulated(t=1,r=2, AZFRayleigh)','Simulated(t=1,r=2,AZFRician)');
xlabel('Average E/N(in dB)');
ylabel('Bit Error Rate');
title('BER for BPSK modulation having 2x2 MIMO & Adaptive Zero Forcing equalizer')


figure
subplot(3,1,1)
a=find(abs(real(yHatra))<=1);
b=real(yHatra(a));
plot(b(1:4));
hold on
for n=2:1:floor((length(b)-4)/4)
  plot(b((n-1)*4+1:n*4));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM AFTER AZF EQUALIZER in Rayleigh Channel')

subplot(3,1,2)
a=find(abs(real(yHatr))<=1);
b=real(yHatr(a));
plot(b(1:4));
hold on
for n=2:1:floor((length(b)-4)/4)
  plot(b((n-1)*4+1:n*4));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM AFTER ZF EQUALIZER in Rayleigh Channel')


subplot(3,1,3)
a=reshape(yr,1,10000);
b=find(abs(real(a))<=1);
c=real(yr(b));
plot(c(1:4));
hold on
for n=2:1:floor((length(c)-4)/4)
    plot(c((n-1)*4+1:n*4));
end
hold off
xlabel('Time(sec)')
title('EYE DIAGRAM WITHOUT EQUALIZER in Rayleigh Channel')


figure
subplot(3,1,1)
a=find(abs(real(yHatria))<=1);
b=real(yHatria(a));
plot(b(1:4));
hold on
for n=2:1:floor((length(b)-4)/4)
  plot(b((n-1)*4+1:n*4));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM AFTER AZF EQUALIZER in Rician Channel')

subplot(3,1,2)
a=find(abs(real(yHatri))<=1);
b=real(yHatri(a));
plot(b(1:4));
hold on
for n=2:1:floor((length(b)-4)/4)
  plot(b((n-1)*4+1:n*4));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM AFTER ZF EQUALIZER in Rician Channel')


subplot(3,1,3)
a=reshape(yri,1,10000);
b=find(abs(real(a))<=1);
c=real(yri(b));
plot(c(1:4));
hold on
for n=2:1:floor((length(c)-4)/4)
    plot(c((n-1)*4+1:n*4));
end
hold off
xlabel('Time(sec)')
title('EYE DIAGRAM WITHOUT EQUALIZER in Rician Channel')