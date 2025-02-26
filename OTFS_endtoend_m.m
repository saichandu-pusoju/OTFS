clc;
clear all ;
close all;

EbdB = -10:2:10;
Eb = 10.^(EbdB/10);
No = 1;
SNR = 2*Eb/No;
SNRdB = 10*log10(SNR);
M=32;
N=16;
nTaps=5;
DelayTaps =0:4;
DopplerTaps = [0 1 2 3 4];
Ncp = max(DelayTaps);
BER_MMSE = zeros(length(Eb),1);
BER_ZF = zeros(length(Eb),1);
ITER = 100;
F_M = 1/sqrt(M)*dftmtx(M);
F_N= 1/sqrt(N)*dftmtx(N);


for ite=1:ITER
    ite
    XddBits = randi([0,1],M,N);
    h = sqrt(1/2)*(randn(1,nTaps)+1j*randn(1,nTaps));
    Hmat = zeros(M*N,M*N);
    omega = exp(1j*2*pi/(M*N));
    for tx = 1:nTaps
        Hmat =Hmat + h(tx)*circshift(eye(M*N),DelayTaps(tx))*...
        (diag(omega.^((0:M*N-1)*DopplerTaps(tx))));
    end
    Heff = kron(F_N,eye(M))*Hmat*kron(F_N',eye(M));
    ChNoise = sqrt(No/2)*(randn(1,M*N)+1j*randn(1,M*N));


    for ix= 1:length(Eb)
        X_DD =sqrt(Eb(ix))*(2*XddBits-1);
        S_mat = eye(M)*X_DD*F_N';
        TxSamples = reshape(S_mat,M*N,1).';
        TxSamplesCP = [TxSamples(M*N-Ncp+1:M*N) TxSamples];
        RxsamplesCP= 0;

        for tx=1:nTaps
            Doppler = exp(1j*2*pi/M*(-Ncp:M*N-1)*DopplerTaps(tx)/N);
            RxsamplesCP = RxsamplesCP+h(tx)*circshift(TxSamplesCP.*Doppler,[1,DelayTaps(tx)]);
        end
        Rxsamples = RxsamplesCP(Ncp+1:M*N+Ncp)+ChNoise;
        R_mat = reshape(Rxsamples.',M,N);
        Y_DD = eye(M)*R_mat*F_N;
        y_DD = reshape(Y_DD,M*N,1);
        
        xhatMMSE = inv(Heff'*Heff + eye(M*N)/Eb(ix))*Heff'*y_DD;
        DecodedBits = (real(xhatMMSE)>=0);
        BER_MMSE(ix)= BER_MMSE(ix)+sum(DecodedBits~=reshape(XddBits,M*N,1));

        xhatZF = pinv(Heff)*y_DD;
        DecodedBits = (real(xhatZF)>=0);
        BER_ZF(ix)=BER_ZF(ix)+ sum(DecodedBits~=reshape(XddBits,M*N,1));
    end
end

BER_MMSE = BER_MMSE/M/N/ITER;
BER_ZF = BER_ZF/M/N/ITER;
semilogy(EbdB,BER_MMSE,'b-s','lineWidth',3.0,'MarkerFaceColor','b','MarkerSize',9.0);
hold on;
grid on;
semilogy(EbdB,BER_ZF,'r-o','lineWidth',3.0,'MarkerFaceColor','b','MarkerSize',9.0);
axis tight;
legend('MMSE','ZF');
title('OTFS BER v/s SNR');
xlabel('SNR(dB)');
ylabel('BER');