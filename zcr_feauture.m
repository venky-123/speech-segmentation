[x , fs] = audioread('wav/noise.wav');
% function: [stft, f, t] = stft(x, wlen, hop, nfft, fs)


wlen = 512;
hop = wlen/4;
nfft = 2^nextpow2(wlen);

% represent x as column-vector
x = x(:);

% length of the signal
xlen = length(x);

% form a periodic hamming window
win = hamming(wlen, 'periodic');

% stft matrix estimation and preallocation
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/hop);      % calculate the total number of columns
stft = zeros(rown, coln);           % form the stft matrix
specAmp = zeros(rown , coln);

% calculate the time and frequency vectors
t = (wlen/2:hop:wlen/2+(coln-1)*hop)/fs;
treal = 1/fs:1/fs:xlen/fs;
f = (0:rown-1)*fs/nfft;
fcol = transpose(f);
lencol = length(fcol);

specCentroid = zeros(1,coln);
specSpread = zeros(1,coln);
specSkewness = zeros(1,coln);
specKurtosis = zeros(1,coln);
specEntropy = zeros(1,coln);
STE = zeros(1,coln);
zcr = zeros(1,coln);

% initialize the signal time segment index
indx = 0;

% perform STFT
for col = 1:coln
    % windowing
    xw = x(indx+1:indx+wlen).*win;
    xw2 = x(indx+2:indx+wlen+1).*win ;
    
    % FFT
    X = fft(xw, nfft);
    Xmag = abs(X);
    S =  Xmag / abs((sum(Xmag)));
    %ste calculation
    STE(:,col) = sum(Xmag.^2)/lencol;
    
    % update the stft matrix
    stft(:, col) = X(1:rown);
    specAmp(:,col) = S(1:rown,1);
    
    specCentroid(:,col) = sum ( fcol .* S(1:lencol) );
    
    %Spectral Spread
    specSpread(:,col) = sqrt(sum( ((fcol - specCentroid(1,col)).^2) .* S(1:lencol) ));
    
    %Spectral Skewness start
    specSkewness(:,col) = (sum( ((fcol - specCentroid(1,col)).^3) .* S(1:lencol) ))/specSpread(:,col)^3;
    
    %removing the mean from skewness.
    specSkewness(:,col) = specSkewness(:,col) - (sum(specSkewness(:,col))/lencol);

    %Spectral entropy
    specEntropy(:,col) = -sum(S .* log2(S));
    zcr(:,col) = 1/2 * (sum(abs(sign(xw) - sign(xw2)))) * fs / wlen;
    % update the index
    indx = indx + hop;
end

KL2 = zeros(1,coln);


% NOW no of windows for which KL2 is calculated 

now = 64;



for i = 1:now/2:coln-now 
    
    klaa = zeros(rown*now/2,1);
    klbb = zeros(rown*now/2,1);
    for j = 1:now/2
        
        klaa((j-1)*rown + 1: j*rown) = specAmp(:,i+j-1);
        klbb((j-1)*rown + 1: j*rown) = specAmp(:,i+(now/2)+j-1);
        
    end
    
    
        
        kla = sum(klaa .* (log2(klaa) - log2(klbb)));
    
        klb = sum(klaa .* ( log2(klaa) - log2(klbb)));
    KL2(:,i) = kla + klb;
end



 
subplot(211);
plot(abs(stft));
subplot(212);
plot(specAmp);
figure(2);
mesh(specAmp);
figure(3);
subplot(411);
plot(t,specCentroid);
subplot(412);
plot(t,specSpread);
subplot(413);
plot(t,specSkewness);
subplot(414);
plot(treal,x);

figure(5);
subplot(311);
plot(treal,x);
subplot(312);
plot(t,STE);
subplot(313);
plot(t,specEntropy);

figure(6);
subplot(211);
plot(treal,x);
subplot(212);
plot(t,KL2);


figure(7);
subplot(211);
plot(x);
subplot(212);
plot(t,zcr);
ylabel('Zcr');
xlabel('time');
title('Graph between ZCR and time');