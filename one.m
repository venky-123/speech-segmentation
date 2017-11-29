
[x , fs] = audioread('123Sam16K1.wav');
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
speaker = zeros(1,coln);
specCentroid = zeros(1,coln);
specSpread = zeros(1,coln);
specSkewness = zeros(1,coln);
specKurtosis = zeros(1,coln);
specEntropy = zeros(1,coln);
STE = zeros(1,coln);
zcr = zeros(1,coln);
unvoiced = zeros(1,coln);
pitch = zeros(1,coln);
specFM = zeros(1,coln);
% initialize the signal time segment index
indx = 0;
silence = zeros(1,coln);
% perform STFT
for col = 1:coln
    % windowing
    xw = x(indx+1:indx+wlen).*win;
 %   xw2 = x(indx+2:indx+wlen+1).*win ;
 
  % FFT
    X = fft(xw, nfft);
    Xmag = abs(X(1:rown));
    Xtot = sum(Xmag);
    
    S =  Xmag / abs((Xtot));
    %ste calculation
    STE(:,col) = sum(Xmag.^2)/lencol;
    
    % update the stft matrix
    stft(:, col) = X(1:rown);
    specAmp(:,col) = S(1:rown,1);

    %Spectral Flat measure Calculation
    
    specFM(:,col) = geomean(Xmag(1:rown)) * rown/Xtot ;  
    
    if(STE(:,col) > 0.5 && specFM(:,col) > 0.2)
        unvoiced(:,col) = 1;
    end
    
    if(STE(:,col) < 0.5)
        silence(:,col) = 1;
    end
    if(col<15)
        speaker(:,col) = 0;
    end
    if(col > 20)
        count = 0;
       for colz = col-20 : col
           if(silence(:,colz) == 1 )
               count = count+1;
           end
       end
           if(count >= 20)
               speaker(:,col) = 1;
           else
               speaker(:,col) = 0;
           end
          
    end

    %Pitch calculation
    
    %Filter of cutoff 900HZ
    [b , a] = butter(2,900/(fs/2));
    xclip = filter(b,a,xw);
    
    %30% of max value of y
    cl = 0.3 * max(xclip);
    
    %center clipping
    ii = 0;
    for ii = 1:wlen
        
     if (abs(xclip(ii)) > cl)
        
         if  xclip(ii) > 0
                 xclip(ii) = xclip(ii) - cl;
         else 
                 xclip(ii) = xclip(ii) + cl;
         end
     else 
         
         xclip(ii) = 0;
     end
        
    end
    
    
    xautocorr = autocorr(xclip,wlen-1);
    
    
    auto=xautocorr(21:240);
    max1=0;
    for uu=1:220
        if(auto(uu)>max1)
            max1=auto(uu);
             sample_no=uu;
        end
    end
    if(speaker(:,col)< 1 || unvoiced(:,col) < 1)
    if( max1 >= 0.3*xautocorr(1))
            pitch(:,col)=1/((20+sample_no)*(1/fs));
    else
            pitch(:,col) = 0;
    end
    end
    if(speaker(:,col) == 1 || unvoiced(:,col) == 1 )
        pitch(:,col) = 0;
    end
 
    
    %ZCR calculation
    zcr(:,col) = 1/2 * (sum(abs(sign(xw(1:wlen-1)) - sign(xw(2:wlen))))) * fs / wlen;
    if(speaker(:,col) < 1 )
    specCentroid(:,col) = sum ( fcol .* S(1:lencol) );
    
    %Spectral Spread
    specSpread(:,col) = sqrt(sum( ((fcol - specCentroid(1,col)).^2) .* S(1:lencol) ));
    
    %Spectral Skewness start
    specSkewness(:,col) = (sum( ((fcol - specCentroid(1,col)).^3) .* S(1:lencol) ))/specSpread(:,col)^3;
    
    %Spectral kurtosis
    specKurtosis(:,col) = (sum( ((fcol - specCentroid(1,col)).^4) .* S(1:lencol) ))/specSpread(:,col)^4;
    
    %removing the mean from skewness.
    specSkewness(:,col) = specSkewness(:,col) - (sum(specSkewness(:,col))/lencol);

    %Spectral entropy
    specEntropy(:,col) = -sum(S .* log2(S));
    end
    
    if(speaker(:,col) == 1 )
        specCentroid(:,col) = 0;
    
    %Spectral Spread
    specSpread(:,col) = 0;
    
    %Spectral Skewness start
    specSkewness(:,col) = 0;
    
    %Spectral kurtosis
    specKurtosis(:,col) = 0;
    
    %removing the mean from skewness.
    specSkewness(:,col) = 0 ;

    %Spectral entropy
    specEntropy(:,col) = 0;
    
    end
       
   
        % update the index
    indx = indx + hop;
end

%KL2 = zeros(1,coln);

wlen1 = 512*1;
hop1 = wlen1/4;
nfft1 = 2^nextpow2(wlen1);


rown1 = ceil((1+nfft1)/2);            % calculate the total number of rows
coln1 = 1+fix((xlen-wlen1)/hop1);

t1 = (wlen1/2:hop1:wlen1/2+(coln1-1)*hop1)/fs;

klwin = hamming(wlen1, 'periodic');
klft = zeros(rown1,coln1);
indxa = 0;

for cola = 1:coln1
    % windowing
    xwa = x(indxa+1:indxa+wlen1).*klwin;
    
    % FFT
    Xa = fft(xwa, nfft1);
    SpecAmp1 =  abs(Xa)/sum(abs(Xa));
    klft(:,cola) = SpecAmp1(1:rown1);
    
    
    % update the index
    indxa = indxa + hop1;
end

    kl2 = zeros(1,coln1);
     klcos = zeros(1,coln1);
     klnew = zeros(1,coln1);
for colb = 1 : coln1-1

    klab = sum (  klft(:,colb) .* (log2(klft(:,colb)) - log2(klft(:,colb+1)))  );
    klba = sum (  klft(:,colb+1) .* (log2(klft(:,colb+1)) - log2(klft(:,colb))) );
    klcos(1,colb) = sum(klft(:,colb) .* klft(:,colb+1));
    kl2(1,colb) = klab + klba;
    klnew(1,colb) = sum( ( klft(:,colb) .* log2(klft(:,colb))) .* ( klft(:,colb+1) .* log2(klft(:,colb+1))));
end


% NOW no of windows for which KL2 is calculated 

% now = 64;
% 
% 
% 
% for i = 1:now/2:coln-now 
%     
%     klaa = zeros(rown*now/2,1);
%     klbb = zeros(rown*now/2,1);
%     for j = 1:now/2
%         
%         klaa((j-1)*rown + 1: j*rown) = specAmp(:,i+j-1);
%         klbb((j-1)*rown + 1: j*rown) = specAmp(:,i+(now/2)+j-1);
%         
%     end
%     
%     
%         
%         kla = sum(klaa .* (log2(klaa) - log2(klbb)));
%     
%         klb = sum(klaa .* ( log2(klaa) - log2(klbb)));
%     KL2(:,i) = kla + klb;
% end



 
figure(3);
subplot(411);
plot(t,specCentroid);
xlabel('time');  ylabel('spectral centroid');
subplot(412);
plot(t,specSpread);
xlabel('time'); ylabel('spectral Spread');
subplot(413);
plot(t,specSkewness);
xlabel('time'); ylabel('Spectral Skewness');
subplot(414);
plot(treal,x);
xlabel('time'); ylabel('Signal');


figure(5);
subplot(311);
plot(treal,x);
xlabel('time'); ylabel('Signal');
subplot(312);
plot(t,STE);
xlabel('time'); ylabel('Short term Energy');
subplot(313);
plot(t,specEntropy);
xlabel('time'); ylabel('SPectral Entropy');

figure(6);
subplot(411);
plot(treal,x);
xlabel('time'); ylabel('Signal');
subplot(412);
plot(t1 , kl2);
xlabel('time'); ylabel('KL distance');
title('KL2 distance');
subplot(413);
plot(t1,klcos);
subplot(414);
plot(t1,klnew);



figure(7);
subplot(211);
plot(treal,x);
xlabel('time'); ylabel('Signal');
subplot(212);
plot(t,zcr);
ylabel('Zcr');
xlabel('time');
title('Graph between ZCR and time');

figure(8);
ftset = [specCentroid; specSkewness; specSpread; specKurtosis; STE; specEntropy; kl2; ];
mesh(ftset);

ft = transpose(ftset);
[IDX, Ck] = kmeans(ft,3);

idt = transpose(IDX);

plot(t,idt,treal,x);


figure(9)
subplot(211)
plot(t,pitch);
title('pitch');
subplot(212);
title('Signal');
plot(treal,x);

figure(10)
subplot(211)
plot(t,specFM);
title('SPectral Flat measure');
subplot(212);
plot(treal,x);

figure(11)
subplot(211);
plot(t,unvoiced);
title('unvoiced');
subplot(212);
plot(treal,x);

figure(12)
subplot(211);
plot(t,silence);
title('silence');
subplot(212);
plot(treal,x);

figure(13)
subplot(211);
plot(t,speaker);
title('speaker tranns');
subplot(212);
plot(treal,x,t,speaker);

figure(18);
ftset = [specCentroid; specSkewness; specSpread; specKurtosis; STE; specEntropy; kl2; ];
mesh(ftset);

