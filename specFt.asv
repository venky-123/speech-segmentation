function [stft,specCentroid, specSpread , specKurtosis , specSkewness, specFM , subSpecCentroid , subMelCentroid , specEntropy, f, t, treal] = specFt(x,wlen, hop , nfft, fs)



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
pitch = zeros(1,coln);
specFM = zeros(1,coln);
specFlux = zeros(1,coln);
unvoiced = zeros(1,coln);
speaker = zeros(1,coln);
silence = zeros(1,coln);
sil_inv = zeros(1,coln);

%Subband spectral centroids

subSpecCentroid = zeros(8,coln); 
subMelCentroid = zeros(8,coln); 

[mFilters  ,~ , fScale ] = melfilters(12,fcol); 
% mFilters = mFilters/sum(mFilters);

%

% initialize the signal time segment index
indx = 0;

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
   
    sil_inv(:,col) = 1 - silence(:,col);
    
    specCentroid(:,col) = sum ( fcol .* S(1:lencol) )   * sil_inv(col);
    
    %Spectral Spread
    specSpread(:,col) = (sqrt(sum( ((fcol - specCentroid(1,col)).^2) .* S(1:lencol) )));
    
    %Spectral Skewness start
    specSkewness(:,col) = ((sum( ((fcol - specCentroid(1,col)).^3) .* S(1:lencol)))/specSpread(:,col)^3 );
    
    %Spectral kurtosis
    specKurtosis(:,col) = ((sum( ((fcol - specCentroid(1,col)).^4) .* S(1:lencol) ))/specSpread(:,col)^4);
    
    %removing the mean from skewness..*
    specSkewness(:,col) = (specSkewness(:,col) - (sum(specSkewness(:,col))/lencol));

    %Spectral entropy
    specEntropy(:,col) = (-sum(S .* log2(S))).*  sil_inv(:,col);
    
   
    %ZCR calculation
    zcr(:,col) = 1/2 * (sum(abs(sign(xw(1:wlen-1)) - sign(xw(2:wlen))))) * fs / wlen;
    
    
    if (col > 1)
    
        specFlux(:,col) = 1 - ((sum( Xmag .* temp1))/(sqrt((sum(Xmag.^2)) * (sum(temp1.^2))  )));
    
    end
    temp1 = Xmag;
    
    
    
    
    %Subband Spectral centroids.
    
    for xtem = 1:32:rown-32
            
        
        subS = (Xmag(xtem:xtem+32)/sum(Xmag(xtem:xtem+32)) ) ;
        
        subSpecCentroid(floor(xtem/32) + 1 ,col) = sum ( fcol(xtem:xtem+32) .* subS);
            
    end
    
    
    for xx = 1:12
        
        subS = Xmag(fScale(xx,1) : fScale(xx,2)) .* transpose(mFilters(xx,fScale(xx,1) : fScale(xx,2))) /sum(Xmag(fScale(xx,1) : fScale(xx,2)));
        subMelCentroid(xx,col) = sum( subS .* fcol (fScale(xx,1) : fScale(xx,2)) );
        
        
        
        
    end
    
    %
    % update the index
    indx = indx + hop;
    
    
    
    
    
    
    
    
end




% calculate the time and frequency vectors
t = (wlen/2:hop:wlen/2+(coln-1)*hop)/fs;
f = (0:rown-1)*fs/nfft;
treal = 1/fs:1/fs:xlen/fs;
end