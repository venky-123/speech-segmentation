function    [pitch t] =   calculate_pitch(x,wlen,hop,nfft,fs)

win=hamming(wlen,'periodic');
xlen=length(x);
rown=ceil((1+nfft)/2);
coln=1+fix((xlen-wlen)/hop);
pitch=zeros(1,coln);
time=xlen/(fs*coln);

for j=1:n
     xw = x(indx+1:indx+wlen).*win;
     t=[0:length(xw)-1]/fs;
     y=fft(xw);
     ms1=fs/1000;
     ms20=fs/50;
     Y=fft(xw.*hamming(length(xw)));
     C=fft(log(abs(Y)+eps));
     q=(ms1:ms20)/fs;
     [Maxamp_at_pitch,fx]=max(abs(C(ms1:ms20)));
     pitch(j) = fs/(ms1+fx-1);

end

