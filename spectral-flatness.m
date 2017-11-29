[y,fs]=audioread('bird.wav');  
x=periodogram(y);  % Power spectral density
num=geomean(x);
denomin=mean(x);
spf=num/denomin;
disp(spf);





   
