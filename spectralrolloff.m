[x,f_s]=audioread('bird.wav');
X=fft(x);
vsr=zeros(1,size(X,2));
Sum=sum(X,1);
kappa=0.85;
for(n=1:length(vsr))
    vsr(n)=find(cumsum(X(:,n)) >=kappa*Sum(n),1);
end
    
 vsr1 =vsr/size(X,1)*f_s/2;
 vrms=sqrt(mean(X.^2));
 vsr2=vsr1.*vrms;
 vsr=sum(vsr2)./sum(vrms);
 disp(vsr);
   
 
   
