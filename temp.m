        N = 50;
        Fs = 16000;
       x = sin(2*pi*110*(0:(1/Fs):5));
       [Pxx,F] = periodogram(x,[],512,Fs);
       
       [fil , melfreq, iFilter] = melfilters(8,F);
       iFilter = iFilter * Fs/512 ;