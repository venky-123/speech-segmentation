[y,f_s]=audioread('foot.wav');
X=fft(y);
 


UseBookDefinition = true;
 
    if (UseBookDefinition)
        % compute mean and standard deviation
        mu_x    = mean(abs(X), 1);
        std_x   = std(abs(X), 1);
 
        % remove mean
        X       = X - repmat(mu_x, size(X,1), 1);
 
        % compute kurtosis
        vsk1     = sum ((X.^4)./(repmat(std_x, size(X,1), 1).^4*size(X,1)));
    else
        % interpret the spectrum as pdf, not as signal
        f       = linspace(0, f_s/2, size(X,1));
        % compute mean and standard deviation
        mu_X    = (f * X) ./ sum(X,1);
        tmp     = repmat(f, size(X,2),1) - repmat(mu_X, size(X,1),1)';
        var_X   = diag (tmp.^2 * X) ./ (sum(X,1)'*size(X,1));
 
        vsk1    = diag (tmp.^4 * X) ./ (var_X.^2 .* sum(X,1)'*size(X,1));
        
    end
    vsk1     = vsk1-3;
    vrms =  sqrt(mean(y.^2));
    vsk2=vsk1.*vrms;
    vsk=sum(vsk2)./sum(vrms);
    disp(vsk);
    
