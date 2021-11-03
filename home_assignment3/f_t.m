function ret = f_t(lambda, t, tau) % conditional distribution for t. 
    d = length(lambda);
    tdiff = t(2:end) - t(1:end-1);
    ni = zeros(d,1);
    for k=1:d
       ni(k)=n_i(tau,t(k),t(k+1)); %counts occurences. 
    end
    ret = exp(sum(log(lambda).*ni + +log(tdiff) - lambda.*tdiff));
end