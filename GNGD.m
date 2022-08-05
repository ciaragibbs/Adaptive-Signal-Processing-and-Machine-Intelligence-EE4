function [xOut,w,err] = GNGD(noise,toPredict,step,epsilon_init,order,rho)
    

    xOut = zeros(length(noise),1);
    err = xOut;
    w = zeros(order,length(noise)+1); % since order can be > 1 (here 2)
    nShift = zeros(order,length(noise)); % x(n-k)
    epsilon = (1/step)*ones(length(noise)+1,1);
    epsilon(1)=epsilon_init;
    beta=1;

    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        nShift(i,:) = [ zeros(1,i), noise(1: length(noise)-i)]; 
    end

    for k = 1: length(noise)
        
       % calculate the prediction
        xOut(k) = w(:,k)'*nShift(:,k);
        err(k) = toPredict(k)-xOut(k);
        % update
        w(:,k+1) = w(:,k) + ((beta*err(k))/(epsilon(k) +(nShift(:,k))'*nShift(:,k)))*nShift(:,k);
        if k>1
            epsilon(k+1) = epsilon(k) - (rho*step)*((err(k)*err(k-1)*(nShift(:,k)')*nShift(:,k-1))/((epsilon(k-1)+(nShift(:,k-1)')*nShift(:,k-1)))^2);
        end
    end
    w =  w(:,2:end);
    epsilon = epsilon(:,2:end);
end

