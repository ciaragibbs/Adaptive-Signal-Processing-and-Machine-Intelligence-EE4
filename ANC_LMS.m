function [noiseEst,w,xHat] = ANC_LMS(xIn,ref,step,leak,order)

    noiseEst = zeros(length(xIn),1);
    xHat = noiseEst;
    w = zeros(order,length(xIn)+1); % since order can be > 1 (here 2)
    uShift = zeros(order,length(xIn)); % x(n-k)
    
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        uShift(i,:) = [ zeros(1,i-1), ref(1: length(xIn)-(i-1))]; 
    end
    
    for k = 1: length(xIn)
        
        % calculate the prediction
        noiseEst(k) = w(:,k)'*uShift(:,k);
        xHat(k) = xIn(k) - noiseEst(k);
        % update
        w(:,k+1) = (1 - step*leak).*w(:,k) + step*xHat(k)*uShift(:,k);
  
    end
    w =  w(:,2:end);


end
