function [xOut,w,err] = LMS(xIn,step,leak,order)

    xOut = zeros(length(xIn),1);
    err = xOut;
    w = zeros(order,length(xIn)+1); % since order can be > 1 (here 2)
    xShift = zeros(order,length(xIn)); % x(n-k)
    
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        xShift(i,:) = [ zeros(1,i), xIn(1: length(xIn)-i)']; 
    end
    
    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut(k) = w(:,k)'*xShift(:,k);
        err(k) = xIn(k)-xOut(k);
        % update
        w(:,k+1) = (1 - step*leak).*w(:,k) + step*err(k)*xShift(:,k);
        % ^^ for part 1b, the leak is equal to 0 and hence it because the
        % conventional linear adapative LMS filter
    end
    w =  w(:,2:end);
end



