function [g, h, err] = aclms(xIn, y, mu, order)
    N = length(xIn);
    g = complex(zeros(order,N));
    h = g;
    err = complex(zeros(1,N));
    
    xShift = zeros(order,length(xIn)); % x(n-k)
    
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        xShift(i,:) = [ zeros(1,i), xIn(1: length(xIn)-i)]; 
    end
    
    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut(k) = h(:,k)'*xShift(:,k) + g(:,k)'*conj(xShift(:,k));
        err(k) = y(k)-xOut(k);
        % update
        h(:,k+1) = h(:,k) + mu*conj(err(k))*xShift(:,k);
        g(:,k+1) = g(:,k) + mu*conj(err(k))*conj(xShift(:,k));
    end
    g = g(:,2:end);
    h =  h(:,2:end);
    
end