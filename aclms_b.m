function [g, h, err] = aclms_b(xIn, y, mu, order)

    N = length(xIn);
    g = complex(zeros(order,N));
    h = g;
    err = complex(zeros(1,N));
    
    xShift = zeros(N+order-1,1);
    xShift(order:N+order-1) = xIn;
    
    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut(k) = h(:,k)'*xShift(k:k+order-1) + g(:,k)'*conj(xShift(k:k+order-1));
        err(k) = y(k)-xOut(k);
        % update
        h(:,k+1) = h(:,k) + mu*conj(err(k))*xShift(k:k+order-1);
        g(:,k+1) = g(:,k) + mu*conj(err(k))*conj(xShift(k:k+order-1));
    end
    
    g = g(:,2:end);
    h =  h(:,2:end);
    

end