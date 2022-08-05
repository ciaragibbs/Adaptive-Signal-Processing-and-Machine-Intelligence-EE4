function [h, err] = clms_b(xIn, y, mu, order)
    N = length(xIn);
    h = complex(zeros(order,N));
    err = complex(zeros(1,N));
    
    xShift = zeros(N+order-1,1);
    xShift(order:N+order-1) = xIn;

    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut(k) = h(:,k)'*xShift(k:k+order-1);
        err(k) = y(k)-xOut(k);
        % update
        h(:,k+1) = h(:,k) + mu*conj(err(k))*xShift(k:k+order-1);
  
    end
    h =  h(:,2:end);
    
end
