function [h, err] = clms(xIn, y, mu, order)
    N = length(xIn);
    h = complex(zeros(order,N));
    err = complex(zeros(1,N));
    
    xShift = zeros(order,length(xIn)); % x(n-k)
    
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        xShift(i,:) = [ zeros(1,i), xIn(1: length(xIn)-i)]; 
    end
    
    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut(k) = h(:,k)'*xShift(:,k);
        err(k) = y(k)-xOut(k);
        % update
        h(:,k+1) = h(:,k) + mu*conj(err(k))*xShift(:,k);
        % ^^ for part 1b, the leak is equal to 0 and hence it because the
        % conventional linear adapative LMS filter
    end
    h =  h(:,2:end);
    
end
