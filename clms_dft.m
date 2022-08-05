
function [h, err] = clms_dft(xIn, y, mu,gamma,L)
   
    h = complex(zeros(L,length(y)));
    err = complex(zeros(1,length(y)));
 
    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut = h(:,k)'*xIn(:,k);
        err(k) = y(k)-xOut;
        % update
        h(:,k+1) = (1 - gamma*mu)*h(:,k) + mu*conj(err(k))*xIn(:,k);
     
    end
    h =  h(:,2:end);
    
end