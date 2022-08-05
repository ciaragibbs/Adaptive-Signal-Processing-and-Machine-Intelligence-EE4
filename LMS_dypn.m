function [xOut,w,err,alphas] = LMS_dypn(xIn,step,order,alpha,a_update,winit,b)
    
    % b is the bias, either equal to 0 or 1

    xOut = zeros(length(xIn),1);
    err = xOut;
    w = zeros(order + b,length(xIn)+1); % since order can be > 1 (here 2)
    w(:,1) = winit;
    xShift = zeros(order,length(xIn)); % x(n-k)
    alphas = alpha;
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        xShift(i,:) = [ zeros(1,i), xIn(1: length(xIn)-i)']; 
    end
    
    if b
        xShift = [ones(1,length(xIn)); xShift];
    end
    
    for k = 1: length(xIn)
        
        % calculate the prediction
        xOut(k) = alpha*tanh(w(:,k)'*xShift(:,k));
        err(k) = xIn(k)-xOut(k);
        act_function = alpha*(1-(xOut(k)/alpha)^2);
        % update
        w(:,k+1)=w(:,k)+(step*act_function*err(k)).*xShift(:,k);
        if a_update
            alpha = alpha + 0.3*err(k)*(xOut(k)/alpha);
            alphas = [alphas,alpha];
        else
            % nothing, alphas remains as alpha
        end
    end
    w =  w(:,2:end);
end




