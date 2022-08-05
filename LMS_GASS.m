function [xOut,w,err,mus] = LMS_GASS(noise,toPredict,step,rho,alpha,order,classify)
    leak =0;
    xOut = zeros(length(noise),1);
    err = xOut;
    w = zeros(order,length(noise)+1); % since order can be > 1 (here 2)
    nShift = zeros(order,length(noise)); % x(n-k)
    psi = w;
    
    if classify < 3
        mus = step*ones(1,length(noise)+1);
    else
        mus = zeros(1,length(noise)+1);
        mus(1) = step;
    end
    
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        nShift(i,:) = [ zeros(1,i), noise(1: length(noise)-i)]; 
    end
    
    for k = 1: length(noise)
        
        % calculate the prediction
        xOut(k) = w(:,k)'*nShift(:,k);
        err(k) = toPredict(k)-xOut(k);
        % update
        w(:,k+1) = (1 - step*leak).*w(:,k) + mus(k)*err(k)*nShift(:,k);
        
        
        switch classify
            case 1
                % don't do anything
            case 2
                % don't do anything
            case 3 % benveniste
                % now adding an update for mu
                mus(k+1)= mus(k) + rho*err(k)*nShift(:,k)'*psi(:,k);
                psi(:,k+1) = (eye(order,order) - (mus(k)*nShift(:,k)'*(nShift(:,k))))*psi(:,k) + err(k)*nShift(:,k);
            case 4 % ang and farhang
                % now adding an update for mu
                mus(k+1)= mus(k) + rho*err(k)*nShift(:,k)'*psi(:,k);
                psi(:,k+1) = alpha*psi(:,k) + err(k)*nShift(:,k);
            otherwise
                % now adding an update for mu
                mus(k+1)= mus(k) + rho*err(k)*nShift(:,k)'*psi(:,k);
                psi(:,k+1) = err(k)*nShift(:,k);
        end
        
    end
    w =  w(:,2:end);
    mus = mus(:,2:end);
end
