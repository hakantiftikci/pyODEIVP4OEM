function [xv,yv] = RK4(fdx,fu,fy,x0,t0,t1,dt,px,pu,py,nstore)
    
    function dx=stateDerivative(t,x,fu,fdx,px,pu)
        u = fu(t,x,pu)
        dx = fdx(t,x,u,px)
    endfunction
    
    a = [[0.0, 0.0, 0.0, 0.0]; ...
        [1.0/2.0,0.0,0.0,0.0]; ...
        [0.0,1.0/2.0,0.0,0.0]; ...
        [0.0,0.0,1.0,0.0]];
    b = [ 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 ];
    c = [ 0.0, 1.0/2.0, 1.0/2.0, 1.0 ];
    
    t = t0
    x = x0
    dx = 0*x0;
    u0 = fu(t0,x0,pu)
    y0 = fy(t0,x0,dx,u0,px,pu,py)
    
    
    xv = zeros(nstore, length(x0))
    yv = zeros(nstore, length(x0))
    y = y0

    istore = 1

    while (t<t1)
        if (istore<nstore) then
            //disp(xv)
            //disp(x)
            xv(istore,:) = x'
            yv(istore,:) = y'
            istore = istore+1
        end 

        deltax = 0
        for i=1:4
            xstage = x
            for j=1:(i-1)
                xstage = xstage+dt*a(i,j)*k(:,j)
            end 
            k(:,i) = stateDerivative(t+c(i)*dt,xstage,fu,fdx,px,pu)
            
            deltax = deltax + dt*b(i)*k(:,i)
        end
    
        x = x + deltax
        t = t + dt
        u = fu(t,x,pu)
        y = fy(t,x,dx,u,px,pu,py)
    end

endfunction
