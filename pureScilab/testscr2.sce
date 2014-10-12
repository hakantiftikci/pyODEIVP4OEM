exec RK4.sci

function y = fy1(t,x,dx,u,px,pu,py)
    y = [x(1)+x(2);x(1)-x(2)];
endfunction
function u = fu1(t,x,pu)
    u = 10*[sin(t);cos(t)]
endfunction
function dx = fdx1(t,x,u,px)
    dx(1) = px(1)*x(1) + px(2)*x(2) + px(3)*x(1)*x(2) + px(7)*u(1)
    dx(2) = px(4)*x(1) + px(5)*x(2) + px(6)*(sin(x(1))+cos(x(2))) + px(8)*u(2)
    dx = dx(:);
endfunction


//px = [-1,0.2,0.0,  -2,0.5,0.0,0.7,0.8];
px = [-1,0.2,-0.0012,  -2,0.5,-0.03,  0.7,0.8];
t0 = 0.0;
t1 = 30.0;
dt = 0.1;
pu = [1];
py = [1];

nu = 2
ny = 2
nstore = int(t1/dt)

x0 = [1;2];
u0 = [1;1];
//function [xv,yv] = RK4(fdx,fu,fy,x0,t0,t1,dtpx,pu,py,nstore)
[xv,yv] = RK4(fdx1,fu1,fy1,x0,t0,t1,dt,px,pu,py,nstore)
