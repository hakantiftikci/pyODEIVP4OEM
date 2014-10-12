function xdot=f(t, x)
    //px = [-1,0.2,0.0,  -2,0.5,0.0,0.7,0.8];
    px = [-1,0.2,-0.0012,  -2,0.5,-0.03,  0.7,0.8];
    u = 10*[sin(t);cos(t)]
    xdot(1) = px(1)*x(1) + px(2)*x(2) + px(3)*x(1)*x(2) + px(7)*u(1)
    xdot(2) = px(4)*x(1) + px(5)*x(2) + px(6)*(sin(x(1))+cos(x(2))) + px(8)*u(2)
endfunction

y0=[1;2];
t0=0;
t=0:0.1:30;
y = ode(y0,t0,t,f);
figure;plot(t,y,'r.-')
figure;plot(y(1,:),y(2,:))

tyvec = [t',y']
//save ty tyvec
unix('del TY.dat');
write('TY.dat',tyvec)
