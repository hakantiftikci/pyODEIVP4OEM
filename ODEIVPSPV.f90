!# -*- coding: utf-8 -*-
!"""
!Created on Fri Dec 06 22:06:59 2013
!
!@author: hakan
!"""

subroutine fdxsample1(t,x,u,p,dx,nx,nu,np)
    implicit none
    real, intent(in) :: t    
    
    integer, intent(in) :: nx,nu,np
    real, intent(in), dimension(nx) :: x
    real, intent(in), dimension(nu) :: u
    real, intent(in), dimension(np) :: p
    real, intent(out), dimension(nx) :: dx

    dx(1) = p(1)*x(1) + p(2)*x(2) + p(3)*x(1)*x(2) + p(7)*u(1)
    dx(2) = p(4)*x(1) + p(5)*x(2) + p(6)*(sin(x(1))+cos(x(2))) + p(8)*u(2)
end subroutine fdxsample1

subroutine fdxsample2(t,x,u,p,dx,nx,nu,np)
    implicit none
    real, intent(in) :: t    
    
    integer, intent(in) :: nx,nu,np
    real, intent(in), dimension(nx) :: x
    real, intent(in), dimension(nu) :: u
    real, intent(in), dimension(np) :: p
    real, intent(out), dimension(nx) :: dx

    real :: m,c,k
    m = p(1)
    c = p(2)
    k = p(3)

    dx(1) = x(2)
    dx(2) = (u(1)-c*x(2)-k*x(1))/m
end subroutine fdxsample2

subroutine stateDerivative(t,x,p,fu,fdx,dx,nx,nu,np)
    integer, intent(in) :: nx,nu,np
    real, intent(in) :: t
    real, intent(in), dimension(nx) :: x    
    real, intent(in), dimension(np) :: p
    real, intent(out), dimension(nx) :: dx

    external fdx
    external fu

    real, dimension(nu) :: u
    
    call fu(t,x,p,u,nx,nu,np)
    call fdx(t,x,u,p,dx,nx,nu,np)
end subroutine stateDerivative


!subroutine EulerIntegration(fdx,fu,fy,x0,t0,t1,dt,px,pu,py,nx,nu,ny,npx,npu,npy)
subroutine EulerIntegration(fdx,fu,fy,x0,tv,xv,yv,p,nx,nu,ny,np,nt)
    implicit none
    integer, intent(in) :: nx,nu,ny,np,nt !,nstore
    real, intent(in), dimension(nx) :: x0
    real, intent(in), dimension(nt)  :: tv
    real, intent(in), dimension(np) :: p
    real, intent(out), dimension(nt, nx) :: xv
    real, intent(out), dimension(nt, ny) :: yv
    external fdx,fu,fy,stateDerivative

    integer i
    real t,dt
    real, dimension(nx) :: x,dx
    real, dimension(nu) :: u
    real, dimension(ny) :: y

    t = tv(1)
    x = x0   

    !call fu(t,x,u,pu,nx,nu,np)
    !call fdx(t,x,u,px,dx,nx,nu,np)
    call stateDerivative(t,x,p,fu,fdx,dx,nx,nu,np)
    call fy(t,x,u,p,dx,y,nx,nu,np,ny)

    xv(1,:) = x
    yv(1,:) = y

    !do while (t<t1)        
    do i=1,(nt-1)
        !call fu(t,x,u,pu,nx,nu,npu)
        !call fdx(t,x,u,px,dx,nx,nu,npx)
        !call stateDerivative(t,x,dx,fu,fdx,px,pu,nx,nu,npx,npu)
        call stateDerivative(t,x,p,fu,fdx,dx,nx,nu,np)

        dt = tv(i+1)-tv(i)
        x = x + dt*dx
        t = t + dt

        !call fu(t,x,u,pu,nx,nu,npu)
        !call fdx(t,x,u,px,dx,nx,nu,npx)
        call stateDerivative(t,x,p,fu,fdx,dx,nx,nu,np)
        call fy(t,x,u,p,dx,y,nx,nu,np,ny)

        !if (i<nt) then
        xv(1+i,:) = x
        yv(1+i,:) = y
        !end if
        
        !print *,t,x,y
        !flush(6)
    end do

end subroutine EulerIntegration 

!subroutine RK4(fdx,fu,fy,x0,t0,t1,dt,px,pu,py,nx,nu,ny,npx,npu,npy)
!subroutine RK4(fdx,fu,fy,x0,tv,xv,yv,px,pu,py,nx,nu,ny,npx,npu,npy,nt)
subroutine RK4(fdx,fu,fy,x0,tv,xv,yv,p,nx,nu,ny,np,nt)
    implicit none
    integer, intent(in) :: nx,nu,ny,np,nt
    real, intent(in), dimension(nx) :: x0
    real, intent(in), dimension(nt)  :: tv
    real, intent(in), dimension(np) :: p
    real, intent(out), dimension(nt, nx) :: xv
    real, intent(out), dimension(nt, ny) :: yv
    external fdx,fu,fy,stateDerivative

    real, dimension(4) :: c 
    real, dimension(4) :: b
    !real, dimension(4,4) :: a = (/ (/0.0, 0.0, 0.0, 0.0/) ,  (/1.0/2.0,0.0,0.0,0.0/) , (/0.0,1.0/2.0,0.0,0.0/) , (/0.0,0.0,1.0,0.0/) /)
    real, dimension(4,4) :: a 

    real, dimension(nx,4) :: k
    real, dimension(nx) :: xstage

    integer i,j,m
    !integer istore
    real t,dt
    real, dimension(nx) :: x,dx,deltax
    real, dimension(nu) :: u
    real, dimension(ny) :: y

    !a = reshape( (/0.0, 0.0, 0.0, 0.0 , 1.0/2.0,0.0,0.0,0.0 , 0.0,1.0/2.0,0.0,0.0 , 0.0,0.0,1.0,0.0/), (/4,4/))
    ! k(i) = f(t+c(i)*h, y(n)+sum(a(i,j)*k(j)))
    ! y(n+1) = y(n) + h*sum(b(i)*k(i))
    a(1,:) = (/0.0, 0.0, 0.0, 0.0/) 
    a(2,:) = (/1.0/2.0,0.0,0.0,0.0/) 
    a(3,:) = (/0.0,1.0/2.0,0.0,0.0/) 
    a(4,:) = (/0.0,0.0,1.0,0.0/)
    b = (/ 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 /)
    c = (/ 0.0, 1.0/2.0, 1.0/2.0, 1.0 /)

    !print *,"a=",a
    !print *,"b=",b
    !print *,"c=",c

    t = tv(1)
    x = x0

    !call fu(t,x,u,pu,nx,nu,npu)
    !call fdx(t,x,u,px,dx,nx,nu,npx)
    call stateDerivative(t,x,p,fu,fdx,dx,nx,nu,np)
    call fy(t,x,u,p,dx,y,nx,nu,np,ny)

    xv(1,:) = x
    yv(1,:) = y

    !do while (t<t1)
    do m=1,(nt-1)
        !call fu(t,x,u,pu,nx,nu,npu)
        !call fdx(t,x,u,px,dx,nx,nu,npx)
        dt = tv(i+1)-tv(i)

        deltax = 0
        do i=1,4
            xstage = x
            do j=1,i-1
                xstage = xstage+dt*a(i,j)*k(:,j)
            end do
            !call stateDerivative(t+c(i)*dt,xstage,k(:,i),fu,fdx,px,pu,nx,nu,npx,npu)
            call stateDerivative(t+c(i)*dt,xstage,p,fu,fdx,k(:,i),nx,nu,np)
            !subroutine stateDerivative(t,x,dx,fu,fdx,px,pu,nx,nu,npx,npu)
            deltax = deltax + dt*b(i)*k(:,i)
        end do
    
        x = x + deltax
        t = t + dt
        
        call stateDerivative(t,x,p,fu,fdx,dx,nx,nu,np)
        call fy(t,x,u,p,dx,y,nx,nu,np,ny)

       !if (istore<=nstore) then
        xv(1+m,:) = x
        yv(1+m,:) = y
        !end if

        !print *,t,x,y
        !flush(6)
    end do

end subroutine RK4 
