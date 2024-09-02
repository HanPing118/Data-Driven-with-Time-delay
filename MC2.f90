program monte_carlo_pdf
    implicit none
    character(14)::filename
    real*8,parameter::PI=3.141592653589793239d0
    integer,parameter::Nx=2000000
    integer,parameter::Nd=100
    real*8,parameter::h=0.01        

    integer,parameter::NL=1000  !tn=500
    integer,parameter::Nt=999
    integer,parameter::Ns=Nx*(NL-Nt)

    integer::m,i,j,k,dx,kk,j1
    real*8::t,t1,noise,hx,hd,x0
    real*8 k1,k2,U1,U2
    real*8::time_begin,time_end
    real*8,parameter::x_min=0d0,x_max=10.0d0

    real*8,dimension(0:NL)::x 
    real*8,dimension(1:Nd)::PDF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    hx=(x_max-x_min)/real(Nx)
    hd=(x_max-x_min)/real(Nd)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call cpu_time(time_begin)
    call random_seed()

    open(unit=100,file='tumor_test.txt')
    PDF=0.0d0
    j1=0
!do j=1,tra,1
    do dx=1,Nx,1
        x0=x_min+(dx-0.5)*hx
        !****初值的处理***************************
        x(0)=x0
        t=0.0d0
        do i=1,NL,1
            call random_number(U1) 
            call random_number(U2) 	
            noise=dsqrt(-2.0d0*dlog(U1)/h)*dsin(2.0d0*PI*U2)  

            call FCN(k1,t,x(i-1),noise)
			call FCN(k2,t+h,x(i-1)+h*k1,noise)

            x(i)=x(i-1)+(k1+k2)*h/2.d0
            t=t+h

            if(i>Nt)then
                if ((x(i)>=x_min).and.(x(i)<x_max)) then  
                    m=floor((x(i)-x_min)/hd)+1
                    PDF(m)=PDF(m)+1.0d0  
                else
                    j1=j1+1
                end if 

            end if 
        end do 
        if(mod(dx,100000)==0)then
        	write(*,*) dx/100000,j1
        end if
    end do   !do dx=1,Nx
!    if(mod(j,1000)==0)then
!	    write(*,*) j/1000,j1
!    end if
!end do
    !********************************************************
    do k=1,Nd 
        PDF(k)=PDF(k)/hd/(Ns-j1)
        write(100,"(2f35.20)") x_min+(k-0.5)*hd,PDF(k) 
    enddo
    !********************************************************
!end do
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    call cpu_time(time_end)
    write(*,*)"计算时间为：",time_end-time_begin
    contains      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine FCN(YPRIME,tt,XX,whitenoise)
        implicit none
        real*8 YPRIME,tt,XX,whitenoise
        real*8::DD=0.05d0,theta=0.1d0,beta=2.0d0
        !***************Original system***********************************************************
        YPRIME=XX*(1-theta*XX)-beta*XX**2/(1+XX**2)+dsqrt(2*DD)*whitenoise
!        YPRIME=XX*(1-YY)-beta+dsqrt(2*DD)*XX*whitenoise
!        YPRIME=XX-XX**2+dsqrt(2*DD)*YY*whitenoi
    end subroutine
end program
