program MC
    implicit none
    real*8,parameter::PI=3.141592653589793239d0
    integer,parameter::Nx=4
    integer,parameter::NL=2500000    !2500s
    integer,parameter::Nt=500000     !500s
    real*8,parameter::h=0.001d0
    integer,parameter::N=1    
    integer,parameter::nbin=100
    integer,parameter::mm=Nx*(NL-Nt)/N
    
    real*8,parameter::x_min=-3.0d0,x_max=3.0d0
    real*8::t,hx,hd,tau       
    real*8::U1,U2,noise,x0   
    integer::i,j,k,ii,j1,m,ij,i1,i2
    real*8::k1,k2,k3,k4,Y1,Y2,count1
    real*8::time_begin,time_end,total_time
    real*8,dimension(-10000:NL)::x
    real*8::delta_x,delta_y,j2
    real*8,dimension(1:mm)::tau_sum,xx_sum,Y1_sum,Y2_sum
    real*8,dimension(1:nbin)::tau_m,xx_m,Y1_m,Y2_m
    real*8,dimension(1:nbin)::tau_all,xx_all,Y1_all,Y2_all
    real*8,dimension(1:nbin+1)::xbin,ybin
    real*8,dimension(1:nbin)::wbin

    hx=(x_max-x_min)/real(Nx)
    hd=(x_max-x_min)/real(nbin)
    open(100,file='Sample_delay2_2.txt')
    open(200,file='Total_time.txt')
    call cpu_time(time_begin)
    call random_seed()
    ! 对x进行划分小块
    do i=1,nbin
        delta_x=(x_max-x_min)/nbin
        xbin(i)=x_min+(i-1)*delta_x    
    end do
    xbin(nbin+1)=x_max
    !*********求漂移和扩散系数****************************
    j1=0
!    do k=1,1   !tau=0.2(6),0.5(21),0.8(36),0.1(1)
!        tau=0.1+(k-1)*0.02d0
        tau=0.5d0
        m=int(tau/h)
        ij=1
        xx_sum=0.0
        Y1_sum=0.0
        Y2_sum=0.0
        count1=0
        wbin=0
        do j=1,Nx,1  
            x0=x_min+(j-0.5)*hx

            t=0.d0
            do ii=-m,0
                x(ii)=x0
            end do
            do i=1,NL,1
                call random_number (U1)  
                call random_number (U2) 
                noise=dsqrt(-2.0d0*dlog(U1)/h)*dcos(2.0d0*PI*U2)

                call FCN(k1,t,x(i-1),x(i-1-m),noise)
                call FCN(k2,t+h/2,x(i-1)+h*k1,x(i-1-m)+h*k1,noise)

                x(i)=x(i-1)+(k1+k2)*h/2.d0
                t=t+h

                if(i>Nt.and.mod(i,N)==0)then
                    if(x(i-N)>x_min.and.x(i-N)<=x_max)then
                        Y1=(x(i)-x(i-N))/h/N
                        Y2=(x(i)-x(i-N))**2/h/N
    !                    write(100,"(4f30.20)") tau,x(i-N),Y1,Y2
                        ij=floor((x(i-N)-x_min)/hd)+1
                        xx_sum(ij)=xx_sum(ij)+x(i-N)
                        Y1_sum(ij)=Y1_sum(ij)+Y1
                        Y2_sum(ij)=Y2_sum(ij)+Y2
                        wbin(ij)=wbin(ij)+1
                        count1=count1+1                
                    else
                        j1=j1+1
                    end if
                endif
            end do   !*****time*******
            write(*,*) j,j1
        end do
        
        ! 对x进行划分小块,求漂移系数和扩散系数   
        do i1=1,nbin            
            xx_m(i1)=xx_sum(i1)/wbin(i1)
            Y1_m(i1)=Y1_sum(i1)/wbin(i1)
            Y2_m(i1)=Y2_sum(i1)/wbin(i1)
            wbin(i1)=wbin(i1)/(count1/nbin)
            write(100,"(4f30.20)") xx_m(i1),Y1_m(i1),Y2_m(i1),wbin(i1)
        !    write(100,"(4f30.20)") tau_all(i)/j2,xx_all(i)/j2,Y1_all(i)/j2，Y2_all(i)/j2
            write(*,*) i1
        end do

!        write(*,*) k,j1
!    end do
    !*************************************************
    call cpu_time(time_end)
    total_time=time_end-time_begin
    write(200,"(f10.6)") total_time
    !******************************************
    contains
    !*******************************************************
    subroutine FCN(YPRIME,tt,XX,YY,whitenoise)
        implicit none
        real*8 YPRIME,tt,XX,YY,whitenoise
        real*8::DD=0.02d0
!        YPRIME=XX*(1-YY)+dsqrt(2*DD)*XX*whitenoise
        YPRIME=XX-XX**3+dsqrt(2*DD)*YY*whitenoise 
    end subroutine    
end program