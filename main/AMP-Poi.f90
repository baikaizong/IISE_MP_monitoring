module constants
    implicit none
    integer, parameter::dim=3,simu_IC=10000,simu_OC=2000,n=20,window_size=50,ewma=50
end module constants
program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real limit,arl,std,limitR,Limitl,lambdaewma
    real*8 dClkStart, dClkFinish
    real c_parameter,lambda_0(dim,dim)
    open(10,file="./results/AMP_Poi-N20lambda-0.10.txt")
    lambdaewma=0.10
    call cpu_time(dClkStart)
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
    !limit=  0.8625         
	limitl=0.80
	limitR=0.88
	call parameter_setting(c_parameter,lambda_0)
	call AMP_limitsearch(lambdaewma,lambda_0,limitl,limitR,Limit,arl,std)   
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl     
    call AMP_Poi(0.05,0.20,lambdaewma,limit,lambda_0) 
    write(10,*),"===============c_parameter is 2.0========================="  
	stop  
    c_parameter=2.0
    !limit=0.8495       
	limitl=0.80
	limitR=0.88
	call parameter_setting(c_parameter,lambda_0)
	call AMP_limitsearch(lambdaewma,lambda_0,limitl,limitR,Limit,arl,std)    
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl
    call AMP_Poi(0.04,0.20,lambdaewma,limit,lambda_0)
	stop 
    write(10,*),"===============c_parameter is 5.0========================="
    c_parameter=5.0
    !limit=0.825       
	limitl=0.80
	limitR=0.88
    call parameter_setting(c_parameter,lambda_0)
	call AMP_limitsearch(lambdaewma,lambda_0,limitl,limitR,Limit,arl,std)       
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl
    call AMP_Poi(0.02,0.20,lambdaewma,limit,lambda_0) 
    call cpu_time(dClkFinish)  
    write(10,*) "Time is", dClkFinish-dClkStart, "seconds" 
    stop
end program main

subroutine AMP_Poi(increment_meanshift,increment_covshift,lambdaewma,limit,lambda_0)
    use IMSL
    use constants
    implicit none
    integer i,j
    real increment_meanshift,increment_covshift,lambdaewma
    real limit,arl,std
    real lambda_0(dim,dim),shift(dim,dim),c_parameter 
   
    write(10,*),"++++++++++++++++shifts in miu_1+++++++++++++++++++++++"
    shift=0.0     
    do i = 1, 10, 1
        if (i==5) then
           cycle
        endif
        shift(1,1)=(lambda_0(1,1)+lambda_0(1,2)+lambda_0(1,3))*increment_meanshift*i
        call AMP_OCARL(lambdaewma,lambda_0,Limit,shift,arl,std)
        write(10,*),"shift is",increment_meanshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
    write(10,*),"++++++++++++++++shifts in theta_12+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 10, 1
        if (i==5) then
            cycle
        endif

        shift(1,2)=lambda_0(1,2)*increment_covshift*i
        shift(2,1)=shift(1,2)
		shift(1,1)=-shift(1,2)
		shift(2,2)=-shift(1,2)
        call AMP_OCARL(lambdaewma,lambda_0,Limit,shift,arl,std)
        write(10,*),"shift is",increment_covshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do    
end subroutine AMP_Poi

subroutine parameter_setting(c_parameter,lambda_0)
    use constants
    implicit none
    real lambda_0(dim,dim),c_parameter
    integer i,j
    do i = 1, dim, 1
        do j = 1, dim, 1
            lambda_0(i,j)=c_parameter/(abs(i-j)+1)
        end do       
    end do 
    return 
end subroutine parameter_setting

subroutine sample_generation(lambda_0,sample)
    use IMSL
    use constants
    implicit none
    real lambda_0(dim,dim)
    integer sample(dim,n),samplegeneration(dim,dim,n)
    integer i,j,k
    samplegeneration=0
    do i = 1, dim, 1
            do j = i, dim, 1
                if ( lambda_0(i,j)/=0.0 ) then
                    call RNPOI(n,lambda_0(i,j),samplegeneration(i,j,:))
                end if
                samplegeneration(j,i,:)=samplegeneration(i,j,:)          
            end do        
    end do
    sample=0.0
    do i = 1, dim, 1
        do j = 1, dim, 1
        sample(i,:)= sample(i,:)+samplegeneration(i,j,:)          
        end do
    end do
    return  
end subroutine sample_generation

subroutine AMP_OCARL(lambdaewma,lambda_0,Limit,shift,arl,std)
use IMSL
use constants
use Chartingmodule
implicit none
    integer i,j
    integer sample(window_size,dim,n)
    real lambda(dim,dim),lambda_0(dim,dim),shift(dim,dim)
    real simuarl(simu_OC),lambdaewma,lambdaco
    real Test
    real Limit,arl,st,std  
    lambda=lambda_0+shift

    i=1
    do while(i<=simu_OC)
    Test=0.0
    do j = 1, window_size, 1 
    call sample_generation(lambda_0,sample(j,:,:))
    end do
    do j = 1, ewma, 1
        sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
        call sample_generation(lambda_0,sample(window_size,:,:))
        call LRT_test(dim,window_size,n,lambdaewma,sample,lambda_0,Test)  
        if (Test>Limit) then
            exit
        endif  
    end do
    if (Test>Limit) then
        cycle
    endif
    simuarl(i)=0.0
    do while (Test<Limit)
    sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
    call sample_generation(lambda,sample(window_size,:,:))
    call LRT_test(dim,window_size,n,lambdaewma,sample,lambda_0,Test)
    simuarl(i)=simuarl(i)+1.0
    end do
    print*,i
    print*,simuarl(i)
  !stop
    print*,"****************"
    i=i+1
    end do
    arl=(sum(simuarl))/(1.0*(simu_OC))

   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu_OC
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu_OC))
  std=std/(sqrt(dble(simu_OC))) 
  return 
end subroutine AMP_OCARL

subroutine AMP_limitsearch(lambdaewma,lambda_0,limitl,limitR,Limit,arl,std)
use IMSL
use constants
use Chartingmodule
implicit none
    integer i,j
    integer sample(window_size,dim,n)
    real shift(dim,dim),lambda_0(dim,dim)
    real lambdaewma,lambdaco
    real simuarl(simu_IC)
    real Limit,arl,st,std,LimitL,LimitR
    real Test
    arl=0.0
    do while (abs(arl-200.0) >=2.0)
    Limit = (LimitL + LimitR) / 2.0
    i=1
    do while(i<=simu_IC)
    Test=0.0
    do j = 1, window_size, 1 
       call sample_generation(lambda_0,sample(j,:,:))
    end do
    do j = 1, ewma, 1
        sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
        call sample_generation(lambda_0,sample(window_size,:,:))
        call LRT_test(dim,window_size,n,lambdaewma,sample,lambda_0,Test)  
        if (Test>Limit) then
            exit
        endif  
    end do

    !print*,Test
    if (Test>Limit) then
        cycle
    endif   
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
    call sample_generation(lambda_0,sample(window_size,:,:))
    call LRT_test(dim,window_size,n,lambdaewma,sample,lambda_0,Test)
    !print*,Test
    simuarl(i)=simuarl(i)+1.0
    end do
    print*,i
    print*,simuarl(i)
  !stop
    print*,"****************"
    i=i+1
    end do
    arl=(sum(simuarl))/(1.0*(simu_IC))
    write(10,*),"limit is",limit
    write(10,*),"arl is",arl
    if (arl < 200.0) then
      LimitL = Limit
     else
      LimitR = Limit
     end if
    enddo  
   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu_IC
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu_IC))
  std=std/(sqrt(dble(simu_IC))) 
  return 
end subroutine AMP_limitsearch