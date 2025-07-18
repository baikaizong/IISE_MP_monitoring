module constants
    implicit none
    integer, parameter::dim=3,simu=10000,simu_OC=2000,n=20,window_size=50,ewma=50
end module constants

program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j,IC_samplesize
    real limit,arl,std,limitR,Limitl,lambdaewma
    real*8 dClkStart, dClkFinish
    real c_parameter,lambda_0(dim,dim),lambda_IC(dim,dim)
    open(10,file="./results/AMP_Poi-N20ICestimate.txt")
    lambdaewma=0.10
    call cpu_time(dClkStart)
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
    limit=  0.865         
    call parameter_setting(c_parameter,lambda_0)
    IC_samplesize=500
    call AMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std 
    IC_samplesize=1000
    call AMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    IC_samplesize=2000
    call AMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    IC_samplesize=5000
    call AMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    IC_samplesize=10000
    call AMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    call cpu_time(dClkFinish)  
    write(10,*) "Time is", dClkFinish-dClkStart, "seconds" 
    stop
end program main

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

subroutine sample_generation(sample_size,lambda_0,sample)
    use IMSL
    use constants
    implicit none
    integer sample_size
    real lambda_0(dim,dim)
    integer sample(dim,sample_size),samplegeneration(dim,dim,sample_size)
    integer i,j,k
    samplegeneration=0
    do i = 1, dim, 1
            do j = i, dim, 1
                if ( lambda_0(i,j)/=0.0 ) then
                    call RNPOI(sample_size,lambda_0(i,j),samplegeneration(i,j,:))
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

subroutine AMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
use IMSL
use constants
use Chartingmodule
implicit none
    integer i,j,k,IC_samplesize
    integer sample(window_size,dim,n),sampleIC(dim,IC_samplesize)
    real lambda(dim,dim),lambda_0(dim,dim),lambda_IC(dim,dim)
    real simuarl(simu),lambdaewma,lambdaco
    real Test
    real Limit,arl,st,std 

    i=1
    do while(i<=simu)	
    call sample_generation(IC_samplesize,lambda_0,sampleIC)	
    lambda_IC=lambda_0
    call MP_Estimation(IC_samplesize,sampleIC,lambda_IC)
    Test=0.0
    do j = 1, window_size, 1 
    call sample_generation(n,lambda_0,sample(j,:,:))
    end do
    do j = 1, ewma, 1
        sample(1:window_size-1,:,:)=sample(2:window_size,:,:)
        call sample_generation(n,lambda_0,sample(window_size,:,:))
        call LRT_test(dim,window_size,n,lambdaewma,sample,lambda_IC,Test)   
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
    call sample_generation(n,lambda_0,sample(window_size,:,:))
    call LRT_test(dim,window_size,n,lambdaewma,sample,lambda_IC,Test) 
    simuarl(i)=simuarl(i)+1.0
    end do
    print*,i
    print*,simuarl(i)
  !stop
    print*,"****************"
    i=i+1
    end do
    arl=(sum(simuarl))/(1.0*(simu-1))

   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu))
  std=std/(sqrt(dble(simu))) 
  return 
end subroutine AMP_estimate