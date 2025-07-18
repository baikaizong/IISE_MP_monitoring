module constants
    implicit none
    integer, parameter::dim=3,simu=10000,n=20,window_size=50,ewma=50
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
    open(10,file="./results/SMP_Poi-N20ICestimate.txt")
    lambdaewma=0.10
    call cpu_time(dClkStart)
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
    limit=  1.313         
    call parameter_setting(c_parameter,lambda_0)
    IC_samplesize=500
    call SMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std 
    IC_samplesize=1000
    call SMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    IC_samplesize=2000
    call SMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    IC_samplesize=5000
    call SMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
    write(10,*),"IC_samplesize is",IC_samplesize
    write(10,*),"IC_arl is",arl    
    write(10,*),"IC_std is",std
    IC_samplesize=10000
    call SMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
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

subroutine SMP_estimate(IC_samplesize,lambdaewma,lambda_0,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j,IC_samplesize
    integer sample(dim,n),sampleIC(dim,IC_samplesize)
    real lambda_0(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu)
    real Limit,Test,arl,st,std,lambdaewma
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    

    i=1
    do while(i<=simu)
    call sample_generation(IC_samplesize,lambda_0,sampleIC)
    call MM_Estimation(dim,IC_samplesize,sampleIC,miu_exp,Cov_exp)
    Inverse_Cov_exp= .i. Cov_exp
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp 
    do j = 1, ewma, 1 
    call sample_generation(n,lambda_0,sample)
    call MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcov(dim,n,miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif     
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    call sample_generation(n,lambda_0,sample)
    call MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcov(dim,n,miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
      !print*,Test
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
  std=sqrt(st/(simu-1.0))
  std=std/(sqrt(dble(simu)-1)) 
  return 
end subroutine SMP_estimate