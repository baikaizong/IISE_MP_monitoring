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
    real limit,arl,std,lambdaewma
    real limitl,limitR,miu_exp(dim),Cov_exp(dim,dim),lambda_0(dim,dim)
    real dClkStart, dClkFinish
    real c_parameter
    open(10,file="./results/SMP_Poi.txt")
    call cpu_time(dClkStart)
    lambdaewma=0.10
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
    !limit=  1.313        
    !arl= 201.673   
    limitl=1.3
    limitR=1.7
    call parameter_setting(c_parameter,lambda_0)
    call Exp_meancov(lambda_0,miu_exp,Cov_exp) 
    call SMP_limitsearch(lambdaewma,lambda_0,miu_exp,Cov_exp,limitl,limitR,limit,arl,std)
    write(10,*),"Limit is",limit
    write(10,*),"IC_arl is",arl   
    call SMP_Poi(lambdaewma,0.05,0.20,limit,lambda_0,miu_exp,Cov_exp) 
    write(10,*),"===============c_parameter is 2.0========================="    
    c_parameter=2.0
    !limit=1.050       
    !arl= 196.331   
    limitl=1.1
    limitR=1.6
    call parameter_setting(c_parameter,lambda_0) 
    call Exp_meancov(lambda_0,miu_exp,Cov_exp) 
    call SMP_limitsearch(lambdaewma,lambda_0,miu_exp,Cov_exp,limitl,limitR,limit,arl,std)
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl
    call SMP_Poi(lambdaewma,0.04,0.20,limit,lambda_0,miu_exp,Cov_exp) 
    write(10,*),"===============c_parameter is 5.0========================="
    c_parameter=5.0
    !limit=1.000        
    !arl= 198.216 
    limitl=1.0
    limitR=1.5
    call parameter_setting(c_parameter,lambda_0)
    call Exp_meancov(lambda_0,miu_exp,Cov_exp) 
    call SMP_limitsearch(lambdaewma,lambda_0,miu_exp,Cov_exp,limitl,limitR,limit,arl,std)    
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl
    call SMP_Poi(lambdaewma,0.02,0.20,limit,lambda_0,miu_exp,Cov_exp) 
    call cpu_time(dClkFinish)  
    write(10,*) "Time is", dClkFinish-dClkStart, "seconds" 
    stop
end program main

subroutine SMP_Poi(lambdaewma,increment_meanshift,increment_covshift,limit,lambda_0,miu_exp,Cov_exp)
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real lambdaewma,increment_meanshift,increment_covshift
    real limit,arl,std
    real lambda_0(dim,dim),miu_exp(dim),Cov_exp(dim,dim),shift(dim,dim)   
    write(10,*),"++++++++++++++++shifts in miu_1+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 9, 1
        if (i==5) then
           cycle
        endif
        shift(1,1)=(lambda_0(1,1)+lambda_0(1,2)+lambda_0(1,3))*increment_meanshift*(i-5)
        call SMP_OCARL(lambdaewma,lambda_0,miu_exp,Cov_exp,Limit,shift,arl,std)
        write(10,*),"shift is",increment_meanshift*(i-5)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
    write(10,*),"++++++++++++++++shifts in theta_12+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 9, 1
         if (i==5) then
            cycle
        endif

        shift(1,2)=lambda_0(1,2)*increment_covshift*(i-5)
        shift(2,1)=shift(1,2)
        shift(1,1)=-shift(1,2)
        shift(2,2)=-shift(1,2)
        call SMP_OCARL(lambdaewma,lambda_0,miu_exp,Cov_exp,Limit,shift,arl,std)
        write(10,*),"shift is",increment_covshift*(i-5)
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do    
end subroutine SMP_Poi

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

subroutine Exp_meancov(lambda_0,miu_exp,Cov_exp)
    use constants
    implicit none
    integer i,j
    real lambda_0(dim,dim),Cov_exp(dim,dim),miu_exp(dim)
    Cov_exp=lambda_0
    do i = 1, dim, 1
        do j = 1, dim, 1
            if (j/=i)then
            Cov_exp(i,i)=Cov_exp(i,i)+Cov_exp(i,j)
            end if          
        end do
    miu_exp(i)=Cov_exp(i,i)
    end do
    return
end subroutine Exp_meancov


subroutine SMP_OCARL(lambdaewma,lambda_0,miu_exp,Cov_exp,Limit,shift,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real lambda(dim,dim),lambda_0(dim,dim),shift(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu_OC)
    real Limit,Test,arl,st,std,lambdaewma
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    
    lambda=lambda_0+shift

    Inverse_Cov_exp= .i. Cov_exp
    i=1
    do while(i<=simu_OC)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call sample_generation(lambda_0,sample)
    call MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcov(dim,n,miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif   
    !print*,Test
  !stop
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    call sample_generation(lambda,sample)
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
    arl=(sum(simuarl))/(1.0*(simu_OC))

   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu_OC
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu_OC))
  std=std/(sqrt(dble(simu_OC))) 
  return 
end subroutine SMP_OCARL

subroutine SMP_limitsearch(lambdaewma,lambda_0,miu_exp,Cov_exp,limitl,limitR,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real lambda_0(dim,dim),shift(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu_IC)
    real Limit,Test,arl,st,std,limitl,limitR,lambdaewma
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim)
    arl=0.0
    Inverse_Cov_exp= .i. Cov_exp
    do while (abs(arl-200.0) >= 2.0)
    Limit = (LimitL + LimitR) / 2.0
    i=1
    do while(i<=simu_IC)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call sample_generation(lambda_0,sample)
    call MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcov(dim,n,miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    !print*,Test
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif   
    !print*,Test
  !stop
    simuarl(i)=0.0
    !print*,"OK"
    do while (Test<Limit)
    call sample_generation(lambda_0,sample)
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
end subroutine SMP_limitsearch