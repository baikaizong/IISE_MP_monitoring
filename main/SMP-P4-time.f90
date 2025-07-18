module constants
    implicit none
    integer, parameter::n=20,simu_time=10000,dim=4,ewma=50
end module constants
program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    real max_time,min_time,ave_time,std_time
    real c_parameter,lambda_0(dim,dim),lambdaewma
    real dClkStart, dClkFinish
    open(10,file="results/SMP_P4_time.txt")
    lambdaewma=0.10
    call cpu_time(dClkStart)
    c_parameter=0.5  
    call parameter_setting(c_parameter,lambda_0)  
    call SMP_Time_count(lambdaewma,lambda_0,max_time,min_time,ave_time,std_time) 
    write(10,*),"c_parameter is",c_parameter
    write(10,*),"max_time is",max_time
    write(10,*),"min_time is",min_time
    write(10,*),"average_time is",ave_time
    write(10,*),"standard erre_time is",std_time
    write(10,*),"========================================"
    c_parameter=2.0 
    call parameter_setting(c_parameter,lambda_0)   
    call SMP_Time_count(lambdaewma,lambda_0,max_time,min_time,ave_time,std_time) 
    write(10,*),"c_parameter is",c_parameter
    write(10,*),"max_time is",max_time
    write(10,*),"min_time is",min_time
    write(10,*),"average_time is",ave_time
    write(10,*),"standard erre_time is",std_time
    write(10,*),"========================================"
    c_parameter=5.0    
    call parameter_setting(c_parameter,lambda_0)
    call SMP_Time_count(lambdaewma,lambda_0,max_time,min_time,ave_time,std_time) 
    write(10,*),"c_parameter is",c_parameter
    write(10,*),"max_time is",max_time
    write(10,*),"min_time is",min_time
    write(10,*),"average_time is",ave_time
    write(10,*),"standard erre_time is",std_time
    write(10,*),"========================================"
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
    
end subroutine Exp_meancov

subroutine SMP_Time_count(lambdaewma,lambda_0,max_time,min_time,ave_time,std_time)
use IMSL
use constants
use Chartingmodule
implicit none
    integer i,j
    integer sample(dim,n)
    real lambda_0(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu_time)
    real Test,lambdaewma
    real dClkStart, dClkFinish
    real max_time,min_time,ave_time,std_time,Tem_time(simu_time),st
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 

    call Exp_meancov(lambda_0,miu_exp,Cov_exp)
    Tem_time=0.0
    Inverse_Cov_exp= .i. Cov_exp
    i=1
    do while(i<=simu_time)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    do j = 1, ewma, 1 
    call sample_generation(lambda_0,sample)
    call MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcov(dim,n,miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    end do 
    call sample_generation(lambda_0,sample)
    call cpu_time(dClkStart)    
    call MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    call Statistic_meanandcov(dim,n,miu_ewma,miu_exp,Cov_ewma,Inverse_Cov_exp,Test)
    call cpu_time(dClkFinish)
    Tem_time(i)=dClkFinish-dClkStart
    print*,i
    i=i+1
    end do
    ave_time=(sum(Tem_time))/(1.0*(simu_time))
    max_time=maxval(Tem_time)
    min_time=minval(Tem_time)
   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu_time
    st=st+(Tem_time(j)-ave_time)**2.0
  enddo 
  std_time=sqrt(st/(simu_time-1.0))
  std_time=std_time/(sqrt(dble(simu_time))) 
  return 
end subroutine  SMP_Time_count