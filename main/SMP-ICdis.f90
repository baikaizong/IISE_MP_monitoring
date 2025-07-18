module constants
    implicit none
    integer, parameter::dim=3,n=20,window_size=50,simu_IC=10000,ewma=50
end module constants
program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real limit,arl,limitR,Limitl,std,sdrl,qtlone,qtlfive,qtlnine,FAR,lambdaewma,qt(5)
    real*8 dClkStart, dClkFinish
    real c_parameter,lambda_0(dim,dim),miu_exp(dim),Cov_exp(dim,dim)
    open(10,file="./results/AMP_Poi-IC-0.5.txt")
	lambdaewma=0.10
    call cpu_time(dClkStart)
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
    limit=  0.826        
	call parameter_setting(c_parameter,lambda_0)
    call Exp_meancov(lambda_0,miu_exp,Cov_exp)
	call SMP_IC(lambdaewma,lambda_0,miu_exp,Cov_exp,Limit,arl,std,sdrl,qt,FAR)   
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl 
	write(10,*),"std is",std
    write(10,*),"sdrl is",sdrl
	write(10,*),"qt0.10 is",qt(1)
    write(10,*),"qt0.25 is",qt(2)
	write(10,*),"qt0.50 is",qt(3)
    write(10,*),"qt0.75 is",qt(4)
    write(10,*),"qt0.90 is",qt(5)
    write(10,*),"FAR is",FAR
   
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
    return
end subroutine Exp_meancov

subroutine SMP_IC(lambdaewma,lambda_0,miu_exp,Cov_exp,Limit,arl,std,sdrl,qt,FAR )
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real lambda_0(dim,dim)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu_IC)
    real Limit,Test,arl,st,lambdaewma,std,sdrl,qt(5),FAR 
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 
    

    Inverse_Cov_exp= .i. Cov_exp
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

   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu_IC
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  sdrl=sqrt(st/(simu_IC))
  std=sdrl/(sqrt(dble(simu_IC)))
  CALL SVRGN (simu_IC, simuarl, simuarl)
  qt(1)=(simuarl(simu_IC*0.10)+simuarl(simu_IC*0.10+1))/2.0
  qt(2)=(simuarl(simu_IC*0.25)+simuarl(simu_IC*0.25+1))/2.0
  qt(3)=(simuarl(simu_IC*0.50)+simuarl(simu_IC*0.50+1))/2.0
  qt(4)=(simuarl(simu_IC*0.75)+simuarl(simu_IC*0.75+1))/2.0
  qt(5)=(simuarl(simu_IC*0.90)+simuarl(simu_IC*0.90+1))/2.0
  FAR=0
  do i = 1, simu_IC, 1
      if ( simuarl(i)>=30 ) then
          FAR=1.0*i/simu_IC
          exit
      end if
  end do
  return 
end subroutine SMP_IC