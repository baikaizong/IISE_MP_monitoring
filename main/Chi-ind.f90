module constants
    implicit none
    integer, parameter::dim=3,simu_IC=10000,simu_OC=2000,n=1,window_size=50,ewma=50,IC_sample_size=20000
end module constants
program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    integer sampleIC(dim,IC_sample_size)
    real limit,arl,std,limitR,limitl,lambdaewma
    real dClkStart, dClkFinish
    real c_parameter
	real lambda_0(dim,dim),shift(dim,dim) 
    real midian(dim),ICpr_cell(2**dim)
    open(10,file="./results/Pei_ind.txt")
    call cpu_time(dClkStart)
    lambdaewma=0.10
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
	limitR=1.3
	limitl=0.7
    !limit=  0.950        
    !arl=  203.500 
    call parameter_setting(c_parameter,lambda_0)
    call sample_generation(IC_sample_size,lambda_0,sampleIC)
    call Median_get(dim,IC_sample_size,sampleIC,midian,ICpr_cell)  
	call Pei_limitsearch(lambdaewma,lambda_0,midian,ICpr_cell,limitL,limitR,Limit,arl,std)
	write(10,*),"midian is",midian  
	write(10,*),"ICpr_cell is",ICpr_cell
    write(10,*),"limit is",limit
    write(10,*),"IC_arl is",arl   
	write(10,*),"IC_std is",std   
    call Pei_Poi(lambdaewma,0.10,0.20,limit,c_parameter,midian,ICpr_cell) 
    call cpu_time(dClkFinish)  
    write(10,*) "Time is", dClkFinish-dClkStart, "seconds" 
    stop
end program main

subroutine Pei_Poi(lambdaewma,increment_meanshift,increment_covshift,limit,c_parameter,midian,ICpr_cell)
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real increment_meanshift,increment_covshift
    real limit,arl,std,lambdaewma
    real lambda_0(dim,dim),shift(dim,dim),c_parameter 
    real midian(dim),ICpr_cell(2**dim)
    call parameter_setting(c_parameter,lambda_0)
    write(10,*),"++++++++++++++++shifts in miu_1+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 10, 1
        shift(1,1)=(lambda_0(1,1)+lambda_0(1,2)+lambda_0(1,3))*increment_meanshift*i
        call Pei_OCARL(lambdaewma,lambda_0,midian,ICpr_cell,Limit,shift,arl,std)
        write(10,*),"shift is",increment_meanshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
    write(10,*),"++++++++++++++++shifts in theta_12+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 10, 1
        shift(1,2)=lambda_0(1,2)*increment_covshift*i
        shift(2,1)=shift(1,2)
        call Pei_OCARL(lambdaewma,lambda_0,midian,ICpr_cell,Limit,shift,arl,std)
        write(10,*),"shift is",increment_covshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do    
end subroutine Pei_Poi

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

subroutine sample_generation(samplesize,lambda_0,sample)
    use IMSL
    use constants
    implicit none
    integer samplesize
    real lambda_0(dim,dim)
    integer sample(dim,samplesize),samplegeneration(dim,dim,samplesize)
    integer i,j,k
    samplegeneration=0
    do i = 1, dim, 1
            do j = i, dim, 1
                if ( lambda_0(i,j)/=0.0 ) then
                    call RNPOI(samplesize,lambda_0(i,j),samplegeneration(i,j,:))
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

subroutine Pei_OCARL(lambdaewma,lambda_0,midian,ICpr_cell,Limit,shift,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real lambda(dim,dim),lambda_0(dim,dim),shift(dim,dim),lambdaewma
    real simuarl(simu_OC)
    real Limit,Test,arl,st,std
    real midian(dim),ICpr_cell(2**dim),sample_cutting_ewma(2**dim) 
    lambda=lambda_0+shift
    i=1
    do while(i<=simu_OC)
    Test=0.0
    sample_cutting_ewma=1.0*n*ICpr_cell
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call sample_generation(n,lambda_0,sample)
    call Statistic_Pei(dim,n,lambdaewma,sample,ICpr_cell,midian,sample_cutting_ewma,test)
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
    call sample_generation(n,lambda,sample)
    call Statistic_Pei(dim,n,lambdaewma,sample,ICpr_cell,midian,sample_cutting_ewma,test)
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
end subroutine Pei_OCARL

subroutine Pei_limitsearch(lambdaewma,lambda_0,midian,ICpr_cell,limitL,limitR,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
  integer sample(dim,n)
    real lambda_0(dim,dim),sample_cutting_ewma(2**dim),lambdaewma
    real simuarl(simu_IC)
    real Limit,Test,arl,st,std,limitl,limitR
    real midian(dim),ICpr_cell(2**dim)
        arl=0.0
    do while (abs(arl-200.0) >= 2.0)
    Limit = (LimitL + LimitR) / 2.0
    i=1
     do while(i<simu_IC)
    Test=0.0
    sample_cutting_ewma=1.0*n*ICpr_cell
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call sample_generation(n,lambda_0,sample)
    call Statistic_Pei(dim,n,lambdaewma,sample,ICpr_cell,midian,sample_cutting_ewma,test)
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
    call sample_generation(n,lambda_0,sample)
    call Statistic_Pei(dim,n,lambdaewma,sample,ICpr_cell,midian,sample_cutting_ewma,test)
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
end subroutine Pei_limitsearch