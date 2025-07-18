module constants
    implicit none
    integer, parameter::dim=3,simu_IC=10000,simu_OC=2000,n=1,window_size=50,ewma=50
end module constants
program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real limit,arl,std,lambdaewma,CL,limitL,limitR,lambda_0(dim,dim)
    real dClkStart, dClkFinish
    real c_parameter
    open(10,file="./results/D_ind.txt")
    call cpu_time(dClkStart)
    lambdaewma=0.10
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============c_parameter is 0.5========================="
    c_parameter=0.5 
    !limit=  0.257        
    !arl= 196.340 
    limitL=1.0
	limitR=1.3
    call parameter_setting(c_parameter,lambda_0)
    call D_CLget(dim,lambda_0,CL)
	call D_limitsearch(lambdaewma,lambda_0,CL,limitL,limitR,Limit,arl,std)
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl     
    call D_Poi(lambdaewma,0.10,0.20,limit,lambda_0,CL) 
    call cpu_time(dClkFinish)  
    write(10,*) "Time is", dClkFinish-dClkStart, "seconds" 
    stop
end program main

subroutine D_Poi(lambdaewma,increment_meanshift,increment_covshift,limit,lambda_0,CL)
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real lambdaewma,increment_meanshift,increment_covshift
    real limit,arl,std,CL
    real lambda_0(dim,dim),shift(dim,dim),c_parameter 
    write(10,*),"++++++++++++++++shifts in miu_1+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 9, 1
        shift(1,1)=(lambda_0(1,1)+lambda_0(1,2)+lambda_0(1,3))*increment_meanshift*i
        call D_OCARL(lambdaewma,lambda_0,CL,Limit,shift,arl,std)
        write(10,*),"shift is",increment_meanshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
    write(10,*),"++++++++++++++++shifts in theta_12+++++++++++++++++++++++"
    shift=0.0        
    do i = 1, 9, 1
        shift(1,2)=lambda_0(1,2)*increment_covshift*i
        shift(2,1)=shift(1,2)
        call D_OCARL(lambdaewma,lambda_0,CL,Limit,shift,arl,std)
        write(10,*),"shift is",increment_covshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do    
end subroutine D_Poi

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

subroutine D_OCARL(lambdaewma,lambda_0,CL,Limit,shift,arl,std)
  use constants
  use IMSL
    implicit none
    integer i,j
    real lambda_0(dim,dim),shift(dim,dim),lambda(dim,dim),lambda_OC(dim,dim),lambdaewma
    real simuarl(simu_OC)
    real Limit,Test,arl,st,std,CL
    integer sample(dim,n)
    real sample_ewma(dim,n)

    lambda_OC=lambda_0+shift

    i=1   
    simuarl=0.0
    do while(i<=simu_OC)
    Test=CL    
    do j = 1, ewma, 1
      call sample_generation(lambda_0,sample)  
      call D_statistic(dim,n,lambdaewma,sample,Test)
      if ( Test>CL+Limit .or. Test<CL-Limit ) then
          exit
      end if
    end do
    if ( Test>CL+Limit .or. Test<CL-Limit  ) then
        cycle
    end if
    !print*,"OK"
    do while (Test<CL+Limit .and. Test>CL-Limit)
      call sample_generation(lambda_OC,sample)
      call D_statistic(dim,n,lambdaewma,sample,Test)
      simuarl(i)=simuarl(i)+1.0
    end do 
    print*,i
    print*,simuarl(i)
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
end subroutine D_OCARL

subroutine D_limitsearch(lambdaewma,lambda_0,CL,limitL,limitR,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    real lambda_0(dim,dim),lambda_IC(dim,dim),lambdaewma
    real simuarl(simu_IC)
    real Limit,Test,arl,st,std,LimitL,LimitR,CL
    integer sample(dim,n)
    arl=0.0
  do while (abs(arl-200.0) >= 2.0)
  Limit = (LimitL + LimitR) / 2.0
    i=1  
    simuarl=0.0
    do while(i<=simu_IC)
    Test=CL    
    do j = 1, ewma, 1
      call sample_generation(lambda_0,sample)
      call D_statistic(dim,n,lambdaewma,sample,Test) 
      if ( Test>CL+Limit .or. Test<CL-Limit  ) then
          exit
      end if
      
    end do
    if ( Test>CL+Limit .or. Test<CL-Limit  ) then
        cycle
    end if
    do while (Test<CL+Limit .and. Test>CL-Limit)
      call sample_generation(lambda_0,sample)
      call D_statistic(dim,n,lambdaewma,sample,Test) 
      simuarl(i)=simuarl(i)+1.0
    end do
    print*,i
    print*,simuarl(i)
    print*,"****************" 
    i=i+1
    end do
    arl=(sum(simuarl))/(1.0*(simu_IC))
    write(10,*),"***********************"
    write(10,*),"Limit is",limit
    write(10,*),"IC ARL is", arl
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
end subroutine D_limitsearch