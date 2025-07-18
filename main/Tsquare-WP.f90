module constants
    implicit none
    integer, parameter::dim=3,simu_IC=10000,simu_OC=2000,n=20,window_size=50,ewma=50,level=25,IC_size=20000
end module constants
program main
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    real limit,arl,std,limitl,limitR,lambdaewma
    real*8 dClkStart, dClkFinish
 
    real nu_matrix(dim,dim),theta_matrix(dim,dim),pro_tensor_IC(dim,dim,level),pro_tensor_OC(dim,dim,level),Cov_exp(dim,dim),miu_exp(dim)
    integer i,j,d,sample(dim,IC_size)
    open(10,file="./results/Tsquare_WP.txt")
    do i = 1, dim, 1
        do  j=1,dim,1
            nu_matrix(i,j)=1.0/(1.0+abs(i-j))
            theta_matrix(i,j)=0.8/(1.0+abs(i-j))
        end do
    end do
    call cpu_time(dClkStart)
    lambdaewma=0.10
    write(10,*),"lambdaewma is",lambdaewma
    write(10,*),"===============nu_matrix is ========================="
    write(10,*),nu_matrix
    write(10,*),"===============theta_matrix is ========================="
    write(10,*),theta_matrix
    call reparameterization(nu_matrix,theta_matrix,pro_tensor_IC)
    call samplegenerationWP(IC_size,pro_tensor_IC,sample)
    Cov_exp=theta_matrix
    do i = 1, dim, 1
        do j = 1, dim, 1
            miu_exp(i)=miu_exp(i)+nu_matrix(i,j)
            if (j/=i)then
            Cov_exp(i,i)=Cov_exp(i,i)+Cov_exp(i,j)
            end if          
        end do
    end do   
    
    write(10,*),"===============Poi_matrix_estimated is ========================="
    write(10,*),Cov_exp
    limitl=0.4
    limitR=1.0
    call Tsquare_limitsearch(lambdaewma,pro_tensor_IC,miu_exp,Cov_exp,limitl,limitR,Limit,arl,std)  
    write(10,*),"Limit is",Limit
    write(10,*),"IC_arl is",arl
    write(10,*),"IC_std is",std
    call Tsquare_Poi(lambdaewma,nu_matrix,theta_matrix,0.05,-0.10,limit,miu_exp,Cov_exp)
    call cpu_time(dClkFinish)  
    write(10,*) "Time is", dClkFinish-dClkStart, "seconds" 
    stop
end program main

subroutine Tsquare_Poi(lambdaewma,nu_matrix,theta_matrix,increment_meanshift,increment_covshift,limit,miu_exp,Cov_exp)
    use IMSL
    use constants
    use Chartingmodule
    implicit none
    integer i,j
    real lambdaewma,increment_meanshift,increment_covshift
    real limit,arl,std
    real pro_tensor_IC(dim,dim,level),pro_tensor_OC(dim,dim,level),nu_matrix(dim,dim),theta_matrix(dim,dim),miu(dim,dim),Cov(dim,dim),Cov_exp(dim,dim),miu_exp(dim)
    call reparameterization(nu_matrix,theta_matrix,pro_tensor_IC)

    write(10,*),"++++++++++++++++shifts in miu_11+++++++++++++++++++++++"
    miu=nu_matrix
    Cov=theta_matrix   
    do i = 1, 9, 1
        !if (i==5) then
           !cycle
        !endif
        miu(1,1)=nu_matrix(1,1)*(1+increment_meanshift*i)
        call reparameterization(miu,Cov,pro_tensor_OC)
        call Tsquare_OCARL(lambdaewma,pro_tensor_IC,pro_tensor_OC,miu_exp,Cov_exp,Limit,arl,std)
        write(10,*),"meanshift in 11 is",increment_meanshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
    write(10,*),"++++++++++++++++shifts in miu_12+++++++++++++++++++++++"
    miu=nu_matrix
    Cov=theta_matrix   
    do i = 1, 9, 1
        !if (i==5) then
           !cycle
        !endif
        miu(1,2)=nu_matrix(1,2)*(1+increment_meanshift*i)
        call reparameterization(miu,Cov,pro_tensor_OC)
        call Tsquare_OCARL(lambdaewma,pro_tensor_IC,pro_tensor_OC,miu_exp,Cov_exp,Limit,arl,std)
        write(10,*),"meanshift in 11 is",increment_meanshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
        write(10,*),"++++++++++++++++shifts in theta_11+++++++++++++++++++++++"
    miu=nu_matrix
    Cov=theta_matrix        
    do i = 1, 9, 1
        !if (i==5) then
            !cycle
        !endif
        Cov(1,1)=theta_matrix(1,1)*(1+increment_covshift*i)
        call reparameterization(miu,Cov,pro_tensor_OC)
        call Tsquare_OCARL(lambdaewma,pro_tensor_IC,pro_tensor_OC,miu_exp,Cov_exp,Limit,arl,std)
        write(10,*),"Covshift in 11 is",increment_covshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do
    write(10,*),"++++++++++++++++shifts in theta_12+++++++++++++++++++++++"
    miu=nu_matrix
    Cov=theta_matrix        
    do i = 1, 9, 1
        !if (i==5) then
            !cycle
        !endif
        Cov(1,2)=theta_matrix(1,2)*(1+increment_covshift*i)
        call reparameterization(miu,Cov,pro_tensor_OC)
        call Tsquare_OCARL(lambdaewma,pro_tensor_IC,pro_tensor_OC,miu_exp,Cov_exp,Limit,arl,std)
        write(10,*),"Covshift in 12  is",increment_covshift*i
        write(10,*),"arl is",arl
        write(10,*),"std is",std
        write(10,*),"--------------------------------------"
    end do    
end subroutine Tsquare_Poi

subroutine WP_prob(omiga,a,pro)
    use IMSL
    use constants
    implicit none   
    integer i
    real omiga,a,pro(level),sum_tem
    sum_tem=0
    pro=0.0
    do i = 1, level, 1
        pro(i)=max(0.0,(omiga**(i-1))*(i-1+a)*exp(-omiga)/(FAC(i-1)*(omiga+a)))
        sum_tem=sum_tem+pro(i)
        if ( sum_tem>1.0 ) then
            pro(i)=max(0.0,pro(i)-(sum_tem-1.0))
            exit
        end if
    end do 
    return
end subroutine WP_prob

subroutine reparameterization(nu_matrix,theta_matrix,pro_tensor)
    use IMSL
    use constants
    implicit none
    integer  i,j,type_dis
    real nu_matrix(dim,dim),theta_matrix(dim,dim),omiga_matrix(dim,dim),beta_matrix(dim,dim),pro_tensor(dim,dim,level)

    do i = 1, dim, 1
        do j = i, dim, 1
            omiga_matrix(i,j)=nu_matrix(i,j)-sqrt(nu_matrix(i,j)-theta_matrix(i,j))
            beta_matrix(i,j)=omiga_matrix(i,j)*(1.0-sqrt(nu_matrix(i,j)-theta_matrix(i,j)))/sqrt(nu_matrix(i,j)-theta_matrix(i,j))
            call WP_prob(omiga_matrix(i,j),beta_matrix(i,j),pro_tensor(i,j,:))  
        end do  
    end do
    return
end subroutine reparameterization

subroutine samplegenerationWP(sample_size,pro_tensor,sample)
    use IMSL
    use constants
    implicit none
    integer  sample_size
    integer sample(dim,sample_size),samplegeneration(dim,dim,sample_size)
    integer i,j,k
    integer  IWK(level)
    real PROBS(level), WK(level),pro_tensor(dim,dim,level)  
    samplegeneration=0
    do i = 1, dim, 1
            do j = i, dim, 1
                CALL RNGDA (sample_size, 0, 0, level, pro_tensor(i,j,:), IWK, WK, samplegeneration(i,j,:))
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
    
end subroutine samplegenerationWP

subroutine Tsquare_OCARL(lambdaewma,pro_tensor_IC,pro_tensor_OC,miu_exp,Cov_exp,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real pro_tensor_IC(dim,dim,level),pro_tensor_OC(dim,dim,level)
    real Cov_exp(dim,dim),Inverse_Cov_exp(dim,dim)
    real simuarl(simu_OC)
    real Limit,Test,arl,st,std,lambdaewma
    real Cov_ewma(dim,dim),miu_ewma(dim),miu_exp(dim) 

    Inverse_Cov_exp= .i. Cov_exp
    i=1
    do while(i<=simu_OC)
    Test=0.0
    Cov_ewma=Cov_exp
    miu_ewma=miu_exp
    !print*,Cov_ewma
    do j = 1, ewma, 1 
    call samplegenerationWP(n,pro_tensor_IC,sample)
    call Mean_EstimationI(dim,n,lambdaewma,sample,miu_ewma)
    call T_square(dim,n,miu_ewma,miu_exp,Inverse_Cov_exp,Test)
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
    call samplegenerationWP(n,pro_tensor_OC,sample)
    call Mean_EstimationI(dim,n,lambdaewma,sample,miu_ewma)
    call T_square(dim,n,miu_ewma,miu_exp,Inverse_Cov_exp,Test)
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
end subroutine Tsquare_OCARL

subroutine Tsquare_limitsearch(lambdaewma,pro_tensor_IC,miu_exp,Cov_exp,limitl,limitR,Limit,arl,std)
use IMSL
use constants
implicit none
    integer i,j
    integer sample(dim,n)
    real pro_tensor_IC(dim,dim,level)
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
    call samplegenerationWP(n,pro_tensor_IC,sample)
    call Mean_EstimationI(dim,n,lambdaewma,sample,miu_ewma)
    call T_square(dim,n,miu_ewma,miu_exp,Inverse_Cov_exp,Test)
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
    call samplegenerationWP(n,pro_tensor_IC,sample)
    call Mean_EstimationI(dim,n,lambdaewma,sample,miu_ewma)
    call T_square(dim,n,miu_ewma,miu_exp,Inverse_Cov_exp,Test)
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
end subroutine Tsquare_limitsearch