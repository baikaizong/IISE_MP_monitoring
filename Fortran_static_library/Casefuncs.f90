module funcsconstants
    implicit none
    real,parameter:: epsilon=0.00000001
end module funcsconstants

subroutine Poisample_generation(lambda_IC,n,sample)
    use IMSL
    use funcsconstants
    implicit none
    integer n
    real*8 lambda_IC(4,4)
    real lambda_0(4,4)
    integer sample(4,n),samplegeneration(4,4,n)
    integer i,j,k
    samplegeneration=0
    lambda_0=lambda_IC
    do i = 1, 4, 1
            do j = i, 4, 1
                if ( lambda_0(i,j)/=0.0 ) then
                    call RNPOI(n,lambda_0(i,j),samplegeneration(i,j,:))
                end if
                samplegeneration(j,i,:)=samplegeneration(i,j,:)          
            end do        
    end do
    sample=0.0
    do i = 1, 4, 1
        do j = 1, 4, 1
        sample(i,:)= sample(i,:)+samplegeneration(i,j,:)          
        end do
    end do
    return  
end subroutine Poisample_generation

subroutine MP_Estimation(allsample_size,allsample,index_leave,lambda_em,LRTvalue)
    use IMSL
    use funcsconstants
    implicit none
    integer allsample_size,sample_size,index_leave
    integer allsample(4,allsample_size),sample(4,allsample_size-1),max_val(4),max_dim(4),tem_index_leave(4)
    real*8 lambda_0(4,4),lambda_last(4,4),LRTvalue,miu_0(4),lambda_em(4,4)
    integer i,j,k,d
    integer iterations
    real*8 FT,likelihoodlast,likelihood,likelihoodmax,sum_meanIC
    real*8 S(4,4)
    real*8, allocatable :: pro_table(:,:,:,:)
    lambda_0=0.1
    sample_size=allsample_size-1
    sample(:,:)=allsample(:,:)
    if ( index_leave<allsample_size ) then
        sample(:,index_leave:sample_size)=allsample(:,index_leave+1:allsample_size)      
    end if

    do i = 1, 4, 1
        max_dim(i)=maxval(allsample(i,:))
    end do
    allocate(pro_table(max_dim(1)+2,max_dim(2)+2,max_dim(3)+2,max_dim(4)+2))
    do i = 1, 4, 1
        miu_0(i)=1.0*sum(sample(i,1:sample_size))/sample_size
        lambda_0(i,i)=miu_0(i)  
    end do
    do i = 1, 4, 1
        max_val(i)=min(6.0*lambda_0(i,i),1.0*max_dim(i),34.0)
    end do
    FT=1.0
    likelihoodlast=-10000000000000000000.0
    iterations=0
    likelihoodmax=-100000000000000000000.0
    !print*,lambda_0
    do while ( FT>0.000001)
        pro_table(:,:,:,:)=0.0
        !print*,"---------2-----------"       
        pro_table(2,2,2,2)=exp(-(lambda_0(1,1)+lambda_0(2,2)+lambda_0(3,3)+lambda_0(4,4)+lambda_0(1,2)+lambda_0(1,3)+lambda_0(1,4)+lambda_0(2,3)+lambda_0(2,4)+lambda_0(3,4))) 
        do i = 3, max_val(1)+2, 1
            pro_table(i,2,2,2)= pro_table(2,2,2,2)*(lambda_0(1,1)**(i-2))/FAC(i-2)
        end do
        do i = 3, max_val(2)+2, 1
            pro_table(2,i,2,2)= pro_table(2,2,2,2)*(lambda_0(2,2)**(i-2))/FAC(i-2)
        end do
        do i = 3, max_val(3)+2, 1
            pro_table(2,2,i,2)= pro_table(2,2,2,2)*(lambda_0(3,3)**(i-2))/FAC(i-2)
        end do 
        do i = 3, max_val(4)+2, 1
            pro_table(2,2,2,i)= pro_table(2,2,2,2)*(lambda_0(4,4)**(i-2))/FAC(i-2)
        end do 
        !print*,"---------3-----------"
        !print*,sum(pro_table)
        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(2)+2, 1
                pro_table(i,j,2,2)=(lambda_0(1,1)*pro_table(i-1,j,2,2)+lambda_0(1,2)*pro_table(i-1,j-1,2,2))/(i-2)
            end do  
        end do 

        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(3)+2, 1
                pro_table(i,2,j,2)=(lambda_0(1,1)*pro_table(i-1,2,j,2)+lambda_0(1,3)*pro_table(i-1,2,j-1,2))/(i-2)
            end do  
        end do 
        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(4)+2, 1
                pro_table(i,2,2,j)=(lambda_0(1,1)*pro_table(i-1,2,2,j)+lambda_0(1,4)*pro_table(i-1,2,2,j-1))/(i-2)
            end do  
        end do   
        do i = 3, max_val(2)+2, 1
            do j = 3, max_val(3)+2, 1
                pro_table(2,i,j,2)=(lambda_0(2,2)*pro_table(2,i-1,j,2)+lambda_0(2,3)*pro_table(2,i-1,j-1,2))/(i-2)
            end do  
        end do 
        do i = 3, max_val(2)+2, 1
            do j = 3, max_val(4)+2, 1
                pro_table(2,i,2,j)=(lambda_0(2,2)*pro_table(2,i-1,2,j)+lambda_0(2,4)*pro_table(2,i-1,2,j-1))/(i-2)
            end do  
        end do 
        do i = 3, max_val(3)+2, 1
            do j = 3, max_val(4)+2, 1
                pro_table(2,2,i,j)=(lambda_0(3,3)*pro_table(2,2,i-1,j)+lambda_0(3,4)*pro_table(2,2,i-1,j-1))/(i-2)
            end do  
        end do 
        !print*,"---------4-----------"
        !print*,sum(pro_table)
        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(2)+2, 1
                do k = 3, max_val(3)+2, 1
                   pro_table(i,j,k,2)=((lambda_0(1,1))*pro_table(i-1,j,k,2)+lambda_0(1,2)*pro_table(i-1,j-1,k,2)+lambda_0(1,3)*pro_table(i-1,j,k-1,2))/(i-2)
                end do  
            end do    
        end do
        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(2)+2, 1
                do k = 3, max_val(4)+2, 1
                   pro_table(i,j,2,k)=((lambda_0(1,1))*pro_table(i-1,j,2,k)+lambda_0(1,2)*pro_table(i-1,j-1,2,k)+lambda_0(1,4)*pro_table(i-1,j,2,k-1))/(i-2)
                end do  
            end do    
        end do
        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(3)+2, 1
                do k = 3, max_val(4)+2, 1
                   pro_table(i,2,j,k)=((lambda_0(1,1))*pro_table(i-1,2,j,k)+lambda_0(1,3)*pro_table(i-1,2,j-1,k)+lambda_0(1,4)*pro_table(i-1,2,j,k-1))/(i-2)
                end do  
            end do    
        end do
        do i = 3, max_val(2)+2, 1
            do j = 3, max_val(3)+2, 1
                do k = 3, max_val(4)+2, 1
                   pro_table(2,i,j,k)=((lambda_0(2,2))*pro_table(2,i-1,j,k)+lambda_0(2,3)*pro_table(2,i-1,j-1,k)+lambda_0(2,4)*pro_table(2,i-1,j,k-1))/(i-2)
                end do  
            end do    
        end do
        !print*,"---------5-----------"
        !print*,sum(pro_table)
        do i = 3, max_val(1)+2, 1
            do j = 3, max_val(2)+2, 1
                do k = 3, max_val(3)+2, 1
                    do d = 3, max_val(4)+2, 1
                   pro_table(i,j,k,d)=((lambda_0(1,1))*pro_table(i-1,j,k,d)+lambda_0(1,2)*pro_table(i-1,j-1,k,d)+lambda_0(1,3)*pro_table(i-1,j,k-1,d)+lambda_0(1,4)*pro_table(i-1,j,k,d-1))/(i-2)
                   end do
                end do  
            end do    
        end do
        do i = 1, max_dim(1)+2, 1
            do j = 1, max_dim(2)+2, 1
                do k = 1, max_dim(3)+2, 1
                    do d = 1, max_dim(4)+2, 1
                        if ( pro_table(i,j,k,d)<=0.000000000000001 ) then
                            pro_table(i,j,k,d)=0.000000000000001
                        end if                    
                    end do
                end do  
            end do    
        end do
        !print*,"---------3-----------"
        !print*,sum(pro_table)
        likelihood=0.0
        do i = 1, sample_size, 1
           likelihood=likelihood+log(pro_table(sample(1,i)+2,sample(2,i)+2,sample(3,i)+2,sample(4,i)+2))
        end do
        S=0.0
        do i = 1, sample_size, 1     
            do k = 1, 3, 1
                do d = k+1, 4, 1
                    tem_index_leave=sample(:,i)+2
                    tem_index_leave(k)=tem_index_leave(k)-1
                    tem_index_leave(d)=tem_index_leave(d)-1
                    S(k,d)=S(k,d)+(lambda_0(k,d))*pro_table(tem_index_leave(1),tem_index_leave(2),tem_index_leave(3),tem_index_leave(4))/pro_table(sample(1,i)+2,sample(2,i)+2,sample(3,i)+2,sample(4,i)+2)
                end do    
            end do
        end do
        FT=abs(likelihood-likelihoodlast)
        if ( likelihood-likelihoodmax>0.0 ) then
          lambda_em = lambda_0
          likelihoodmax=likelihood
          !exit
        end if 
        lambda_last=lambda_0
        do i = 1, 3, 1
            do j = i+1, 4, 1
               lambda_0(i,j)=min(S(i,j)/sample_size,miu_0(i),miu_0(j))
               lambda_0(j,i)=lambda_0(i,j)
             end do    
        end do 
           
        do i = 1, 4, 1
            do j = 1, 4, 1
            if (j/=i)then
                lambda_0(i,i)=miu_0(i)-lambda_0(i,j)
            end if          
            end do
            lambda_0(i,i)=max(0.0,lambda_0(i,i))
        end do        
        !print*,"---------5-----------"
        !print*,lambda_0
     
        likelihoodlast=likelihood
        iterations=iterations+1
        if ( iterations>50 ) then
            exit
        end if
    end do
    sum_meanIC=lambda_em(1,1)+lambda_em(2,2)+lambda_em(3,3)+lambda_em(4,4)+lambda_em(1,2)+lambda_em(1,3)+lambda_em(1,4)+lambda_em(2,3)+lambda_em(2,4)+lambda_em(3,4)
    if ( sum(allsample(:,index_leave))<sum_meanIC ) then
        LRTvalue=0.0
    else
        LRTvalue=-log(pro_table(allsample(1,index_leave)+2,allsample(2,index_leave)+2,allsample(3,index_leave)+2,allsample(4,index_leave)+2))
    end if
    return        
end subroutine MP_Estimation

subroutine LRT_test(lambdaewma,window_size,n,sample,lambda_IC,globaltest)
    use funcsconstants
    use IMSL
    implicit none
    integer window_size,n
    integer sample(window_size,4,n),max_dim(4),max_val(4),tem_index_leave(4)
    real*8 lambda_IC(4,4),lambda_0(4,4),lambda_last(4,4),tem_frac,miu_tem(4),lambdaewma,lambdaco,miu_0(4),lambda_em(4,4)
    real*8 globaltest,sum_meanIC
    integer i,j,k,d
    integer iterations
    real*8 FT,likelihoodlast,likelihood,likelihoodnull,likelihoodmax
    real*8 likelihood_array,S(4,4),S_ewma(4,4)
    real*8, allocatable :: pro_table(:,:,:,:)
    integer tmp(4,window_size)
    lambdaco=1.0-lambdaewma
    lambda_0=lambda_IC
    do i = 1, 4, 1
        do j = 1, window_size, 1
            tmp(i,j)=maxval(sample(j,i,:))
        end do
        max_dim(i)=maxval(tmp(i,:))
    end do    
    do i = 1, 4, 1
        miu_0(i)=1.0*sum(sample(1,i,:))/n
    end do  
    do i = 1, 4, 1
        do j = 2, window_size, 1
            miu_0(i)=lambdaco*miu_0(i)+lambdaewma*sum(sample(j,i,:))/n   
        end do        
    end do
    sum_meanIC=lambda_0(1,1)+lambda_0(2,2)+lambda_0(3,3)+lambda_0(4,4)+lambda_0(1,2)+lambda_0(1,3)+lambda_0(1,4)+lambda_0(2,3)+lambda_0(2,4)+lambda_0(3,4)
    !if ( sum(miu_0)<sum_meanIC ) then
        !globaltest=0
        !return
    !end if
    !print*,"*****one*************"
    allocate(pro_table(max_dim(1)+2,max_dim(2)+2,max_dim(3)+2,max_dim(4)+2))

    !print*,"---------1-----------"
    do i = 1, 4, 1
        max_val(i)=min(6.0*miu_0(i),1.0*max_dim(i),34.0)
    end do    
    pro_table=0.0
    
  
    pro_table(2,2,2,2)=exp(-sum_meanIC) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2,2)= pro_table(2,2,2,2)*(lambda_0(1,1)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2,2)= pro_table(2,2,2,2)*(lambda_0(2,2)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i,2)= pro_table(2,2,2,2)*(lambda_0(3,3)**(i-2))/FAC(i-2)
    end do 
    do i = 3, max_val(4)+2, 1
        pro_table(2,2,2,i)= pro_table(2,2,2,2)*(lambda_0(4,4)**(i-2))/FAC(i-2)
    end do 
     ! print*,"---------2-----------"
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2,2)=(lambda_0(1,1)*pro_table(i-1,j,2,2)+lambda_0(1,2)*pro_table(i-1,j-1,2,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j,2)=(lambda_0(1,1)*pro_table(i-1,2,j,2)+lambda_0(1,3)*pro_table(i-1,2,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(i,2,2,j)=(lambda_0(1,1)*pro_table(i-1,2,2,j)+lambda_0(1,4)*pro_table(i-1,2,2,j-1))/(i-2)
        end do  
    end do   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j,2)=(lambda_0(2,2)*pro_table(2,i-1,j,2)+lambda_0(2,3)*pro_table(2,i-1,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,i,2,j)=(lambda_0(2,2)*pro_table(2,i-1,2,j)+lambda_0(2,4)*pro_table(2,i-1,2,j-1))/(i-2)
        end do  
    end do 
    do i = 3, max_val(3)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,2,i,j)=(lambda_0(3,3)*pro_table(2,2,i-1,j)+lambda_0(3,4)*pro_table(2,2,i-1,j-1))/(i-2)
        end do  
    end do 
    !print*,"--------------------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k,2)=((lambda_0(1,1))*pro_table(i-1,j,k,2)+lambda_0(1,2)*pro_table(i-1,j-1,k,2)+lambda_0(1,3)*pro_table(i-1,j,k-1,2))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,j,2,k)=((lambda_0(1,1))*pro_table(i-1,j,2,k)+lambda_0(1,2)*pro_table(i-1,j-1,2,k)+lambda_0(1,4)*pro_table(i-1,j,2,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,2,j,k)=((lambda_0(1,1))*pro_table(i-1,2,j,k)+lambda_0(1,3)*pro_table(i-1,2,j-1,k)+lambda_0(1,4)*pro_table(i-1,2,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(2,i,j,k)=((lambda_0(2,2))*pro_table(2,i-1,j,k)+lambda_0(2,3)*pro_table(2,i-1,j-1,k)+lambda_0(2,4)*pro_table(2,i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
                do d = 3, max_val(4)+2, 1
               pro_table(i,j,k,d)=((lambda_0(1,1))*pro_table(i-1,j,k,d)+lambda_0(1,2)*pro_table(i-1,j-1,k,d)+lambda_0(1,3)*pro_table(i-1,j,k-1,d)+lambda_0(1,4)*pro_table(i-1,j,k,d-1))/(i-2)
               end do
            end do  
        end do    
    end do
    do i = 1, max_dim(1)+2, 1
        do j = 1, max_dim(2)+2, 1
            do k = 1, max_dim(3)+2, 1
                do d = 1, max_dim(4)+2, 1
                    if ( pro_table(i,j,k,d)<=0.000000000000001 ) then
                        pro_table(i,j,k,d)=0.000000000000001
                    end if                    
                end do
            end do  
        end do    
    end do
     ! print*,"---------3-----------"
    likelihood=0.0
    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       !print*,pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2)
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2,sample(j,4,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lambdaewma*likelihood_array+lambdaco*likelihood 
    end if         
    end do
    likelihoodnull=likelihood 
    print*,"*****likelihoodnull*************"
    print*,likelihoodnull    
    lambda_em=lambda_0
    miu_tem=miu_0
     do i = 1, 4, 1
        do j = 1, 4, 1
        if (j/=i)then
        miu_tem(i)=miu_tem(i)-lambda_0(i,j)
        end if          
        end do
        lambda_0(i,i)=max(0.0,miu_tem(i))
    end do
    !print*,"*****two*************"
    FT=1.0
    likelihoodlast=likelihood
    iterations=0
    likelihoodmax=likelihood

    do while ( FT>epsilon)
    pro_table=0.0
    
   ! print*,"lambda_0 is", lambda_0
    pro_table(2,2,2,2)=exp(-(lambda_0(1,1)+lambda_0(2,2)+lambda_0(3,3)+lambda_0(4,4)+lambda_0(1,2)+lambda_0(1,3)+lambda_0(1,4)+lambda_0(2,3)+lambda_0(2,4)+lambda_0(3,4))) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2,2)= pro_table(2,2,2,2)*(lambda_0(1,1)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2,2)= pro_table(2,2,2,2)*(lambda_0(2,2)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i,2)= pro_table(2,2,2,2)*(lambda_0(3,3)**(i-2))/FAC(i-2)
    end do 
    do i = 3, max_val(4)+2, 1
        pro_table(2,2,2,i)= pro_table(2,2,2,2)*(lambda_0(4,4)**(i-2))/FAC(i-2)
    end do 
    !print*,"---------1-----------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2,2)=(lambda_0(1,1)*pro_table(i-1,j,2,2)+lambda_0(1,2)*pro_table(i-1,j-1,2,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j,2)=(lambda_0(1,1)*pro_table(i-1,2,j,2)+lambda_0(1,3)*pro_table(i-1,2,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(i,2,2,j)=(lambda_0(1,1)*pro_table(i-1,2,2,j)+lambda_0(1,4)*pro_table(i-1,2,2,j-1))/(i-2)
        end do  
    end do   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j,2)=(lambda_0(2,2)*pro_table(2,i-1,j,2)+lambda_0(2,3)*pro_table(2,i-1,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,i,2,j)=(lambda_0(2,2)*pro_table(2,i-1,2,j)+lambda_0(2,4)*pro_table(2,i-1,2,j-1))/(i-2)
        end do  
    end do 
    do i = 3, max_val(3)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,2,i,j)=(lambda_0(3,3)*pro_table(2,2,i-1,j)+lambda_0(3,4)*pro_table(2,2,i-1,j-1))/(i-2)
        end do  
    end do 
    !print*,"---------2-----------"
    !print*,sum(pro_table)
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k,2)=((lambda_0(1,1))*pro_table(i-1,j,k,2)+lambda_0(1,2)*pro_table(i-1,j-1,k,2)+lambda_0(1,3)*pro_table(i-1,j,k-1,2))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,j,2,k)=((lambda_0(1,1))*pro_table(i-1,j,2,k)+lambda_0(1,2)*pro_table(i-1,j-1,2,k)+lambda_0(1,4)*pro_table(i-1,j,2,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,2,j,k)=((lambda_0(1,1))*pro_table(i-1,2,j,k)+lambda_0(1,3)*pro_table(i-1,2,j-1,k)+lambda_0(1,4)*pro_table(i-1,2,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(2,i,j,k)=((lambda_0(2,2))*pro_table(2,i-1,j,k)+lambda_0(2,3)*pro_table(2,i-1,j-1,k)+lambda_0(2,4)*pro_table(2,i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do
     !print*,"---------3-----------"
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
                do d = 3, max_val(4)+2, 1
               pro_table(i,j,k,d)=((lambda_0(1,1))*pro_table(i-1,j,k,d)+lambda_0(1,2)*pro_table(i-1,j-1,k,d)+lambda_0(1,3)*pro_table(i-1,j,k-1,d)+lambda_0(1,4)*pro_table(i-1,j,k,d-1))/(i-2)
               end do
            end do  
        end do    
    end do
    do i = 1, max_dim(1)+2, 1
        do j = 1, max_dim(2)+2, 1
            do k = 1, max_dim(3)+2, 1
                do d = 1, max_dim(4)+2, 1
                    if ( pro_table(i,j,k,d)<=0.0000000000000001 ) then
                        pro_table(i,j,k,d)=0.000000000000001
                    end if                    
                end do
            end do  
        end do    
    end do    
    !print*,"---------4-----------"
    likelihood=0.0
    do i = 1, n, 1
       likelihood_array=likelihood_array+log(pro_table(sample(1,1,i)+2,sample(1,2,i)+2,sample(1,3,i)+2,sample(1,4,i)+2))
    end do
    do j = 2, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2,sample(j,4,i)+2))
    end do
    likelihood=lambdaewma*likelihood_array+lambdaco*likelihood     
    end do
    !print*,"---------5-----------"
    S_ewma=0.0
    do j = 1, window_size, 1
    S=0.0
    do i = 1, n, 1     
        tem_frac=pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2,sample(j,4,i)+2)
        do k = 1, 4-1, 1
            do d = k+1, 4, 1
                tem_index_leave=sample(j,:,i)+2
                tem_index_leave(k)=tem_index_leave(k)-1
                tem_index_leave(d)=tem_index_leave(d)-1
                S(k,d)=S(k,d)+(lambda_0(k,d))*pro_table(tem_index_leave(1),tem_index_leave(2),tem_index_leave(3),tem_index_leave(4))/tem_frac
            end do    
        end do
    end do
    !print*,"---------6-----------"
    if ( j==1 ) then
    do k = 1, 4-1, 1
        do d = k+1, 4, 1
           S_ewma(k,d)=S(k,d)
         end do    
    end do
    else
    do k = 1, 4-1, 1
        do d = k+1, 4, 1
           S_ewma(k,d)=lambdaewma*S(k,d)+lambdaco*S_ewma(k,d)
         end do    
    end do
    end if
    enddo
    !print*,"---------7-----------"
    lambda_last=lambda_0
    !print*,"lambda_0 is",lambda_0
    do k = 1, 4-1, 1
        do d = k+1, 4, 1
           lambda_0(k,d)=min(S_ewma(k,d)/n,miu_0(k),miu_0(d))
           lambda_0(d,k)=lambda_0(k,d)
         end do    
    end do  

    do i = 1, 4, 1
        do j = 1, 4, 1
        if (j/=i)then
            lambda_0(i,i)=miu_0(i)-lambda_0(i,j)
        end if          
        end do
        lambda_0(i,i)=max(0.0,lambda_0(i,i))
    end do  
    print*,"---------8-----------"    
    print*,"likelihood is",likelihood
    print*,lambda_0
    FT=abs(likelihood-likelihoodlast)
    if ( likelihood-likelihoodmax>0.0 ) then
      lambda_em = lambda_0
      likelihoodmax=likelihood
      !exit
    end if      
    likelihoodlast=likelihood
    iterations=iterations+1
    print*,"iterations is",iterations
    if ( iterations>50 ) then
        exit
    end if
    end do
    globaltest=2.0*(likelihoodmax-likelihoodnull)
    return
end subroutine LRT_test

subroutine MM_Estimation(allsample_size,allsample,index_leave,miu_0,Cov_0,MMTest)
    use IMSL
    use funcsconstants
    implicit none
    integer allsample_size,sample_size,index_leave
    integer allsample(4,allsample_size),sample(4,allsample_size-1)
    real miu_0(1,4),Cov_0(4,4),Test_matrix(1,1),sample_leave(1,4),sample_tem(4,1)
    real MMTest
    integer i,j,k,d
    real Cov_test(4,4)
    real Inverse_Cov_exp(4,4)
    real Tem_matrix(4,4),identy(4,4)
    sample_size=allsample_size-1
    sample(:,:)=allsample(:,:)
    sample_tem=allsample(:,index_leave:index_leave)
	sample_leave=.t.sample_tem
    if ( index_leave<allsample_size ) then
        sample(:,index_leave:sample_size)=allsample(:,index_leave+1:allsample_size)      
    end if
    do i = 1, 4, 1
        miu_0(1,i)=1.0*sum(sample(i,1:sample_size))/sample_size
    end do
    Cov_0=0.0
    do i = 1, 4, 1
        do j=1,4,1
            do d = 1, sample_size, 1
                Cov_0(i,j)=Cov_0(i,j)+1.0*(sample(i,d)-miu_0(1,i))*(sample(j,d)-miu_0(1,j))
            end do
            Cov_0(i,j)=Cov_0(i,j)/sample_size
            Cov_0(j,i)=Cov_0(i,j)
        end do
    end do
    Cov_test=0.0
    do i = 1, 4, 1
        do j=1,4,1
            Cov_test(i,j)=1.0*(allsample(i,index_leave)-miu_0(1,i))*(allsample(j,index_leave)-miu_0(1,j))
            Cov_test(j,i)=Cov_test(i,j)
        end do
    end do
    Inverse_Cov_exp=.i.Cov_0
    Test_matrix=(sample_leave-miu_0) .x. (Inverse_Cov_exp) .xt.  (sample_leave-miu_0)
    identy=0.0
    do i = 1, 4, 1
        identy(i,i)=1.0
    end do    
    Tem_matrix=Cov_test .x. Inverse_Cov_exp
    Tem_matrix=Tem_matrix-identy
    Tem_matrix=Tem_matrix .x. Tem_matrix
    MMtest=0.0
    do i = 1, 4, 1
        MMtest=MMtest+Tem_matrix(i,i)
    end do
    MMtest=1.0*(0.5*MMtest+Test_matrix(1,1))    
    return
end subroutine MM_Estimation

subroutine Statistic_meanandcov(lambdaewma,n,sample,miu_ewma,Cov_ewma,miu_exp,Inverse_Cov_exp,Test)
        use IMSL
        use funcsconstants
        implicit none
        integer i,j,d,n
        integer sample(4,n)
        real Cov_ewma(4,4),Cov(4,4),lambdaewma
        real Inverse_Cov_exp(4,4)
        real Tem_matrix(4,4),identy(4,4)
        real Test
        real miu_ewma(1,4),miu_exp(1,4),Test_matrix(1,1)
        do i = 1, 4, 1
            miu_ewma(1,i)=lambdaewma*sum(sample(i,1:n))/n+(1.0-lambdaewma)*miu_ewma(1,i)
        end do
        !if ( sum(miu_ewma(1,:))<sum(miu_exp(1,:)) ) then
            !Test=0
            !return
        !end if
        Cov=0.0
        do i = 1, 4, 1
            do j=1,4,1
                do d = 1, n, 1
                    Cov(i,j)=Cov(i,j)+1.0*(sample(i,d)-miu_ewma(1,i))*(sample(j,d)-miu_ewma(1,j))
                end do
                Cov(i,j)=Cov(i,j)/n
                Cov_ewma(i,j)=lambdaewma*Cov(i,j)+(1.0-lambdaewma)*Cov_ewma(i,j)
                Cov_ewma(j,i)=Cov_ewma(i,j)
            end do
        end do        

        Test_matrix=0.0
        Test_matrix=(miu_ewma-miu_exp) .x. (Inverse_Cov_exp) .xt.  (miu_ewma-miu_exp)

        identy=0.0
        do i = 1, 4, 1
            identy(i,i)=1.0
        end do    
        Tem_matrix=Cov_ewma .x. Inverse_Cov_exp
        Tem_matrix=Tem_matrix-identy
        Tem_matrix=Tem_matrix .x. Tem_matrix
        Test=0.0
        do i = 1, 4, 1
            Test=Test+Tem_matrix(i,i)
        end do
        Test=1.0*n*(0.5*Test+Test_matrix(1,1))
        return
end subroutine Statistic_meanandcov
subroutine T_Estimation(allsample_size,allsample,index_leave,miu_0,Cov_0,MMTest)
    use IMSL
    use funcsconstants
    implicit none
    integer allsample_size,sample_size,index_leave
    integer allsample(4,allsample_size),sample(4,allsample_size-1)
    real miu_0(1,4),Cov_0(4,4),Test_matrix(1,1),sample_leave(1,4),sample_tem(4,1)
    real MMTest
    integer i,j,k,d
    real Inverse_Cov_exp(4,4)
    sample_size=allsample_size-1
    sample(:,:)=allsample(:,:)
    sample_tem=allsample(:,index_leave:index_leave)
    sample_leave=.t.sample_tem
    if ( index_leave<allsample_size ) then
        sample(:,index_leave:sample_size)=allsample(:,index_leave+1:allsample_size)      
    end if
    do i = 1, 4, 1
        miu_0(1,i)=1.0*sum(sample(i,1:sample_size))/sample_size
    end do
    Cov_0=0.0
    do i = 1, 4, 1
        do j=1,4,1
            do d = 1, sample_size, 1
                Cov_0(i,j)=Cov_0(i,j)+1.0*(sample(i,d)-miu_0(1,i))*(sample(j,d)-miu_0(1,j))
            end do
            Cov_0(i,j)=Cov_0(i,j)/sample_size
            Cov_0(j,i)=Cov_0(i,j)
        end do
    end do
    Inverse_Cov_exp=.i.Cov_0
    Test_matrix=(sample_leave-miu_0) .tx. (Inverse_Cov_exp) .x.  (sample_leave-miu_0)
    MMtest=Test_matrix(1,1)    
    return
end subroutine T_Estimation
subroutine T_square(lambdaewma,n,sample,miu_ewma,miu_exp,Inverse_Cov_exp,Test)
        use IMSL
        use funcsconstants
        implicit none
        integer i,j,d,n
        integer sample(4,n)
        real lambdaewma
        real Inverse_Cov_exp(4,4)
        real Test
        real miu_ewma(1,4),miu_exp(1,4),Test_matrix(1,1)
        do i = 1, 4, 1
            miu_ewma(1,i)=lambdaewma*sum(sample(i,1:n))/n+(1.0-lambdaewma)*miu_ewma(1,i)
        end do  
        !if ( sum(miu_ewma(1,:))<sum(miu_exp(1,:)) ) then
            !Test=0
            !return
        !end if
        Test_matrix=0.0
        Test_matrix=(miu_ewma-miu_exp) .x. (Inverse_Cov_exp) .xt.  (miu_ewma-miu_exp)
        Test=Test_matrix(1,1)
        return
end subroutine T_square

subroutine Chi_estimation(allsample_size,allsample,index_leave,midian,ICpr_cell,Chi_test)
    use IMSL
    use funcsconstants
    implicit none  
    integer i,j,k,count_tem,sample_size,allsample_size,index_leave
    integer allsample(4,allsample_size),sample(4,allsample_size-1),sample_order(4,allsample_size-1),cell(4),indef(12),level_dim(4),factor_dim(6),cells
    real midian(4),ICpr_cell(16),sample_cutting(16),Chi_test
    sample=allsample
	sample_size=allsample_size-1
    if ( index_leave<allsample_size ) then
        sample(:,index_leave:sample_size)=allsample(:,index_leave+1:allsample_size)      
    end if
    cells=16
    level_dim=2
    factor_dim=2
    indef=(/1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4/)
    count_tem=0
    do i = 1, 4, 1
        CALL SVIGN (sample_size, sample(i,:), sample_order(i,:))
        count_tem=sample_size/2
        midian(i)=1.0*(sample_order(i,count_tem)+sample_order(i,count_tem+1))/2.0
    end do    
    !print*,"---------1---------------"
    sample_cutting=0
    do i = 1, sample_size, 1
       count_tem=0
       cell=0
       do j = 1, 4, 1
           if ( sample(j,i)>midian(j) ) then
              cell(j)=1
           else
              cell(j)=0  
           end if 
           count_tem=count_tem+cell(j)*(2**(j-1))
        end do
        sample_cutting(count_tem+1)=sample_cutting(count_tem+1)+1     
    end do
        !print*,"---------2---------------"
    CALL PRPFT (4, level_dim, sample_cutting, 6, factor_dim, indef, 0.01,10, sample_cutting)
    ICpr_cell=sample_cutting/sample_size
    sample_cutting=0.0
    count_tem=0
    cell=0
    !    print*,"---------3---------------"

   do j = 1, 4, 1
       if ( allsample(j,index_leave)>=midian(j) ) then
          cell(j)=1
       else
          cell(j)=0  
       end if 
       count_tem=count_tem+cell(j)*(2**(j-1))
    end do 
    sample_cutting(count_tem+1)=1
    Chi_test=0.0
          !print*,"---------4---------------"
    do i = 1, cells, 1
       Chi_test=Chi_test+(((sample_cutting(i)-ICpr_cell(i))**2.0)/ICpr_cell(i))
    end do  
          !print*,"---------5---------------" 
    return
end subroutine Chi_estimation

subroutine Statistic_Chi(lambdaewma,n,sample,ICpr_cell,midian,sample_cutting_ewma,test)
    use IMSL
    use funcsconstants
    implicit none 
    integer i,j,count_tem,n
    integer cells
    integer sample(4,n),cell(4),sample_cutting(16)
    real midian(4),ICpr_cell(16),sample_cutting_ewma(16)
    real sample_exp(16),lambdaewma
    real Test
    cells=16
    !print*,sample
    sample_exp=ICpr_cell*n
    sample_cutting=0
    do i = 1, n, 1
       count_tem=0
       cell=0
       do j = 1, 4, 1
           if ( sample(j,i)>=midian(j) ) then
              cell(j)=1
           else
              cell(j)=0  
           end if 
           count_tem=count_tem+cell(j)*(2**(j-1))
        end do 
        sample_cutting(count_tem+1)=sample_cutting(count_tem+1)+1     
    end do
    !print*,"sample_cutting is",sample_cutting
    sample_cutting_ewma=lambdaewma*sample_cutting+(1.0-lambdaewma)*sample_cutting_ewma
    Test=0.0
    do i = 1, cells, 1
       Test=Test+(((sample_cutting_ewma(i)-(sample_exp(i)))**2.0)/sample_exp(i))
    end do
    return
end subroutine Statistic_Chi

subroutine D_estimation(allsample_size,allsample,index_leave,CL,D_test)
    use IMSL
    use funcsconstants
    implicit none  
    integer i,j,k,count_tem,sample_size,allsample_size,index_leave
    integer allsample(4,allsample_size),sample(4,allsample_size-1)
    real CL,D_test,mean_value(4)
    sample=allsample
    sample_size=allsample_size-1
    if ( index_leave<allsample_size ) then
        sample(:,index_leave:sample_size)=allsample(:,index_leave+1:allsample_size)      
    end if
    CL=0.0
    do i = 1, 4, 1
        mean_value(i)=1.0*sum(sample(i,:))/sample_size
        CL=CL+mean_value(i)
    end do 
    D_test=abs(sum(allsample(:,index_leave))-CL)
    return
end subroutine D_estimation

subroutine D_statistic(lambdaewma,n,sample,Test)
    use funcsconstants
    implicit none
    integer n
    integer sample(4,n)
    real mean_value(4)
    real Test,Test_tmp,lambdaewma
    integer i,j
    Test_tmp=0.0
    do i = 1, 4, 1
        mean_value(i)=1.0*sum(sample(i,:))/n
        Test_tmp=Test_tmp+mean_value(i)
    end do
  Test=lambdaewma*Test_tmp+(1.0-lambdaewma)*Test
  return
 
end subroutine D_statistic