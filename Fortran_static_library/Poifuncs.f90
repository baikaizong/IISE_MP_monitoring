module funcsconstants
    implicit none
    real,parameter:: epsilon=0.0000001
end module funcsconstants

!**************************AMP*****************************************************

subroutine LRT_testthreedim(window_size,n,lambdaewma,sample,lambda_IC,globaltest)
    use funcsconstants
    use IMSL
    implicit none
    integer window_size,n
    integer sample(window_size,3,n),max_val(3),max_dim(3)
    real lambda_IC(3,3),lambda_est(3,3),lambda_last(3,3),lambdaewma,lambdaco,miu(3)
    real globaltest
    integer i,j,k,d
    integer iterations
    real FT,likelihoodlast,likelihood,likelihoodnull,likelihoodmax
    real likelihood_array,S(3),S_ewma(3)
    real*8, allocatable :: pro_table(:,:,:)
    integer tmp(3,window_size)
    lambdaco=1.0-lambdaewma
    do i = 1, 3, 1
        do j = 1, window_size, 1
            tmp(i,j)=maxval(sample(j,i,:))
        end do
        max_dim(i)=maxval(tmp(i,:))
    end do

    allocate(pro_table(max_dim(1)+2,max_dim(2)+2,max_dim(3)+2))
    lambda_est=lambda_IC
    
    pro_table=0.000000000000001
    do i = 1, 3, 1
        max_val(i)=min(1.0*max_dim(i),34.0)
    end do    
    pro_table(2,2,2)=exp(-(lambda_est(1,1)+lambda_est(2,2)+lambda_est(3,3)+lambda_est(1,2)+lambda_est(1,3)+lambda_est(2,3))) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2)= pro_table(2,2,2)*((lambda_est(1,1)**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2)= pro_table(2,2,2)*((lambda_est(2,2)**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i)= pro_table(2,2,2)*((lambda_est(3,3)**(i-2))/FAC(i-2))
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2)=(lambda_est(1,1)*pro_table(i-1,j,2)+lambda_est(1,2)*pro_table(i-1,j-1,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j)=(lambda_est(1,1)*pro_table(i-1,2,j)+lambda_est(1,3)*pro_table(i-1,2,j-1))/(i-2)
        end do  
    end do 
   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j)=(lambda_est(2,2)*pro_table(2,i-1,j)+lambda_est(2,3)*pro_table(2,i-1,j-1))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k)=(lambda_est(1,1)*pro_table(i-1,j,k)+lambda_est(1,2)*pro_table(i-1,j-1,k)+lambda_est(1,3)*pro_table(i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 1, max_val(1)+2, 1
        do j = 1, max_val(2)+2, 1
            do k = 1, max_val(3)+2, 1
                if ( pro_table(i,j,k)<=0 ) then
                    pro_table(i,j,k)=0.000000000000001
                end if
            end do  
        end do    
    end do
    likelihood=0.0
    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lambdaewma*likelihood_array+lambdaco*likelihood 
    end if         
    end do
    likelihoodnull=likelihood 
    S=0.0
    do j = 1, window_size, 1
        do i = 1, 3, 1
           S(i)=1.0*sum(sample(j,i,:))/n         
           if ( j==1 ) then
            miu(i)=S(i)
           else
           miu(i)=lambdaewma*S(i)+lambdaco*miu(i) 
        end if   
        end do
    end do
    do j = 1, window_size, 1
        S=0.0
        do k = 1, n, 1
            S(1)=S(1)+(sample(j,1,k)-miu(1))*(sample(j,2,k)-miu(2))
            S(2)=S(2)+(sample(j,1,k)-miu(1))*(sample(j,3,k)-miu(3))
            S(3)=S(3)+(sample(j,2,k)-miu(2))*(sample(j,3,k)-miu(3))
        end do
        if ( j==1 ) then
            lambda_est(1,2)=S(1)
            lambda_est(1,3)=S(2)
            lambda_est(2,3)=S(3)
        else
            lambda_est(1,2)=lambdaewma*S(1)+lambdaco*lambda_est(1,2)
            lambda_est(1,3)=lambdaewma*S(2)+lambdaco*lambda_est(1,3)
            lambda_est(2,3)=lambdaewma*S(3)+lambdaco*lambda_est(2,3)
        end if
    end do
    lambda_est(1,2)=max(0.0,lambda_est(1,2)/n)
    lambda_est(1,3)=max(0.0,lambda_est(1,3)/n)
    lambda_est(2,3)=max(0.0,lambda_est(2,3)/n)
    lambda_est(2,1)=lambda_est(1,2)
    lambda_est(3,1)=lambda_est(1,3)
    lambda_est(3,2)=lambda_est(2,3)    
    lambda_est(1,1)=max(0.0,miu(1)-lambda_est(1,2)-lambda_est(1,3))
    lambda_est(2,2)=max(0.0,miu(2)-lambda_est(1,2)-lambda_est(2,3))
    lambda_est(3,3)=max(0.0,miu(3)-lambda_est(1,3)-lambda_est(2,3))
    FT=1.0
    likelihoodlast=likelihoodnull
    likelihoodmax=likelihoodlast
    iterations=0
    do while ( FT>epsilon)
    pro_table=0.000000000000001    
    pro_table(2,2,2)=exp(-(lambda_est(1,1)+lambda_est(2,2)+lambda_est(3,3)+lambda_est(1,2)+lambda_est(1,3)+lambda_est(2,3))) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2)= pro_table(2,2,2)*((lambda_est(1,1)**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2)= pro_table(2,2,2)*((lambda_est(2,2)**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i)= pro_table(2,2,2)*((lambda_est(3,3)**(i-2))/FAC(i-2))
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2)=(lambda_est(1,1)*pro_table(i-1,j,2)+lambda_est(1,2)*pro_table(i-1,j-1,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j)=(lambda_est(1,1)*pro_table(i-1,2,j)+lambda_est(1,3)*pro_table(i-1,2,j-1))/(i-2)
        end do  
    end do 
   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j)=(lambda_est(2,2)*pro_table(2,i-1,j)+lambda_est(2,3)*pro_table(2,i-1,j-1))/(i-2)
        end do  
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k)=(lambda_est(1,1)*pro_table(i-1,j,k)+lambda_est(1,2)*pro_table(i-1,j-1,k)+lambda_est(1,3)*pro_table(i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do

    likelihood=0.0
    do i = 1, max_val(1)+2, 1
        do j = 1, max_val(2)+2, 1
            do k = 1, max_val(3)+2, 1
                if ( pro_table(i,j,k)<=0 ) then
                    pro_table(i,j,k)=0.000000000000001
                end if
            end do  
        end do    
    end do

    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lambdaewma*likelihood_array+lambdaco*likelihood 
    end if         
    end do
    FT=likelihood-likelihoodlast
    if ( likelihoodmax-likelihood<0.0 ) then
      likelihoodmax=likelihood
    end if      
    likelihoodlast=likelihood    
    iterations=iterations+1
    if ( iterations>30 ) then
        exit
    end if
    S_ewma=0.0
    do j = 1, window_size, 1
    S=0.0
    do i = 1, n, 1
         S(1)=S(1)+(lambda_est(1,2)*pro_table(sample(j,1,i)+1,sample(j,2,i)+1,sample(j,3,i)+2)/pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
         S(2)=S(2)+(lambda_est(1,3)*pro_table(sample(j,1,i)+1,sample(j,2,i)+2,sample(j,3,i)+1)/pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
         S(3)=S(3)+(lambda_est(2,3)*pro_table(sample(j,1,i)+2,sample(j,2,i)+1,sample(j,3,i)+1)/pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2))
    end do
    if ( j==1 ) then
       S_ewma(1)=S(1)
       S_ewma(2)=S(2)
       S_ewma(3)=S(3)
    else
       S_ewma(1)=lambdaewma*S(1)+lambdaco*S_ewma(1)
       S_ewma(2)=lambdaewma*S(2)+lambdaco*S_ewma(2)
       S_ewma(3)=lambdaewma*S(3)+lambdaco*S_ewma(3)
    end if
    enddo
    lambda_last=lambda_est
    lambda_est(1,2)=max(0.0,S_ewma(1)/n)
    lambda_est(2,1)=lambda_est(1,2)
    lambda_est(1,3)=max(0.0,S_ewma(2)/n)
    lambda_est(3,1)=lambda_est(1,3)
    lambda_est(2,3)=max(0.0,S_ewma(3)/n)
    lambda_est(3,2)=lambda_est(2,3)
    lambda_est(1,1)=max(0.0,miu(1)-lambda_est(1,2)-lambda_est(1,3))
    lambda_est(2,2)=max(0.0,miu(2)-lambda_est(1,2)-lambda_est(2,3))
    lambda_est(3,3)=max(0.0,miu(3)-lambda_est(1,3)-lambda_est(2,3))
    end do
    globaltest=2.0*(likelihoodmax-likelihoodnull)
    return
end subroutine LRT_testthreedim

subroutine LRT_testfourdim(window_size,n,lambdaewma,sample,lambda_IC,globaltest)
    use funcsconstants
    use IMSL
    implicit none
    integer window_size,n
    integer sample(window_size,4,n),max_val(4),max_dim(4),tem_index(4)
    real lambda_IC(4,4),lambda_est(4,4),lambda_last(4,4),tem_frac,miu(4),miu_tem(4),lambdaewma,lambdaco,lambda_em(4,4)
    real globaltest
    integer i,j,k,d
    integer iterations
    real FT,likelihoodlast,likelihood,likelihoodnull,likelihoodmax
    real likelihood_array,S(4,4),S_ewma(4,4),S_sum(4)
    real, allocatable :: pro_table(:,:,:,:)
    integer tmp(4,window_size)
    lambdaco=1.0-lambdaewma
    lambda_est=lambda_IC
    do i = 1, 4, 1
        do j = 1, window_size, 1
            tmp(i,j)=maxval(sample(j,i,:))
        end do
        max_dim(i)=maxval(tmp(i,:))
    end do

    allocate(pro_table(max_dim(1)+2,max_dim(2)+2,max_dim(3)+2,max_dim(4)+2))
    pro_table=0.000000000000001
    do i = 1, 4, 1
        max_val(i)=min(1.0*max_dim(i),34.0)
    end do     
  
    pro_table(2,2,2,2)=exp(-(lambda_est(1,1)+lambda_est(2,2)+lambda_est(3,3)+lambda_est(4,4)+lambda_est(1,2)+lambda_est(1,3)+lambda_est(1,4)+lambda_est(2,3)+lambda_est(2,4)+lambda_est(3,4))) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2,2)= pro_table(2,2,2,2)*(lambda_est(1,1)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2,2)= pro_table(2,2,2,2)*(lambda_est(2,2)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i,2)= pro_table(2,2,2,2)*(lambda_est(3,3)**(i-2))/FAC(i-2)
    end do 
    do i = 3, max_val(4)+2, 1
        pro_table(2,2,2,i)= pro_table(2,2,2,2)*(lambda_est(4,4)**(i-2))/FAC(i-2)
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2,2)=(lambda_est(1,1)*pro_table(i-1,j,2,2)+lambda_est(1,2)*pro_table(i-1,j-1,2,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j,2)=(lambda_est(1,1)*pro_table(i-1,2,j,2)+lambda_est(1,3)*pro_table(i-1,2,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(i,2,2,j)=(lambda_est(1,1)*pro_table(i-1,2,2,j)+lambda_est(1,4)*pro_table(i-1,2,2,j-1))/(i-2)
        end do  
    end do   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j,2)=(lambda_est(2,2)*pro_table(2,i-1,j,2)+lambda_est(2,3)*pro_table(2,i-1,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,i,2,j)=(lambda_est(2,2)*pro_table(2,i-1,2,j)+lambda_est(2,4)*pro_table(2,i-1,2,j-1))/(i-2)
        end do  
    end do 
    do i = 3, max_val(3)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,2,i,j)=(lambda_est(3,3)*pro_table(2,2,i-1,j)+lambda_est(3,4)*pro_table(2,2,i-1,j-1))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k,2)=((lambda_est(1,1))*pro_table(i-1,j,k,2)+lambda_est(1,2)*pro_table(i-1,j-1,k,2)+lambda_est(1,3)*pro_table(i-1,j,k-1,2))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,j,2,k)=((lambda_est(1,1))*pro_table(i-1,j,2,k)+lambda_est(1,2)*pro_table(i-1,j-1,2,k)+lambda_est(1,4)*pro_table(i-1,j,2,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,2,j,k)=((lambda_est(1,1))*pro_table(i-1,2,j,k)+lambda_est(1,3)*pro_table(i-1,2,j-1,k)+lambda_est(1,4)*pro_table(i-1,2,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(2,i,j,k)=((lambda_est(2,2))*pro_table(2,i-1,j,k)+lambda_est(2,3)*pro_table(2,i-1,j-1,k)+lambda_est(2,4)*pro_table(2,i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
                do d = 3, max_val(4)+2, 1
               pro_table(i,j,k,d)=((lambda_est(1,1))*pro_table(i-1,j,k,d)+lambda_est(1,2)*pro_table(i-1,j-1,k,d)+lambda_est(1,3)*pro_table(i-1,j,k-1,d)+lambda_est(1,4)*pro_table(i-1,j,k,d-1))/(i-2)
               end do
            end do  
        end do    
    end do
    likelihood=0.0
    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2,sample(j,4,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lambdaewma*likelihood_array+lambdaco*likelihood 
    end if         
    end do
    likelihoodnull=likelihood 
    S_sum=0.0
    do j = 1, window_size, 1
        do i = 1, 4, 1
           S_sum(i)=1.0*sum(sample(j,i,:))/n         
           if ( j==1 ) then
            miu(i)=S_sum(i)
           else
           miu(i)=lambdaewma*S_sum(i)+lambdaco*miu(i) 
        end if   
        end do
    end do
    miu_tem=miu
     do i = 1, 4, 1
        do j = 1, 4, 1
        if (j/=i)then
        miu_tem(i)=miu_tem(i)-lambda_est(i,j)
        end if          
        end do
        lambda_est(i,i)=max(0.0,miu_tem(i))
    end do
    FT=1.0
    likelihoodlast=likelihoodnull
    likelihoodmax=likelihoodlast
    iterations=0
    do while ( FT>epsilon)
    pro_table=0.000000000000001

    pro_table(2,2,2,2)=exp(-(lambda_est(1,1)+lambda_est(2,2)+lambda_est(3,3)+lambda_est(4,4)+lambda_est(1,2)+lambda_est(1,3)+lambda_est(1,4)+lambda_est(2,3)+lambda_est(2,4)+lambda_est(3,4))) 

    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2,2)= pro_table(2,2,2,2)*(lambda_est(1,1)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2,2)= pro_table(2,2,2,2)*(lambda_est(2,2)**(i-2))/FAC(i-2)
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i,2)= pro_table(2,2,2,2)*(lambda_est(3,3)**(i-2))/FAC(i-2)
    end do 
    do i = 3, max_val(4)+2, 1
        pro_table(2,2,2,i)= pro_table(2,2,2,2)*(lambda_est(4,4)**(i-2))/FAC(i-2)
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2,2)=(lambda_est(1,1)*pro_table(i-1,j,2,2)+lambda_est(1,2)*pro_table(i-1,j-1,2,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j,2)=(lambda_est(1,1)*pro_table(i-1,2,j,2)+lambda_est(1,3)*pro_table(i-1,2,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(i,2,2,j)=(lambda_est(1,1)*pro_table(i-1,2,2,j)+lambda_est(1,4)*pro_table(i-1,2,2,j-1))/(i-2)
        end do  
    end do   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j,2)=(lambda_est(2,2)*pro_table(2,i-1,j,2)+lambda_est(2,3)*pro_table(2,i-1,j-1,2))/(i-2)
        end do  
    end do 
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,i,2,j)=(lambda_est(2,2)*pro_table(2,i-1,2,j)+lambda_est(2,4)*pro_table(2,i-1,2,j-1))/(i-2)
        end do  
    end do 
    do i = 3, max_val(3)+2, 1
        do j = 3, max_val(4)+2, 1
            pro_table(2,2,i,j)=(lambda_est(3,3)*pro_table(2,2,i-1,j)+lambda_est(3,4)*pro_table(2,2,i-1,j-1))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k,2)=((lambda_est(1,1))*pro_table(i-1,j,k,2)+lambda_est(1,2)*pro_table(i-1,j-1,k,2)+lambda_est(1,3)*pro_table(i-1,j,k-1,2))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,j,2,k)=((lambda_est(1,1))*pro_table(i-1,j,2,k)+lambda_est(1,2)*pro_table(i-1,j-1,2,k)+lambda_est(1,4)*pro_table(i-1,j,2,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(i,2,j,k)=((lambda_est(1,1))*pro_table(i-1,2,j,k)+lambda_est(1,3)*pro_table(i-1,2,j-1,k)+lambda_est(1,4)*pro_table(i-1,2,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            do k = 3, max_val(4)+2, 1
               pro_table(2,i,j,k)=((lambda_est(2,2))*pro_table(2,i-1,j,k)+lambda_est(2,3)*pro_table(2,i-1,j-1,k)+lambda_est(2,4)*pro_table(2,i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
                do d = 3, max_val(4)+2, 1
               pro_table(i,j,k,d)=((lambda_est(1,1))*pro_table(i-1,j,k,d)+lambda_est(1,2)*pro_table(i-1,j-1,k,d)+lambda_est(1,3)*pro_table(i-1,j,k-1,d)+lambda_est(1,4)*pro_table(i-1,j,k,d-1))/(i-2)
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
    likelihood=0.0
    do j = 1, window_size, 1
    likelihood_array=0.0
    do i = 1, n, 1
       likelihood_array=likelihood_array+log(pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2,sample(j,4,i)+2))
    end do
    if ( j==1 ) then
       likelihood=likelihood_array
    else
       likelihood=lambdaewma*likelihood_array+lambdaco*likelihood 
    end if         
    end do
    FT=likelihood-likelihoodlast
    if ( likelihoodmax-likelihood<0.0 ) then
      lambda_em = lambda_est
      likelihoodmax=likelihood
    end if      
    likelihoodlast=likelihood
    iterations=iterations+1
    if ( iterations>30 ) then
        exit
    end if    
    S_ewma=0.0
    do j = 1, window_size, 1
    S=0.0
    do i = 1, n, 1     
        tem_frac=pro_table(sample(j,1,i)+2,sample(j,2,i)+2,sample(j,3,i)+2,sample(j,4,i)+2)
        do k = 1, 4-1, 1
            do d = k+1, 4, 1
                tem_index=sample(j,:,i)+2
                tem_index(k)=tem_index(k)-1
                tem_index(d)=tem_index(d)-1
                S(k,d)=S(k,d)+(lambda_est(k,d))*pro_table(tem_index(1),tem_index(2),tem_index(3),tem_index(4))/tem_frac
            end do    
        end do
    end do
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
    lambda_last=lambda_est
   do k = 1, 4-1, 1
        do d = k+1, 4, 1
           lambda_est(k,d)=max(0.0,S_ewma(k,d)/n)
           lambda_est(d,k)=lambda_est(k,d)
         end do    
    end do  
    miu_tem=miu
     do i = 1, 4, 1
        do j = 1, 4, 1
        if (j/=i)then
        miu_tem(i)=miu_tem(i)-lambda_est(i,j)
        end if          
        end do
        lambda_est(i,i)=max(0.0,miu_tem(i))
    end do  
    end do
    globaltest=2.0*(likelihoodmax-likelihoodnull)
    return
end subroutine LRT_testfourdim

subroutine LRT_test(dim,window_size,n,lambdaewma,sample,lambda_IC,globaltest)
    use funcsconstants
    implicit none
    integer dim,window_size,n
    integer sample(window_size,dim,n)
    real lambda_IC(dim,dim),lambdaewma
    real globaltest
    select case (dim)
        case (3)
          call LRT_testthreedim(window_size,n,lambdaewma,sample,lambda_IC,globaltest)
        case (4)  
        call LRT_testfourdim(window_size,n,lambdaewma,sample,lambda_IC,globaltest)
    end select
    return
end subroutine LRT_test

subroutine MP_Estimation(sample_size,sample,lambda_est)
    use IMSL
    use funcsconstants
    implicit none
    integer sample_size
    integer sample(3,sample_size),max_val(3)
    real lambda_est(3,3),lambda_last(3,3),miu(3)
    integer i,j,k
    integer iterations
    real*8 FT,likelihoodlast,likelihood
    real*8 S(3)
    real*8, allocatable :: pro_table(:,:,:)
  
    do i = 1, 3, 1
        max_val(i)=maxval(sample(i,:))
    end do

    allocate(pro_table(max_val(1)+2,max_val(2)+2,max_val(3)+2))

    do i = 1, 3, 1
       miu(i)=1.0*sum(sample(i,:))/sample_size          
    end do

    S=0.0
    do j = 1, sample_size, 1
        S(1)=S(1)+(sample(1,j)-miu(1))*(sample(2,j)-miu(2))
        S(2)=S(2)+(sample(1,j)-miu(1))*(sample(3,j)-miu(3))
        S(3)=S(3)+(sample(2,j)-miu(2))*(sample(3,j)-miu(3))
    end do
    lambda_est(1,2)=max(0.0,S(1)/sample_size)
    lambda_est(1,3)=max(0.0,S(2)/sample_size)
    lambda_est(2,3)=max(0.0,S(3)/sample_size)
    lambda_est(2,1)=lambda_est(1,2)
    lambda_est(3,1)=lambda_est(1,3)
    lambda_est(3,2)=lambda_est(2,3) 
    lambda_est(1,1)=max(0.0,miu(1)-lambda_est(1,2)-lambda_est(1,3))
    lambda_est(2,2)=max(0.0,miu(2)-lambda_est(1,2)-lambda_est(2,3))
    lambda_est(3,3)=max(0.0,miu(3)-lambda_est(1,3)-lambda_est(2,3))
    FT=1.0
    likelihoodlast=-10000000000000.0
    iterations=0
    do while ( FT>0.00000000001)
    pro_table=0.0

    pro_table(2,2,2)=exp(-(lambda_est(1,1)+lambda_est(2,2)+lambda_est(3,3)+lambda_est(1,2)+lambda_est(1,3)+lambda_est(2,3))) 
    do i = 3, max_val(1)+2, 1
        pro_table(i,2,2)= pro_table(2,2,2)*(1.0*((lambda_est(1,1))**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(2)+2, 1
        pro_table(2,i,2)= pro_table(2,2,2)*(1.0*((lambda_est(2,2))**(i-2))/FAC(i-2))
    end do
    do i = 3, max_val(3)+2, 1
        pro_table(2,2,i)= pro_table(2,2,2)*(1.0*((lambda_est(3,3))**(i-2))/FAC(i-2))
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            pro_table(i,j,2)=1.0*((lambda_est(1,1))*pro_table(i-1,j,2)+lambda_est(1,2)*pro_table(i-1,j-1,2))/(i-2)
        end do  
    end do 

    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(i,2,j)=1.0*((lambda_est(1,1))*pro_table(i-1,2,j)+lambda_est(1,3)*pro_table(i-1,2,j-1))/(i-2)
        end do  
    end do 
   
    do i = 3, max_val(2)+2, 1
        do j = 3, max_val(3)+2, 1
            pro_table(2,i,j)=1.0*((lambda_est(2,2))*pro_table(2,i-1,j)+lambda_est(2,3)*pro_table(2,i-1,j-1))/(i-2)
        end do  
    end do 
    do i = 3, max_val(1)+2, 1
        do j = 3, max_val(2)+2, 1
            do k = 3, max_val(3)+2, 1
               pro_table(i,j,k)=((lambda_est(1,1))*pro_table(i-1,j,k)+lambda_est(1,2)*pro_table(i-1,j-1,k)+lambda_est(1,3)*pro_table(i-1,j,k-1))/(i-2)
            end do  
        end do    
    end do
    likelihood=0.0
    do i = 1, sample_size, 1
       likelihood=likelihood+log(pro_table(sample(1,i)+2,sample(2,i)+2,sample(3,i)+2))
    end do
      S=0.0
    do i = 1, sample_size, 1
         S(1)=S(1)+(lambda_est(1,2)*pro_table(sample(1,i)+1,sample(2,i)+1,sample(3,i)+2)/pro_table(sample(1,i)+2,sample(2,i)+2,sample(3,i)+2))
         S(2)=S(2)+(lambda_est(1,3)*pro_table(sample(1,i)+1,sample(2,i)+2,sample(3,i)+1)/pro_table(sample(1,i)+2,sample(2,i)+2,sample(3,i)+2))
         S(3)=S(3)+(lambda_est(2,3)*pro_table(sample(1,i)+2,sample(2,i)+1,sample(3,i)+1)/pro_table(sample(1,i)+2,sample(2,i)+2,sample(3,i)+2))
    end do
    lambda_last=lambda_est
    !print*,"lambda_est is",lambda_est
    lambda_est(1,2)=max(0.0,S(1)/sample_size)
    lambda_est(2,1)=lambda_est(1,2)
    lambda_est(1,3)=max(0.0,S(2)/sample_size)
    lambda_est(3,1)=lambda_est(1,3)
    lambda_est(2,3)=max(0.0,S(3)/sample_size)
    lambda_est(3,2)=lambda_est(2,3)
    lambda_est(1,1)=max(0.0,miu(1)-lambda_est(1,2)-lambda_est(1,3))
    lambda_est(2,2)=max(0.0,miu(2)-lambda_est(1,2)-lambda_est(2,3))
    lambda_est(3,3)=max(0.0,miu(3)-lambda_est(1,3)-lambda_est(2,3))            
    FT=abs(likelihood-likelihoodlast)
    if ( likelihood-likelihoodlast<0.0 ) then
      lambda_est = lambda_last
      exit
    end if      
    likelihoodlast=likelihood
    iterations=iterations+1
    if ( iterations>100 ) then
        exit
    end if
    end do
    return        
end subroutine MP_Estimation

!**************************SMP*****************************************************

subroutine Statistic_meanandcov(dim,n,miu_test,miu_exp,Cov_test,Inverse_Cov_exp,Test)
        use IMSL
        use funcsconstants
        implicit none
        integer i,dim,n
        real Cov_test(dim,dim)
        real Inverse_Cov_exp(dim,dim)
        real Tem_matrix(dim,dim),identy(dim,dim)
        real Test,Trace
        real miu_test(1,dim),miu_exp(1,dim),Test_matrix(1,1)


        Test_matrix=0.0
        Test_matrix=(miu_test-miu_exp) .x. (Inverse_Cov_exp) .xt.  (miu_test-miu_exp)

        identy=0.0
        do i = 1, dim, 1
            identy(i,i)=1.0
        end do    
        Tem_matrix=Cov_test .x. Inverse_Cov_exp
        Tem_matrix=Tem_matrix-identy
        Tem_matrix=Tem_matrix .x. Tem_matrix
        Test=Trace(dim,Tem_matrix)
        Test=1.0*n*(0.5*Test+Test_matrix(1,1))
        return
end subroutine Statistic_meanandcov

function Trace(dim,matrixcu)
        use funcsconstants
        implicit none        
        integer i,dim
        real matrixcu(dim,dim),trace
        Trace=0.0
        do i = 1, dim, 1
            Trace=Trace+matrixcu(i,i)
        end do
        return
end function Trace

subroutine MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    use funcsconstants
    implicit none
    integer dim,n
    integer i,j,d
    integer sample(dim,n)
    real miu_ewma(dim),Cov_ewma(dim,dim),Cov(dim,dim),lambdaewma,lambdaco
   lambdaco=1.0-lambdaewma
    do i = 1, dim, 1
        miu_ewma(i)=lambdaco*miu_ewma(i)+lambdaewma*sum(sample(i,:))/n
    end do
    Cov=0.0
    do i = 1, dim, 1
    do j=i,dim,1
        do d = 1, n, 1
            Cov(i,j)=Cov(i,j)+1.0*(sample(i,d)-miu_ewma(i))*(sample(j,d)-miu_ewma(j))
        end do
        Cov(i,j)=Cov(i,j)/n
        Cov(j,i)=Cov(i,j)
    end do
    end do
    Cov_ewma=lambdaewma*Cov+lambdaco*Cov_ewma 
    do i = 1, dim, 1
        Cov_ewma(i,i)=miu_ewma(i)
    end do
    return
end subroutine MM_EstimationI

subroutine MM_EstimationII(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    use funcsconstants
    implicit none
    integer dim,n
    integer i,j,d
    integer sample(dim,n)
    real miu_ewma(dim),Cov_ewma(dim,dim),Cov(dim,dim),lambdaewma,lambdaco
   lambdaco=1.0-lambdaewma
    do i = 1, dim, 1
        miu_ewma(i)=lambdaco*miu_ewma(i)+lambdaewma*sum(sample(i,:))/n
    end do
    Cov=0.0
    do i = 1, dim, 1
    do j=i,dim,1
        do d = 1, n, 1
            Cov(i,j)=Cov(i,j)+1.0*(sample(i,d)-miu_ewma(i))*(sample(j,d)-miu_ewma(j))
        end do
        Cov(i,j)=Cov(i,j)/n
        Cov(j,i)=Cov(i,j)
    end do
    end do
    Cov_ewma=lambdaewma*Cov+lambdaco*Cov_ewma 
    return
end subroutine MM_EstimationII

subroutine MM_Estimation(dim,sample_size,sample,miu,Cov)
    use funcsconstants
    implicit none
    integer dim,sample_size
    integer i,j,d
    integer sample(dim,sample_size)
    real miu(dim),Cov(dim,dim)
    do i = 1, dim, 1
        miu(i)=1.0*sum(sample(i,:))/sample_size
    end do
    Cov=0.0
    do i = 1, dim, 1
    do j=i,dim,1
        do d = 1, sample_size, 1
            Cov(i,j)=Cov(i,j)+1.0*(sample(i,d)-miu(i))*(sample(j,d)-miu(j))
        end do
        Cov(i,j)=Cov(i,j)/sample_size
        Cov(j,i)=Cov(i,j)
    end do
    end do
    do i = 1, dim, 1
        Cov(i,i)=miu(i)
    end do
    return
end subroutine MM_Estimation

!*******************************T-square**************************************

subroutine Mean_EstimationI(dim,n,lambdaewma,sample,miu_ewma)
    use funcsconstants
    implicit none
    integer i,j,dim,n
    integer sample(dim,n)
    real miu_ewma(dim),lambdaewma
    do i = 1, dim, 1
        miu_ewma(i)=(1.0-lambdaewma)*miu_ewma(i)+lambdaewma*sum(sample(i,:))/n
    end do
    return
end subroutine Mean_EstimationI

subroutine T_square(dim,n,miu_test,miu_exp,Inverse_Cov_exp,Test)
        use IMSL
        use funcsconstants
        implicit none
        integer dim,i,n
        real Inverse_Cov_exp(dim,dim)
        real Tem_matrix(dim,dim),identy(dim,dim)
        real Test
        real miu_test(1,dim),miu_exp(1,dim),Test_matrix(1,1)
        Test_matrix=0.0
        Test_matrix=(miu_test-miu_exp) .x. (Inverse_Cov_exp) .xt.  (miu_test-miu_exp)
        Test=1.0*n*(Test_matrix(1,1))
        return
end subroutine T_square

!*****************************************D-sum********************************************************

subroutine D_statistic(dim,n,lambdaewma,sample,Test)
    use funcsconstants
    implicit none
    integer dim,n
    integer sample(dim,n)
    real mean(dim)
    real Test,Test_tmp,lambdaewma
    integer i,j
    Test_tmp=0.0
    do i = 1, dim, 1
        mean(i)=1.0*sum(sample(i,:))/n
        Test_tmp=Test_tmp+mean(i)
    end do
  Test=lambdaewma*Test_tmp+(1.0-lambdaewma)*Test
  return
 
end subroutine D_statistic

subroutine D_CLget(dim,miu_matrix,CL)
    use funcsconstants
    implicit none
    integer i,j,dim      
    real miu_matrix(dim,dim),CL
    CL=0
    do i = 1, dim, 1
      do j = 1, dim, 1
        CL=CL+miu_matrix(i,j)
      end do 
    end do 
    return   
end subroutine D_CLget

!*********************************************Chi-square************************************************

subroutine Statistic_Pei(dim,n,lambdaewma,sample,ICpr_cell,midian,sample_cutting_ewma,test)
    use IMSL
    use funcsconstants
    implicit none 
    integer i,j,count_tem,dim,n
    integer cells
    integer sample(dim,n),cell(dim),sample_cutting(2**dim)
    real midian(dim),ICpr_cell(2**dim),sample_cutting_ewma(2**dim)
    real sample_obs(1,2**dim),sample_exp(2**dim),diag_matrix(2**dim,2**dim),lambdaewma
    real Test
    cells=2**dim
    !print*,sample
    sample_exp=ICpr_cell*n
    sample_cutting=0
    do i = 1, n, 1
       count_tem=0
       cell=0
       do j = 1, dim, 1
           if ( sample(j,i)>=midian(j) ) then
              cell(j)=1
           else
              cell(j)=0  
           end if 
           count_tem=count_tem+cell(j)*(2**(j-1))
        end do 
        sample_cutting(count_tem+1)=sample_cutting(count_tem+1)+1     
    end do
    sample_cutting_ewma=lambdaewma*sample_cutting+(1.0-lambdaewma)*sample_cutting_ewma
    Test=0.0
    do i = 1, cells, 1
       Test=Test+(((sample_cutting_ewma(i)-(sample_exp(i)))**2.0)/sample_exp(i))
    end do
    return
end subroutine Statistic_Pei

subroutine Median_get(dim,sample_size,sample,midian,ICpr_cell)
    use IMSL
    use funcsconstants
    implicit none  
    integer i,j,k,count_tem,dim,sample_size
    integer sample(dim,sample_size),sample_order(dim,sample_size),cell(dim),indef(2*dim),level_dim(dim)
    real midian(dim),ICpr_cell(2**dim),sample_cutting(2**dim)
    level_dim=2
    indef=(/1, 2, 1, 3, 2, 3/)
    count_tem=0
    do i = 1, dim, 1
        CALL SVIGN (sample_size, sample(i,:), sample_order(i,:))
        count_tem=sample_size/2
        midian(i)=1.0*(sample_order(i,count_tem)+sample_order(i,count_tem+1))/2.0
    end do    
    !print*,midian
    sample_cutting=0
    do i = 1, sample_size, 1
       count_tem=0
       cell=0
       do j = 1, dim, 1
           if ( sample(j,i)>=midian(j) ) then
              cell(j)=1
           else
              cell(j)=0  
           end if 
           count_tem=count_tem+cell(j)*(2**(j-1))
        end do
        sample_cutting(count_tem+1)=sample_cutting(count_tem+1)+1     
    end do
    CALL PRPFT (dim, level_dim, sample_cutting, 3, level_dim, indef, 0.01,5, sample_cutting)
    ICpr_cell=sample_cutting/sample_size
    return
end subroutine Median_get




