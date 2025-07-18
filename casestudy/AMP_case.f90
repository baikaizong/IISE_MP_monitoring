module constants
    implicit none
    integer, parameter::dim=4,all_size=360,IC_size=260,on_size=100,leave_size=260*0.9,window_size=50,n=1
    real*8,parameter:: lambdaewma=0.20
end module constants
program main
    use constants
	use MPmodule
    implicit none 
    integer allsample(dim,all_size),sample(dim,IC_size),index_leave,testsample(window_size,dim,n),onsample(dim,on_size)
    real*8 lambda_0(dim,dim),Testmax,lambda_IC(dim,dim),LRTvalue,globaltest
    real Limit,arl,std
    integer i,j,k
    open(unit=11, file="./data/countdata.txt")
    do i = 1, all_size, 1
      read(11,*), allsample(:,i)
     end do
    close(11)
    open(unit=10, file="./data/AMP_statistic.txt")
	open(unit=12, file="./data/AMP_lambda.txt")    
    sample=allsample(:,1:IC_size)
    onsample=allsample(:,IC_size+1:IC_size+on_size) 
    !i=IC_size  
    !index_leave=16    
	!write(12,*),"+-+-+-+-+-+-+-+--++-"
    !write(12,*),sample(:,index_leave)
    !sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    !i=i-1
    !index_leave=111     
	!write(12,*),"+-+-+-+-+-+-+-+--++-"
	!write(12,*),sample(:,index_leave)
    !sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    !i=i-1
    !index_leave=111    
	!write(12,*),"+-+-+-+-+-+-+-+--++-"
    !write(12,*),sample(:,index_leave)
    !sample(:,index_leave:i-1)=sample(:,index_leave+1:i)

    !i=i-1
    !index_leave=269    
	!write(12,*),"+-+-+-+-+-+-+-+--++-"
    !write(12,*),sample(:,index_leave)
    !sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    !i=i-1   
    !do while ( i>leave_size )
       ! Testmax=0.0
        !do j = 1, i, 1		
		   ! print*,j
           ! CALL MP_Estimation(i,sample(:,1:i),j,lambda_0,LRTvalue)
		!	print*,"**********************"	
			!print*,LRTvalue
           ! if ( Testmax<LRTvalue ) then
                !Testmax=LRTvalue
                !lambda_IC=lambda_0
                !index_leave=j
           !end if
        !end do
		!print*,"+-+-+-+-+-+-+-+--++-"
		!print*,sample(:,index_leave)      
        !write(12,*),"+-+-+-+-+-+-+-+--++-"
        !write(12,*),sample(:,index_leave)
		!print*,lambda_IC
        !if ( index_leave<i ) then
           ! sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
        !end if 
        !i=i-1
    !end do
	!write(12,*),lambda_IC
    lambda_IC=reshape((/5.67480598130887,0.101665431695977,7.916901720573267E-004,6.451025800736700E-002,0.101665431695977,2.76148549946504,2.349309644620400E-002,5.902732104778009E-002,7.916901720573267E-004,2.349309644620400E-002,3.07448095617222,0.126373744682486,6.451025800736700E-002,5.902732104778009E-002,0.126373744682486,4.06166044335170/),(/4,4/))   
    j=1
    do i = 1, window_size, 1
        call Poisample_generation(lambda_IC,n,testsample(i,:,:))
    end do
    CALL LRT_test(lambdaewma,window_size,n,testsample,lambda_IC,globaltest) 
    print*,globaltest
    do k = 1, on_size, 1        
        testsample(1:window_size-1,:,:)=testsample(2:window_size,:,:)
        testsample(window_size,:,:)=onsample(:,j:j+n-1)
        CALL LRT_test(lambdaewma,window_size,n,testsample,lambda_IC,globaltest) 
        print*,globaltest
        write(10,*),globaltest
        j=j+n
    end do
    stop
end program main

