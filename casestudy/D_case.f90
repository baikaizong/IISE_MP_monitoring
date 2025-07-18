module constants
    implicit none
    integer, parameter::dim=4,all_size=360,IC_size=260,on_size=100,leave_size=260*0.9,n=1
    real,parameter:: lambdaewma=0.20
end module constants
program main
    use constants
	use IMSL
	use MPmodule
    implicit none 
    integer allsample(dim,all_size),sample(dim,IC_size),index_leave,onsample(dim,on_size),testsample(dim,n)
    real Limit,arl,std,globaltest
    real D_test,D_testmax,CL,CL_IC
    integer i,j,k
    open(unit=11, file="./data/countdata.txt")
    open(unit=12, file="./data/D_IC.txt")

    do i = 1, all_size, 1
      read(11,*), allsample(:,i)
     end do
     !print*,allsample
    close(11)
    open(unit=10, file="./data/D_statistic.txt")
  
    sample=allsample(:,1:IC_size)
    onsample=allsample(:,IC_size+1:IC_size+on_size)   
    i=IC_size
    index_leave=16    
	write(12,*),"+-+-+-+-+-+-+-+--++-"
    write(12,*),sample(:,index_leave)
    sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    i=i-1
    index_leave=111     
	write(12,*),"+-+-+-+-+-+-+-+--++-"
	write(12,*),sample(:,index_leave)
    sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    i=i-1
    index_leave=111    
	write(12,*),"+-+-+-+-+-+-+-+--++-"
    write(12,*),sample(:,index_leave)
    sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    i=i-1
    index_leave=269    
	write(12,*),"+-+-+-+-+-+-+-+--++-"
    write(12,*),sample(:,index_leave)
    sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
    i=i-1 
    do while ( i>leave_size )
        D_testmax=0.0
        do j = 1, i, 1		
		    print*,j
            CALL D_estimation(i,sample,j,CL,D_test)
			print*,"**********************"	
			print*,D_test
            if ( D_testmax<D_test ) then
                D_testmax=D_test
                CL_IC=CL
                index_leave=j
            end if
        end do
		print*,sample(:,index_leave)
	    write(12,*),"+-+-+-+-+-+-+-+--++-"
        write(12,*),sample(:,index_leave)
        if ( index_leave<i ) then
            sample(:,index_leave:i-1)=sample(:,index_leave+1:i)
        end if 
        i=i-1
    end do
	write(12,*),"***************************"
	write(12,*),"CL_IC is", CL_IC
    j=1
    globaltest=CL
    do k = 1, on_size, 1 
        testsample(:,:)=onsample(:,j:j+n-1)
        call D_statistic(lambdaewma,n,testsample,globaltest)
        print*,globaltest
        write(10,*),globaltest     
        j=j+n
    end do
    stop
end program main

