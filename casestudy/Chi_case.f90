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
    real midian(4),pr_cell(16),sample_cutting_ewma(16),Chi_test,midian_IC(4),ICpr_cell(16),Chi_testmax
    integer i,j,k
    open(unit=11, file="./data/countdata.txt")
    open(unit=12, file="./data/Chi_IC.txt")
	print*,CHIIN(0.999,15.0)*0.2*0.8
	stop
    do i = 1, all_size, 1
      read(11,*), allsample(:,i)
     end do
     !print*,allsample
    close(11)
    open(unit=10, file="./data/Chi_statistic.txt")
  
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
        Chi_testmax=0.0
        do j = 1, i, 1		
		    print*,j
            CALL Chi_estimation(i,sample,j,midian,pr_cell,Chi_test)
			print*,"**********************"	
			print*,Chi_test
			print*,pr_cell
			print*,sum(pr_cell)
            if ( Chi_testmax<Chi_test ) then
                Chi_testmax=Chi_test
                midian_IC=midian
                ICpr_cell=pr_cell
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
	write(12,*),"midian_IC is", midian_IC
	write(12,*),"ICpr_cell is", ICpr_cell
    j=1
    sample_cutting_ewma=ICpr_cell*n
    do k = 1, on_size, 1 
        testsample(:,:)=onsample(:,j:j+n-1)
        call Statistic_Chi(lambdaewma,n,testsample,ICpr_cell,midian_IC,sample_cutting_ewma,globaltest)
        print*,globaltest
        write(10,*),globaltest     
        j=j+n
    end do
    stop
end program main

