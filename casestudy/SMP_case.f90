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
    real miu_0(1,4),Cov_0(dim,dim),miu_IC(1,4),Cov_IC(dim,dim),Testmax,MMTest,globaltest,Inverse_Cov_exp(dim,dim)
    real Limit,arl,std
    integer i,j,k
    open(unit=11, file="./data/countdata.txt")
    open(unit=12, file="./data/SMP_IC.txt")
    do i = 1, all_size, 1
      read(11,*), allsample(:,i)
     end do
     !print*,allsample
    close(11)
    sample=allsample(:,1:IC_size)
    onsample=allsample(:,IC_size+1:IC_size+on_size)     
    open(unit=10, file="./data/SMP_statistic.txt")
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
        Testmax=0.0
        do j = 1, i, 1		
		    print*,j
            CALL MM_Estimation(i,sample,j,miu_0,Cov_0,MMTest)
			print*,"**********************"	
			print*,MMTest
            if ( Testmax<MMTest ) then
                Testmax=MMTest
                miu_IC=miu_0
                Cov_IC=Cov_0
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
	write(12,*),"miu_IC is", miu_IC
	write(12,*),"Cov_IC is", Cov_IC
    Inverse_Cov_exp=.i.Cov_IC
    j=1 
    do k = 1, on_size, 1 
        testsample(:,:)=onsample(:,j:j+n-1)
        call Statistic_meanandcov(lambdaewma,n,testsample,miu_0,Cov_0,miu_IC,Inverse_Cov_exp,globaltest) 
        print*,globaltest
        write(10,*),globaltest   
        j=j+n
    end do
    stop
end program main

