module MPmodule
	!dec$objcomment lib:"Casefuncs.lib"
	implicit none
	interface 
subroutine Poisample_generation(lambda_0,n,sample)
    use IMSL
    implicit none
    integer,intent(in):: n
    real*8,intent(in):: lambda_0(4,4)
    integer,intent(out):: sample(4,n)
end subroutine Poisample_generation

subroutine MP_Estimation(allsample_size,allsample,index_leave,lambda_0,LRTvalue)
    use IMSL
    implicit none
    integer,intent(in)::allsample_size, index_leave
	integer,intent(in)::allsample(4,allsample_size)
    real*8,intent(out):: lambda_0(4,4),LRTvalue  
end subroutine MP_Estimation

subroutine LRT_test(lambdaewma,window_size,n,sample,lambda_IC,globaltest)
    use IMSL
    implicit none
    integer,intent(in)::window_size, n
	real*8,intent(in):: lambdaewma,lambda_IC(4,4)
	integer,intent(in)::sample(window_size,4,n)
    real*8,intent(out):: globaltest
end subroutine LRT_test   

subroutine MM_Estimation(allsample_size,allsample,index_leave,miu_0,Cov_0,MMTest)
    use IMSL
    implicit none
    integer,intent(in):: allsample_size,index_leave
    integer,intent(in):: allsample(4,allsample_size)
    real,intent(out):: miu_0(1,4),Cov_0(4,4),MMTest
end subroutine MM_Estimation

subroutine Statistic_meanandcov(lambdaewma,n,sample,miu_ewma,Cov_ewma,miu_exp,Inverse_Cov_exp,Test)
        use IMSL
        implicit none
        integer,intent(in)::n
        integer,intent(in):: sample(4,n)
        real,intent(in):: lambdaewma,Inverse_Cov_exp(4,4),miu_exp(1,4)
        real,intent(inout):: Cov_ewma(4,4),miu_ewma(1,4)
        real,intent(out):: Test  
end subroutine Statistic_meanandcov

subroutine T_Estimation(allsample_size,allsample,index_leave,miu_0,Cov_0,MMTest)
    use IMSL
    implicit none
    integer,intent(in):: allsample_size,index_leave
    integer,intent(in):: allsample(4,allsample_size)
    real,intent(out):: miu_0(1,4),Cov_0(4,4), MMTest 
end subroutine T_Estimation

subroutine T_square(lambdaewma,n,sample,miu_ewma,miu_exp,Inverse_Cov_exp,Test)
        use IMSL
        implicit none
        integer,intent(in):: n
        integer,intent(in):: sample(4,n)
        real,intent(in):: lambdaewma,Inverse_Cov_exp(4,4),miu_exp(1,4)
        real,intent(out):: Test
        real,intent(inout):: miu_ewma(1,4)
end subroutine T_square

subroutine Chi_estimation(allsample_size,allsample,index_leave,midian,ICpr_cell,Chi_test)
    use IMSL
    implicit none  
    integer,intent(in):: allsample_size,index_leave
    integer,intent(in):: allsample(4,allsample_size)
    real,intent(out):: midian(4),ICpr_cell(16),Chi_test
end subroutine Chi_estimation

subroutine Statistic_Chi(lambdaewma,n,sample,ICpr_cell,midian,sample_cutting_ewma,Test)
    use IMSL
    implicit none 
    integer,intent(in):: n
    integer,intent(in):: sample(4,n)
    real,intent(in):: midian(4),ICpr_cell(16)
    real,intent(inout)::sample_cutting_ewma(16)
    real,intent(in):: lambdaewma
    real,intent(out):: Test
end subroutine Statistic_Chi

subroutine D_estimation(allsample_size,allsample,index_leave,CL,D_test)
    use IMSL
    implicit none  
    integer,intent(in):: allsample_size,index_leave
    integer,intent(in)::allsample(4,allsample_size)
    real,intent(out):: CL,D_test

end subroutine D_estimation

subroutine D_statistic(lambdaewma,n,sample,Test)
    implicit none
    integer,intent(in):: n
    integer,intent(in):: sample(4,n)
    real,intent(in):: lambdaewma
    real,intent(out):: Test
end subroutine D_statistic
	end interface ! name

end module MPmodule