module Chartingmodule
	!dec$objcomment lib:"Charting_statistics.lib"
	implicit none
	interface 
subroutine LRT_test(dim,window_size,n,lambdaewma,sample,lambda_IC,globaltest)
    implicit none
    integer, intent(in):: dim,window_size,n             !dim=3 or 4
    integer, intent(in):: sample(window_size,dim,n)
    real, intent(in):: lambda_IC(dim,dim),lambdaewma
    real, intent(out):: globaltest
end subroutine LRT_test

subroutine MP_Estimation(sample_size,sample,lambda_est)
    implicit none
    integer, intent(in):: sample_size            
    integer, intent(in):: sample(3,sample_size)
    real, intent(out):: lambda_est(3,3)
end subroutine MP_Estimation


subroutine Statistic_meanandcov(dim,n,miu_test,miu_exp,Cov_test,Inverse_Cov_exp,Test)
        use IMSL
        implicit none
        integer, intent(in):: dim,n
        real, intent(in):: Cov_test(dim,dim)
        real, intent(in):: Inverse_Cov_exp(dim,dim)        
        real, intent(in):: miu_test(1,dim),miu_exp(1,dim)
        real,intent(out):: Test
end subroutine Statistic_meanandcov

subroutine MM_EstimationI(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    implicit none
    integer, intent(in):: dim,n
    integer, intent(in):: sample(dim,n)
    real, intent(in):: lambdaewma
    real, intent(inout):: miu_ewma(dim),Cov_ewma(dim,dim)
end subroutine MM_EstimationI

subroutine MM_EstimationII(dim,n,lambdaewma,sample,miu_ewma,Cov_ewma)
    implicit none
    integer, intent(in):: dim,n
    integer, intent(in):: sample(dim,n)
    real, intent(in):: lambdaewma
    real, intent(inout):: miu_ewma(dim),Cov_ewma(dim,dim)
end subroutine MM_EstimationII

subroutine MM_Estimation(dim,sample_size,sample,miu,Cov)
    implicit none
    integer, intent(in):: dim,sample_size
    integer, intent(in):: sample(dim,sample_size)
    real, intent(out):: miu(dim),Cov(dim,dim)
end subroutine MM_Estimation

subroutine Mean_EstimationI(dim,n,lambdaewma,sample,miu_ewma)
    implicit none
    integer, intent(in):: dim,n
    integer, intent(in):: sample(dim,n)
    real, intent(in):: lambdaewma
    real, intent(inout):: miu_ewma
end subroutine Mean_EstimationI

subroutine T_square(dim,n,miu_test,miu_exp,Inverse_Cov_exp,Test)
        use IMSL
        implicit none
        integer, intent(in):: dim,n
        real, intent(in):: Inverse_Cov_exp(dim,dim),miu_test(1,dim),miu_exp(1,dim)
        real, intent(out):: Test
end subroutine T_square

!*****************************************D-sum********************************************************

subroutine D_statistic(dim,n,lambdaewma,sample,Test)
    implicit none
    integer, intent(in):: dim,n
    integer, intent(in):: sample(dim,n)
    real, intent(in):: lambdaewma
    real, intent(out):: Test
end subroutine D_statistic

subroutine D_CLget(dim,miu_matrix,CL)
    implicit none
    integer, intent(in):: dim      
    real, intent(in):: miu_matrix(dim,dim)
    real ,intent(out):: CL
end subroutine D_CLget

!*********************************************Chi-square************************************************

subroutine Statistic_Pei(dim,n,lambdaewma,sample,ICpr_cell,midian,sample_cutting_ewma,test)
    use IMSL
    implicit none 
    integer, intent(in):: dim,n
    integer, intent(in):: sample(dim,n)
    real, intent(in):: midian(dim),ICpr_cell(2**dim),lambdaewma
    real, intent(inout):: sample_cutting_ewma(2**dim)
    real, intent(out):: Test
end subroutine Statistic_Pei

subroutine Median_get(dim,sample_size,sample,midian,ICpr_cell)
    use IMSL
    implicit none  
    integer, intent(in):: dim,sample_size
    integer, intent(in):: sample(dim,sample_size)
    real, intent(out):: midian(dim),ICpr_cell(2**dim)
end subroutine Median_get


    end interface ! name

end module Chartingmodule
