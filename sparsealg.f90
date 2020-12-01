module sparsealg
  use precision
  use matdef
  implicit none
  private

  public :: matvec

contains

  subroutine matvec(A,x,y)
    type(rCSR), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)

    integer :: k, i

    y = 0.0_dp
    
    do k = 1, A%nrow
       do i = A%rowpnt(k), A%rowpnt(k+1)-1
          y(k) = y(k) + A%nzval(i)*x(A%colind(i))
       end do
    end do     
  end subroutine matvec

end module sparsealg
       

    
    
    
    
  
