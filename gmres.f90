module gmres_solver
  use precision
  use matdef
  use sparsealg
  implicit none
  
  public :: gmres, allocate_gmres_stuff 
  
  real(dp), dimension(:,:), allocatable    :: mQ, mH 
  
  contains

  subroutine allocate_gmres_stuff(vecsize, max_iter)
    integer, intent(in) :: vecsize
    integer, intent(in) :: max_iter 
        
    allocate(         mQ(vecsize,max_iter))
    allocate(         mH(vecsize,max_iter))
   
    print*,"Mem: ",(size(mQ)+size(mH))*8," bytes" 
  end subroutine
  
  
  ! --------------------------------------------------------------------
  ! Generalized Minimal Residual
  subroutine gmres(A_in, b_in, x_out, Iter_in, threshold_in, error)
    
    real(dp), dimension(:), intent(inout) :: x_out
    real(dp), dimension(:), intent(in)  :: b_in
    type(rCSR), intent(in)              :: A_in
    integer, intent(in)                 :: Iter_in
    real(dp), intent(in)                :: threshold_in
    real(dp), intent(inout)             :: error


    integer                             :: k, i, info 
    real(dp)                            :: b_norm, r_norm
    integer, allocatable                :: ipv(:)
    real(dp), allocatable               :: beta(:), tmp_v(:), r(:)
    real(dp), allocatable               :: tmp_m(:,:) 
    real(dp), allocatable               :: cs(:), sn(:), y_gmres(:)
    
    allocate(beta(Iter_in+1))
    allocate(cs(Iter_in))
    allocate(sn(Iter_in))
    allocate(tmp_v(size(mQ,1)))
    allocate(r(size(mQ,1)))

    print *, "======================="
    print *, "starting GMRES solver  "
    print *, "======================="
    print *,""
    
    mH = 0.0_dp
    mQ = 0.0_dp
    
    cs = 0.0_dp
    sn = 0.0_dp
    beta = 0.0_dp
    
    call matvec(A_in,x_out,tmp_v)
    r = b_in - tmp_v
    b_norm = norm2(b_in)
    b_norm = 1.0_dp
    r_norm = norm2(r)
    error = r_norm / b_norm
    print*,'residual error:', error

    do while (error > threshold_in)

       mQ(:,1)      = r/r_norm
       beta(1)      = r_norm
  
       do k = 1, Iter_in-1
        
         ! operates on mH and mQ 
         call arnoldi(A_in, k)
         
         ! eliminate the last element in H ith row and update the rotation matrix
         call apply_givens_rotation(cs, sn, k)
        
         !if (isnan(cs(k))) stop 
         ! Residual vector update
         beta(k+1) = -sn(k)*beta(k)
         beta(k)   =  cs(k)*beta(k)
         error = abs(beta(k+1))/b_norm
       
         print*,k, 'error:', error 
         if (error <= threshold_in) then
           exit
         end if
         
       end do
 
       if (error >= threshold_in) then 
          k = k - 1
       end if

       print *, "done"
       print *,""
       print *, "solving system with DGESV"
       
       allocate(tmp_m(k, k))
       allocate(ipv(k))
       allocate(y_gmres(k))
       tmp_m = mH(1:k,1:k)
       y_gmres = beta(1:k)
       call dgesv(k, 1, tmp_m, k, ipv, y_gmres, k, info) 
       
       if (info /= 0) then
          print*, "Error in solver: ",info
          error = info
          return 
       end if 
  
       ! !$OMP PARALLEL, public(mQ, y_gmres, x_out, Iter_in), private(i)
       ! !$OMP DO
       do i = 1, size(x_out)
         x_out(i)= x_out(i) + dot_product(mQ(i,1:k), y_gmres(1:k)) 
       end do
       ! !$OMP END DO
       ! !$OMP END PARALLEL
       deallocate(tmp_m)
       deallocate(ipv)
       deallocate(y_gmres)
      
       ! cross check the residual: 
       call matvec(A_in,x_out,tmp_v)
       r = b_in - tmp_v
       r_norm = norm2(r)
       error = r_norm / b_norm
       print*,'residual error:', error
    
    end do   

    deallocate(sn,cs)
    deallocate(beta)
    deallocate(tmp_v)
    deallocate(r)
    print *, "done!"
    print *, ""
    
  end subroutine


  ! --------------------------------------------------------------------- 
  subroutine arnoldi(A_in, k_in)
    
    type(rCSR), intent(in)                 :: A_in
    integer, intent(in)                    :: k_in
    
    real(dp), dimension(:), allocatable    :: q
    real(dp)                               :: q_norm
    integer                                :: i
    integer                                :: Qsize
    
    Qsize = size(mQ, 1)
    allocate(q(Qsize))
   
    ! new Krylov vector 
    call matvec(A_in, mQ(:,k_in), q)
    
    ! Modified Gram-Schmidt, keeping the Hessenberg matrix
    do i = 1, k_in
      mH(i, k_in) = dot_product(q(:), mQ(:,i))
      q(:)    = q(:) - mH(i, k_in)*mQ(:,i)
    end do
    
    q_norm = norm2(q)
    mH(k_in+1, k_in) = q_norm
    mQ(:, k_in+1) = q(:) / q_norm
    
    deallocate(q)
    
  end subroutine
  ! ---------------------------------------------------------------------
  subroutine givens_rotation(v1,v2,cs_out,sn_out)
  
    real(dp), intent(in)  :: v1,v2
    real(dp), intent(out) :: cs_out,  sn_out
    
    real(dp)              :: t
    t=sqrt(v1**2 + v2**2)
    
    cs_out = v1/t
    sn_out = v2/t
    
  end subroutine
  
  ! ---------------------------------------------------------------------
  subroutine apply_givens_rotation(cs, sn, k)
    
    real(dp), dimension(:), intent(inout)  :: cs, sn
    integer, intent(in)                 :: k
    
    real(dp)               :: temp
    integer                :: i

    !applied to i-th column    
    do i = 1, k-1
      temp   =  cs(i)*mH(i,k) + sn(i)*mH(i+1,k)
      mH(i+1,k) = -sn(i)*mH(i,k) + cs(i)*mH(i+1,k)
      mH(i,k)   = temp
    end do
    
    ! update the next sin cos values for rotation
    call givens_rotation(mH(k,k), mH(k+1,k), cs(k), sn(k))
    
    ! eliminate H(i,i-1)
    mH(k,k)   = cs(k)*mH(k,k) + sn(k)*mH(k+1,k)
    mH(k+1,k) = 0.0_dp
    
  end subroutine
  ! ---------------------------------------------------------------------

 
end module
