module gmres_solver
  use precision
  use matdef
  use sparsealg
  implicit none
  
  public :: gmres, allocate_gmres_stuff 
  
  real(dp), dimension(:), allocatable      :: r, cs, sn, e1, beta_gmres, error_gmres, y_gmres, tmp_v
  real(dp), dimension(:,:), allocatable    :: mQ, mH, temp_var, tmp_m
  
  contains

  subroutine allocate_gmres_stuff(vecsize, max_iter)
    integer, intent(in) :: vecsize
    integer, intent(in) :: max_iter 
        
    allocate(          r(vecsize))
    allocate(      tmp_v(vecsize))
    allocate(         cs(max_iter))
    allocate(         sn(max_iter))
    allocate(error_gmres(max_iter))
    allocate(    y_gmres(max_iter))
    allocate(         e1(max_iter+1))
    allocate( beta_gmres(max_iter+1))
    allocate(   temp_var(max_iter+1,2))
    allocate(         mQ(vecsize,max_iter))
    allocate(         mH(vecsize,max_iter))
   
    print*,"Mem: ",(size(mQ)+size(mH))*8," bytes" 
  end subroutine
  
  
  ! --------------------------------------------------------------------
  ! Generalized Minimal Residual
  subroutine gmres(A_in, b_in, x_out, Iter_in, threshold_in, error)
    
    real(dp), dimension(:), intent(out) :: x_out
    real(dp), dimension(:), intent(in)  :: b_in
    type(rCSR), intent(in)              :: A_in
    integer, intent(in)                 :: Iter_in
    real(dp), intent(in)                :: threshold_in
    real(dp), intent(inout)             :: error


    integer                             :: k, i, info 
    real(dp)                            :: b_norm, r_norm
    integer, allocatable                :: ipv(:)
    
    
    print *, "======================="
    print *, "starting GMRES solver  "
    print *, "======================="
    print *,""
    
    cs = 0.0_dp
    e1 = 0.0_dp
    e1(1) = 1.0_dp
    
    call matvec(A_in,x_out,tmp_v)
    
    r = b_in - tmp_v
    
    b_norm = norm2(b_in)
    
    error_gmres = norm2(r) / b_norm
    
    r_norm       = norm2(r)
    mQ(:,1)      = r/r_norm
    beta_gmres(1) = r_norm
    
    mH = 0.0_dp
    mQ = 0.0_dp
    
    print*, "start iterations and error minimization"
    
    do k = 1, Iter_in
      
      call arnoldi(A_in, mQ, k, temp_var)
      
      print*,k,maxval(abs(temp_var))

      mH(1:k+1,k)   = temp_var(:,1)
      mQ(    :,k+1) = temp_var(:,2)
      
      
      ! eliminate the last element in H ith row and update the rotation matrix
      call apply_givens_rotation(mH(1:k+1,k), cs, sn, k , cs(k), sn(k))
     
      print*,'cs sn:',cs(k), sn(k) 
      
      ! Residual vector update
      beta_gmres(k+1) = -sn(k)*beta_gmres(k)
      beta_gmres(k)   =  cs(k)*beta_gmres(k)
      error_gmres(k+1) = abs(beta_gmres(k+1))/b_norm
    
      error = error_gmres(k+1)
      print*,'error:',error 
      if (error <= threshold_in) then
        exit
      end if
      
    end do
    
    print *, "done"
    print *,""
    
    print *, "solving system with DGESV"
    
    allocate(tmp_m(Iter_in, Iter_in))
    allocate(ipv(Iter_in))
    tmp_m = mH(1:Iter_in,1:Iter_in)
    y_gmres = beta_gmres
    call dgesv(Iter_in, 1, tmp_m, Iter_in, ipv, y_gmres, Iter_in, info) 
    
    if (info /= 0) then
       print*, "Error in solver: ",info
       error = info
       return 
    end if 

    deallocate(ipv)
    deallocate(tmp_m)

    print *, "done"
    print *, ""

    print *, "outputting result"
    
    
    !$OMP PARALLEL, public(mQ, y_gmres, x_out, Iter_in), private(i)
    !$OMP DO
    do i = 1, size(x_out)
      x_out(i)= x_out(i) + dot_product(mQ(i,1:Iter_in), y_gmres)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    print *, "done!"
    print *, ""
    
  end subroutine
  
  subroutine arnoldi(A_in,Q_in,k_in, Mout)
    
    real(dp), dimension(:,:), intent(out)  :: Mout
    real(dp), dimension(:,:), intent(in)   :: Q_in
    type(rCSR), intent(in)                 :: A_in
    integer, intent(in)                    :: k_in
    real(dp), dimension(:), allocatable    :: q, h
    integer                                :: i
    integer, dimension(0:1)                :: Qsize
    
    Qsize = shape(Q_in)
    print*,'(arnoldi) Qsize:',Qsize 
    allocate(q(Qsize(0)))
    allocate(h(Qsize(0)))
    
    ! Krylov vector
    
    call matvec(A_in, Q_in(:,k_in), q)
    
    ! Modified Gram-Schmidt, keeping the Hessenberg matrix
    do i = 1, k_in
      h(i) = dot_product(q, Q_in(:,i))
      q    = q - h(i)*Q_in(:,i)
    end do
    
    h(k_in) = norm2(q)
    q = q/norm2(q)
    
    Mout(:,1) = h
    Mout(:,2) = q
    
    deallocate(q)
    deallocate(h)
    
  end subroutine

  subroutine givens_rotation(v1,v2,cs_out,sn_out)
  
    real(dp)              :: t
    real(dp), intent(out) :: cs_out,  sn_out
    real(dp), intent(in)  :: v1,v2
    
    t=sqrt(v1**2 + v2**2)
    
    cs_out = v1/t
    sn_out = v2/t
    
  end subroutine
  
  subroutine apply_givens_rotation(h_io, cs_in, sn_in, k_in, cs_k, sn_k)
    
    real(dp), dimension(:), intent(out) :: h_io
    real(dp), intent(out)               :: cs_k,  sn_k
    real(dp), dimension(:), intent(in)  :: cs_in, sn_in
    integer, intent(in)                 :: k_in
    real(dp)                            :: temp
    integer                             :: i
    
    !applied to i-th column
    
    do i = 1, k_in-1
      
      temp      =  cs_in(i)*h_io(i) + sn_in(i)*h_io(i + 1)
      h_io(i+1) = -sn_in(i)*h_io(i) + cs_in(i)*h_io(i + 1)
      h_io(i)   = temp
      
    end do
    
    ! update the next sin cos values for rotation
    call givens_rotation(h_io(k_in-1),h_io(k_in), cs_k, sn_k)
    
    ! eliminate H(i,i-1)
    h_io(k_in) = cs_k*h_io(k_in) + sn_k*h_io(k_in+1)
    h_io(k_in+1)   = 0.0_dp
    
  end subroutine
 

 
end module
