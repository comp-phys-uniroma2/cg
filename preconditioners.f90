module preconditioners
  use precision
  use matdef
  implicit none
  private

  public :: getULD_jacobi
  public :: getULD_ssor
  public :: getULD_ilu0
  public :: solveLDU
  public :: solveL
  public :: solveU

  interface solveU
    module procedure solveU_csr
    module procedure solveU_dns
  end interface solveU     

  interface solveL
    module procedure solveL_csr
    module procedure solveL_dns
  end interface solveL  

  contains
  ! -----------------------------------------------------------------
  ! Solve the linear system LDUx=b
  subroutine solveLDU(L,D,U,b,x)
    type(rCSR), intent(in) :: L, U
    real(dp), intent(in) :: D(:), b(:)
    real(dp), intent(inout) :: x(:)
  
    real(dp), allocatable :: s(:)
    integer :: i

    allocate(s(size(x)))

    call solveL(L,b,s)
    do i = 1, size(x)
      s(i) = D(i)*s(i)
    end do  
    call solveU(U,s,x) 

    deallocate(s)

  end subroutine solveLDU 

  ! -----------------------------------------------------------------
  ! Solve a lower triangular linear system by substitution
  subroutine solveL_csr(L, b, x)
    type(rCSR), intent(in) :: L
    real(dp), intent(in) :: b(:)
    real(dp), intent(inout) :: x(:)

    integer :: i,j,k,N
    real(dp) :: lij, lii

    N = size(b)
    if (L%nzval(1) == 0) then
      stop 'error l11=0'
    end if
    x(1) = b(1) / L%nzval(1)

    do i = 2, N
      x(i) = b(i)
      do j = L%rowpnt(i), L%rowpnt(i+1)-1
        k = L%colind(j)
        lij = L%nzval(j)
        if (k == i) then
          lii=lij
        else
          x(i) = x(i) - lij * x(k)  
        end if
      end do
      if (lii == 0) then
        print*, 'error lii=0; i=',i
        stop
      end if
      x(i) = x(i) / lii
    end do

  end subroutine solveL_csr

  ! -----------------------------------------------------------------
  ! Solve a lower triangular linear system by substitution
  subroutine solveL_dns(L, b, x)
    real(dp), intent(in) :: L(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(inout) :: x(:)

    integer :: i,j,k,N

    N = size(b)
    if (L(1,1) == 0) then
      stop 'error l11=0'
    end if
    x(1) = b(1) / L(1,1)

    do i = 2, N
      x(i) = b(i)
      do j = 1, i-1
        x(i) = x(i) - L(i,j) * x(j)  
      end do
      if (L(i,i) == 0) then
        print*, 'error lii=0; i=',i
        stop
      end if
      x(i) = x(i) / L(i,i)
    end do

  end subroutine solveL_dns
  ! -----------------------------------------------------------------
  ! Solve a upper triangular linear system by back-substitution
  subroutine solveU_csr(U, b, x)
    type(rCSR), intent(in) :: U
    real(dp), intent(in) :: b(:)
    real(dp), intent(inout) :: x(:)

    integer :: i,j,k,N
    real(dp) :: uij, uii

    N = size(b)

    x(N) = b(N) / U%nzval(size(U%nzval))

    do i = N-1, 1, -1
      x(i) = b(i)
      do j = U%rowpnt(i), U%rowpnt(i+1)-1
        k = U%colind(j)
        uij = U%nzval(j)
        if (k == i) then
          uii=uij
        else
          x(i) = x(i) - uij * x(k)  
        end if
      end do
      if (uii == 0) then
        print*, 'error uii=0; i=',i
        stop
      end if
      x(i) = x(i) / uii
    end do

  end subroutine solveU_csr

  ! -----------------------------------------------------------------
  ! Solve a upper triangular linear system by back-substitution
  subroutine solveU_dns(U, b, x)
    real(dp), intent(in) :: U(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(inout) :: x(:)

    integer :: i,j,k,N
    real(dp) :: uij, uii

    N = size(b)

    x(N) = b(N) / U(N,N)

    do i = N-1, 1, -1
      x(i) = b(i)
      do j = i+1, N
        x(i) = x(i) - U(i,j) * x(j)  
      end do
      if (U(i,i) == 0) then
        print*, 'error uii=0; i=',i
        stop
      end if
      x(i) = x(i) / U(i,i)
    end do

  end subroutine solveU_dns 



  ! JACOBI Preconditioner (diagonal of A)
  ! M_j = D(A) 
  ! Can be implemented either:
  ! L=I; D= A^-1; U=I
  ! or
  ! L=D(A); D=D(A); U=D(A) 
  subroutine getULD_jacobi(A, w, U, L, D)
    type(rCSR), intent(in) :: A   
    real(dp), intent(in) :: w 
    type(rCSR), intent(inout) :: U
    type(rCSR), intent(inout) :: L
    real(dp), intent(inout) :: D(:) 

    integer :: N, i, j

    N=A%nrow
    ! create L=Id U=Id and D=A^-1
    call create(U, N, N)
    call create(L, N, N)
    do i = 1, N
       U%rowpnt(i) = i
       L%rowpnt(i) = i
       U%nzval(i) = 1.0_dp 
       L%nzval(i) = 1.0_dp 
       U%colind(i) = i
       L%colind(i) = i
       do j = A%rowpnt(i), A%rowpnt(i+1) - 1
         if (A%colind(j) == i) then
           D(i) = 1.0_dp/A%nzval(j)
         end if
       end do
    end do
    U%rowpnt(N+1) = N+1
    L%rowpnt(N+1) = N+1

  end subroutine getULD_jacobi

  ! SSOR preconditioner
  ! M_ssor = w/(2-w) [L+D/w] D^-1 [U+D/w]
  ! extract U = U+D/w; L = L+D/w and D matrices
  subroutine getULD_ssor(A, w, U, L, D)
    type(rCSR), intent(in) :: A   
    real(dp), intent(in) :: w 
    type(rCSR), intent(inout) :: U
    type(rCSR), intent(inout) :: L
    real(dp), intent(inout) :: D(:) 

    integer :: i,j,kl,ku,N,colindj

    N=A%nrow
    call create(U, N, (A%nnz-N)/2 + N)
    call create(L, N, (A%nnz-N)/2 + N)
    U%rowpnt(1) = 1
    L%rowpnt(1) = 1
    kl = 0; ku = 0;
    do i = 1, N
      do j = A%rowpnt(i), A%rowpnt(i+1) - 1
         colindj = A%colind(j)
         if (colindj > i) then
           ku = ku + 1    
           U%nzval(ku) = A%nzval(j)
           U%colind(ku) = colindj
         end if
         if (colindj == i) then
           D(i) = (2.0_dp-w) * A%nzval(j)/w
           ku = ku + 1    
           U%nzval(ku) = A%nzval(j)/w
           U%colind(ku) = colindj
           kl = kl + 1    
           L%nzval(kl) = A%nzval(j)/w
           L%colind(kl) = colindj
         end if
         if (colindj < i) then
           kl = kl + 1    
           L%nzval(kl) = A%nzval(j)
           L%colind(kl) = colindj
         end if 
      end do
      L%rowpnt(i+1) = kl + 1
      U%rowpnt(i+1) = ku + 1
    end do
         
  end subroutine getULD_ssor

  ! ILU0 preconditioner
  subroutine getULD_ilu0(A, U, L, D)
    type(rCSR), intent(in) :: A   
    type(rCSR), intent(inout) :: U
    type(rCSR), intent(inout) :: L
    real(dp), intent(inout) :: D(:) 

    type(rMSR) :: LU
    integer :: i, j, N, kl, ku, colindj

    call create(LU, A%nrow, A%nnz+1)

    call ilu0(A, LU)

    N = A%nrow
    do i = 1, N
      D(i) = 1.0_dp/LU%nzval(i)
    end do

    call create(U, N, (A%nnz-N)/2 + N)
    call create(L, N, (A%nnz-N)/2 + N)

    U%rowpnt(1) = 1
    L%rowpnt(1) = 1
    kl = 0; ku = 0;
    do i = 1, N
      do j = LU%colind(i), LU%colind(i+1) - 1
         colindj = LU%colind(j)
         if (colindj > i) then
           ku = ku + 1    
           U%nzval(ku) = LU%nzval(j)
           U%colind(ku) = colindj
         end if
         if (colindj < i) then
           kl = kl + 1    
           L%nzval(kl) = LU%nzval(j)
           L%colind(kl) = colindj
         end if 
      end do
      kl = kl + 1    
      L%nzval(kl) = D(i)
      L%colind(kl) = i
      L%rowpnt(i+1) = kl + 1
      ku = ku + 1
      U%nzval(ku) = D(i)
      U%colind(ku) = i
      U%rowpnt(i+1) = ku + 1
    end do

  end subroutine getULD_ilu0

  !----------------------------------------------------------------------
  ! ILU0 adapted from SPARSKIT
  !
  ! IMPORTANT: The subroutine assumes that A has sorted indices 
  !----------------------------------------------------------------------
  subroutine ilu0(A, LU) 
    type(rCSR) :: A
    type(rMSR) :: LU
   
    integer, allocatable :: iw(:), ju(:)
    integer ::  n, ju0, i, ii, j, jj, jrow, jcol, js, jf, jm, jw
    real(dp) :: tl
   
    n = A%nrow
    allocate(iw(n))
    allocate(ju(n))
    iw = 0
    ju0 = n+2
    LU%colind(1) = ju0
   
    do ii = 1, n
       js = ju0
       do j = A%rowpnt(ii), A%rowpnt(ii+1)-1
          ! copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
          jcol = A%colind(j)
          if (jcol .eq. ii) then
             LU%nzval(ii) = A%nzval(j)
             iw(jcol) = ii
             ju(ii)  = ju0 
          else
             LU%nzval(ju0) = A%nzval(j)
             LU%colind(ju0) = A%colind(j)
             iw(jcol) = ju0
             ju0 = ju0+1
          endif
       end do 
       LU%colind(ii+1) = ju0 !MSR: colind(1:n+1) are the rowpnt
       jf = ju0-1
       jm = ju(ii)-1
            
       do j = js, jm
          jrow = LU%colind(j)
          tl = LU%nzval(j) * LU%nzval(jrow)
          LU%nzval(j) = tl
          ! perform  linear combination
          do jj = ju(jrow), LU%colind(jrow+1)-1
             jw = iw(LU%colind(jj))
             if (jw .ne. 0) then
                LU%nzval(jw) = LU%nzval(jw) - tl*LU%nzval(jj)
             end if   
          end do 
       end do
   
       !  invert  and store diagonal element.
   
       if (LU%nzval(ii) .eq. 0.0d0) then
          print*,'ERROR: zero pivot at row',ii
          stop
       end if    
       LU%nzval(ii) = 1.0_dp/LU%nzval(ii)
       iw(ii) = 0
       do i = js, jf
          iw(LU%colind(i)) = 0
       end do   
    end do 
   
    deallocate(iw, ju)

  end subroutine ilu0

end module preconditioners  
