module cg
  use precision
  use matdef
  use preconditioners
  use sparsealg
  implicit none
  private

  public :: conjgrads
  public :: pconjgrads

contains

  subroutine conjgrads(A,b,x0,x,tol)
    type(rCSR), intent(in) :: A
    real(dp), intent(in) :: b(:)
    real(dp), intent(in) :: x0(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(inout) :: tol

    real(dp), allocatable :: r(:), y(:), d(:), d0(:)
    integer :: k, maxiter
    real(dp) :: alfa, beta, norm1, norm2
    
    maxiter = size(b)
    allocate(r(size(b)))
    allocate(y(size(b)))
    allocate(d(size(b)))
    allocate(d0(size(b)))
    
    call matvec(A,x0,y)

    r = b - y
    d = r
    d0= d
    x = x0
    norm1 = dot_product(r,r)

    do k = 1, maxiter
       call matvec(A,d,y)
       alfa = norm1/dot_product(d,y)
       x = x + alfa*d
       r = r - alfa*y
       norm2 = dot_product(r,r)
       write(*,'(a,I5,a,E10.3,a,E10.3)') 'CG iteration:',k,'  err:',sqrt(norm2),'  ort:',dot_product(d0,y)
       if (sqrt(norm2) < tol) then
          exit
       end if
       beta = norm2/norm1
       d = r + beta*d
       norm1 = norm2
    end do

    tol = sqrt(norm2)
    
  end subroutine conjgrads
 
  ! Preconditioned CG 
  subroutine pconjgrads(A,b,x0,x,ip,w,tol)

    type(rCSR), intent(in) :: A
    real(dp), intent(in) :: b(:)
    real(dp), intent(in) :: x0(:)
    real(dp), intent(inout) :: x(:)
    character, intent(in) :: ip
    real(dp), intent(in) :: w 
    real(dp), intent(inout) :: tol

    real(dp), allocatable :: r(:), y(:), p(:), D(:)
    type(rCSR) :: U, L
    integer :: i, j, k, maxiter, N
    real(dp) :: alfa, beta, norm1, norm2
    
    N = size(b)
    maxiter = N
    allocate(r(N))
    allocate(y(N))
    allocate(p(N))
    allocate(D(N))
    
    call matvec(A,x0,y)
    r = b - y

    ! select preconditioner
    select case(ip)
    case('J')
      call getULD_jacobi(A, w, U, L, D)
    case('S')
      call getULD_ssor(A, w, U, L, D)
    case('I')
      call getULD_ilu0(A, U, L, D)
    end select    
    
    ! p = M^-1 r
    ! Solve System: M p = r
    ! => solve [L+D/w] y = r, solve [U+D/w] p = (2-w)/w D y = s
    call solveLDU(L,D,U,r,p)

    x = x0
    norm1 = dot_product(r,p)

    do k = 1, maxiter
       call matvec(A,p,y)
       alfa = norm1/dot_product(p,y)
       x = x + alfa*p
       r = r - alfa*y
       ! Solve System: M y = r
       call solveLDU(L,D,U,r,y)
       
       norm2 = dot_product(r,y)

       write(*,'(a,I5,a,E10.3)') 'CG iteration:',k,'  err:',sqrt(abs(norm2))
       if (sqrt(abs(norm2)) < tol) then
          exit
       end if
       beta = norm2/norm1
       p = y + beta*p
       norm1 = norm2
    end do

    tol = sqrt(norm2)
    
  end subroutine pconjgrads

end module cg
