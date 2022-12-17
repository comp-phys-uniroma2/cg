program poisson
  use precision
  use sparsealg
  use matdef
  use adj_map
  use cg
  implicit none

  integer :: N 
  type(adjlist), allocatable :: graph(:)
  real(dp), allocatable :: phi(:), rhs(:)
  real(dp) :: L, dx, tol, w
  real(dp), parameter :: Pi = 3.141593_dp
  integer :: node, nnodes, nnz, i, j, k, nn, fd
  type(rCSR) :: A_csr
  character(64) :: arg
  character(1) :: PC ! preconditioner type
  
  !|---+---+---+---|
  !1   2   3   4   5

  if (iargc()<2) then
    write(*,*) 'poisson N tol [Precond w]'
    write(*,*) 'N: punti griglia 3d NxNxN'    
    write(*,*) 'Precond: J:=Jacobi; S=SSOR; I=ILU0'
    stop
  end if

  call getarg(1,arg)
  read(arg,*) N
  call getarg(2,arg)
  read(arg,*) tol
  ! Read optional preconditioner
  if (iargc()>2) then
    call getarg(3,arg)
    read(arg,*) PC
    select case(PC)
    case('J')
    case('S')
      if (iargc()<4) then 
        stop 'SSOR requires the weight parameter'
      end if    
      call getarg(4,arg)
      read(arg,*) w
    case('I')
    case default
      stop 'invalid preconditioner'
    end select    
  end if

  L = 1.0_dp
  dx = L/real(N-1,dp)

  nnodes = N*N*N
  
  call create_3D(graph, N, N, N)

  nnz = count_nnz(graph) 

  call create(A_csr, nnodes, nnz)

  print*,'Init Matrices'
  k = 0
  A_csr%rowpnt(1) = 1
  do i = 1, nnodes
    k = k + 1
    A_csr%nzval(k) = -6.0_dp
    A_csr%colind(k) = i
    nn = size(graph(i)%neignode)
    do j = 1, nn
       k = k + 1
       A_csr%nzval(k) = 1.0_dp
       A_csr%colind(k) = graph(i)%neignode(j)
    end do
    A_csr%rowpnt(i+1) = k + 1 
 end do

 print*,k,nnz

 allocate(rhs(nnodes))
 allocate(phi(nnodes))

 rhs = 0.0_dp
 phi = 0.0_dp

 print*, nint(L/10/dx)
 rhs(coo2node(N/2, N/2, N/2+nint(L/5/dx))) = -4.0_dp * Pi * 14.4 / dx
 rhs(coo2node(N/2, N/2, N/2+nint(L/5/dx))) = 3.0_dp *4.0_dp * Pi * 14.4 / dx 
 rhs(coo2node(N/2, N/2-nint(L/5/dx),N/2)) = -4.0_dp * Pi * 14.4 / dx
 rhs(coo2node(N/2, N/2-nint(L/5/dx),N/2)) = -4.0_dp * Pi * 14.4 / dx 

 print*,'Solve'

 select case(PC)
 case('J', 'S', 'I')
   call pconjgrads(A_csr, rhs, phi, phi, PC, w, tol) 
 case default
   call conjgrads(A_csr, rhs, phi, phi, tol) 
 end select 
 
 open(newunit=fd, file='sol.dat')

 do i = 1, N
    do j = 1, N
       write(fd,*) phi(coo2node(N/2,i,j)) 
    end do
    write(fd,*)
 end do

 close(fd)
    
contains

  function coo2node(i,j,k) result(node)
    integer, intent(in) :: i,j,k
    integer :: node

   node = (k-1)*N*N + (j-1)*N + i

  end function coo2node
 
 
 
end program poisson
 
 
