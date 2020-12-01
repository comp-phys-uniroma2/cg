module adj_map
  implicit none
  private

  public :: adjlist
  public :: create_2D
  public :: create_3D
  public :: count_nnz

  type adjlist
    integer, allocatable :: neignode(:)
  end type adjlist

  contains

  subroutine create_2D(graph,N,M)
    type(adjlist), allocatable :: graph(:)
    integer, intent(in) :: N, M

    integer :: i, j, node, nneig

    allocate(graph(N*M))

    do i = 1, N
      do j = 1, M
         node = coo2node(i,j) 
         nneig = 4 
         if (i==1 .or. i==N) nneig = nneig - 1
         if (j==1 .or. j==M) nneig = nneig - 1
         allocate(graph(node)%neignode(nneig))
         nneig = 0
         if (i>1) then
           nneig = nneig +1
           graph(node)%neignode(nneig) = coo2node(i-1,j)
         end if  
         if (i<N) then 
           nneig = nneig +1
           graph(node)%neignode(nneig) = coo2node(i+1,j)
         end if 
         if (j>1) then 
           nneig = nneig +1
           graph(node)%neignode(nneig) = coo2node(i, j-1)  
         end if  
         if (j<M) then
           nneig = nneig +1
           graph(node)%neignode(nneig) = coo2node(i, j+1)  
         end if  
      end do
    end do
  
    contains

    function coo2node(i,j) result(node)
      integer, intent(in) :: i, j
      integer :: node

      node = (j-1)*N + i

    end function coo2node

  end subroutine create_2D

  subroutine create_3D(graph, N, M, P)
    type(adjlist), intent(inout), allocatable :: graph(:)
    integer, intent(in) :: N, M, P

    integer :: i, j, k, node, nneig, l

    allocate(graph(N*M*P))

    do i = 1, N
       do j = 1, M
          do k = 1, P
             node = coo2node(i,j,k)
             nneig = 6
             if (i==1 .or. i==N) nneig = nneig - 1
             if (j==1 .or. j==M) nneig = nneig - 1            
             if (k==1 .or. k==P) nneig = nneig - 1
             allocate(graph(node)%neignode(nneig))
             l = 0
             if (i > 1) then
                l = l + 1
                graph(node)%neignode(l) = coo2node(i-1,j,k)
             end if
             
             if (i < N) then
                l = l + 1
                graph(node)%neignode(l) = coo2node(i+1,j,k)
             end if
             
             if (j > 1) then
                l = l + 1
                graph(node)%neignode(l) = coo2node(i,j-1,k)
             end if
             
             if (j < M) then
                l = l + 1
                graph(node)%neignode(l) = coo2node(i,j+1,k)
             end if
             
              if (k > 1) then
                l = l + 1
                graph(node)%neignode(l) = coo2node(i,j,k-1)
             end if
             
             if (k < P) then
                l = l + 1
                graph(node)%neignode(l) = coo2node(i,j,k+1)
             end if
          end do
       end do
    end do

  contains

    function coo2node(i,j,k) result(node)
      integer, intent(in) :: i,j,k
      integer :: node

      node = (k-1)*N*M + (j-1)*N + i

    end function coo2node
                         
  end subroutine create_3D
             
  function count_nnz(graph) result(nnz)
    type(adjlist), intent(in) :: graph(:)
    integer :: nnz

    integer :: i

    nnz = 0
    do i = 1, size(graph)
       nnz = nnz + size(graph(i)%neignode) + 1
    end do   

  end function count_nnz
    
end module adj_map
