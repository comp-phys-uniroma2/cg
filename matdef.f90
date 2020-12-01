module matdef
  use precision
  implicit none
  private

  public :: rCSR, rCOO, rMSR
  public :: create
  
  type rCSR
     integer :: nrow
     integer :: nnz    
     real(dp), allocatable :: nzval(:)
     integer, allocatable :: colind(:)
     integer, allocatable :: rowpnt(:)    
  end type rCSR
     
  type rCOO
     integer :: nrow
     integer :: nnz
     real(dp), allocatable :: nzval(:)
     integer, allocatable :: colind(:)
     integer, allocatable :: rowind(:)     
  end type rCOO
 
  ! Modified Compressed Sparse Raw 
  ! Diagonal on elements nzval(1:n)
  ! rowpnt on colind(1:n+1)
  type rMSR
     integer :: nrow
     integer :: nnz    
     real(dp), allocatable :: nzval(:)
     integer, allocatable :: colind(:)
  end type rMSR

  interface create
     module procedure create_rCSR
     module procedure create_rCOO
     module procedure create_rMSR
  end interface create
  

contains

  subroutine create_rCSR(M, nrow, nnz)
    type(rCSR), intent(inout) :: M
    integer, intent(in) :: nrow
    integer, intent(in) :: nnz

    M%nrow = nrow
    M%nnz = nnz
    
    allocate(M%nzval(nnz))
    allocate(M%colind(nnz))
    allocate(M%rowpnt(nrow+1))
    
  end subroutine create_rCSR

  subroutine create_rCOO(M, nrow, nnz)
    type(rCOO), intent(inout) :: M
    integer, intent(in) :: nrow
    integer, intent(in) :: nnz

    M%nrow = nrow
    M%nnz = nnz
    
    allocate(M%nzval(nnz))
    allocate(M%colind(nnz))
    allocate(M%rowind(nnz))
    
  end subroutine create_rCOO
  
  ! nnz: ONLY the non-zero off-diagonal entries
  subroutine create_rMSR(M, nrow, nnz)
    type(rMSR), intent(inout) :: M
    integer, intent(in) :: nrow
    integer, intent(in) :: nnz

    M%nrow = nrow
    M%nnz = nrow+nnz+1
    
    allocate(M%nzval(nrow+nnz+1))
    allocate(M%colind(nrow+nnz+1))
    
  end subroutine create_rMSR


!-----------------------------------------------------------------------
! This routine sorts the elements of  a matrix (stored in Compressed
! Sparse Row Format) in increasing order of their column indices within 
! each row. It uses a form of bucket sort with a cost of O(nnz) where
! nnz = number of nonzero elements. 
! requires an integer work array of length 2*nnz.  
!-----------------------------------------------------------------------
  subroutine sort_csr(A,values)
    type(rCSR), intent(inout) :: A
    logical, intent(in) :: values    

    integer :: i,j,k,n, nnz, next, ifirst, ko, krow, irow
    integer, allocatable :: iwork(:)

    n = A%nrow
    allocate(iwork(2*A%nnz))

    do i = 1, n+1
         iwork(i) = 0
    end do
    do i=1, n
      do k=A%rowpnt(i), A%rowpnt(i+1)-1 
         j = A%colind(k)+1
         iwork(j) = iwork(j)+1
      end do 
    end do
    
    ! compute pointers from lengths. 
    iwork(1) = 1
    do i=1,n
      iwork(i+1) = iwork(i) + iwork(i+1)
    end do

    !get the positions of the nonzero elements in order of columns.
    ifirst = A%rowpnt(1) 
    nnz = A%rowpnt(n+1)-ifirst
    do i = 1, n
       do k=A%rowpnt(i), A%rowpnt(i+1)-1 
          j = A%colind(k)
          next = iwork(j) 
          iwork(nnz+next) = k
          iwork(j) = next+1
       end do
     end do
 
     ! convert to coordinate format
     do i = 1, n
        do k = A%rowpnt(i), A%rowpnt(i+1)-1 
           iwork(k) = i
        end do
     end do

     ! loop to find permutation: for each element find the correct 
     ! position in (sorted) arrays a, ja. Record this in iwork. 
     do k=1, nnz
        ko = iwork(nnz+k) 
        irow = iwork(ko)
        next = A%rowpnt(irow)

        ! the current element should go in next position in row. iwork
        ! records this position. 
        iwork(ko) = next
        A%rowpnt(irow) = next+1
     end do
 
     ! perform an in-place permutation of the  arrays.
     !call ivperm(nnz, A%colind(ifirst), iwork) 
     !if (values) call dvperm(nnz, A%nzval(ifirst), iwork) 

     ! reshift the pointers of the original matrix back.
     do i = n, 1, -1
       A%rowpnt(i+1) = A%rowpnt(i)
     end do
     A%rowpnt(1) = ifirst 
  
  end subroutine sort_csr
      
  !-----------------------------------------------------------------------
  ! this subroutine performs an in-place permutation of an integer vector 
  ! ix according to the permutation array perm(*), i.e., on return, 
  ! the vector x satisfies,
  !
  ! ix(perm(j)) :== ix(j), j=1,2,.., n
  !
  !-----------------------------------------------------------------------
  subroutine ivperm(ix, perm) 
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: ix(:)

    integer :: init, n, ii, j, k, tmp, tmp1, next

    init      = 1
    tmp       = ix(init)  
    ii        = perm(init)
    perm(init)= -perm(init)
    n         = size(ix)
    
    outer: do k = 1, n 
      tmp1    = ix(ii) 
      ix(ii)  = tmp
      next    = perm(ii) 
      if (next .lt. 0) then
        do      
          init = init+1
          if (init .gt. n) exit outer
          if (perm(init) .ge. 0) exit
        end do  
        tmp = ix(init)
        ii  = perm(init)
        perm(init) = -perm(init)
      else
        if (k .gt. n) exit outer
        tmp       = tmp1
        perm(ii)  = -perm(ii)
        ii        = next
      end if   
    end do outer 

    do j=1, n
      perm(j) = -perm(j)
    end do
 
  end subroutine ivperm

  subroutine bsort_csr(A,values)
    type(rCSR), intent(inout) :: A
    logical, intent(in) :: values    

    integer :: i, j, k, js, jf, nrow, itmp
    real(dp) :: dtmp

    do i = 1, nrow
      js = A%rowpnt(i)
      jf = A%rowpnt(i+1)-1
      do j = js, jf
        do k = js+1, jf
          if (A%colind(k) < A%colind(j)) then
             itmp = A%colind(k)
             A%colind(k) = A%colind(j)
             A%colind(j) = itmp
             if (values) then
               dtmp = A%nzval(k)
               A%nzval(k) = A%nzval(j)
               A%nzval(j) = dtmp
             end if  
          end if
        end do
      end do
      !call isort(A%colind, js, jf, perm)
      !call iapply(A%colind, js, jf, perm)
      !call dapply(A%nzval, js, jf, perm)
    end do

  end subroutine bsort_csr  


  


end module matdef
