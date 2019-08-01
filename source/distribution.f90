!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module distribution

!BOP
! !MODULE: distribution
!
! !DESCRIPTION:
!  This module provides data types and routines for distributing
!  blocks across processors.
!
!  Added by Junmin, Feb 1, 2016
!  Port routines necessary to create space-filling curves.   
!
! !REVISION HISTORY:
!  CVS:$Id: distribution.F90,v 1.11 2003/12/23 22:11:40 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use communicate
   use broadcast
   use blocks
   use exit_mod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: & ! 
      nblock_mic_pp = 1 ,&! block# per MIC process/rank
      nproc_cpu_pn = 24 ,&! MPI process# on CPU per node
      nproc_mic_pn = 0 ,&! MPI process# on MIC per node
      nproc_pn = 24     ,&! MPI process# per node
      nnode = 1          ! total node number

   real(r8), parameter, public :: & ! 
      k_cpu_mic = 3.0_r8   ! quick hack on block distribution

   integer , parameter, public :: & ! 
      availf = -100000  ! avail flag for effective non-land block to be allocated

! !PUBLIC TYPES:

   type, public :: distrb  ! distribution data type
      integer (int_kind) :: &
         nprocs            ,&! number of processors in this dist
         communicator        ! communicator to use in this dist

      integer (int_kind), dimension(:), pointer :: &
         proc              ,&! processor location for this block
         local_block         ! block position in local array on proc

      integer (int_kind), dimension(:), pointer :: &
         local_block_ids     ! global id for blocks in local process

      integer (int_kind) :: &
         local_block_num     ! number of blocks for local process
   end type

   ! from POP_SpaceCurveMod, by Junmin, Feb 1, 2016
   type, public :: factor_t
        integer(int_kind)        :: numfact ! The # of factors for a value
        integer(int_kind), dimension(:), pointer :: factors ! The factors
        integer(int_kind), dimension(:), pointer :: used
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution, &
             create_distribution_tiny, &
             create_local_block_ids, &
             create_local_block_ids_tiny, &
             init_mic_proc_flag

   ! Fine-grained opt, (nproc_cpu_pn=0), by Junmin, Apr 28, 2016
   integer (int_kind), public :: & ! 
      loc_nproc_cpu_pn, &! MPI process# on CPU per node
      loc_nproc_mic_pn   ! MPI process# on MIC per node

   ! from POP_SpaceCurveMod, by Junmin, Feb 1, 2016
   public :: GenSpaceCurve,     &
	     IsLoadBalanced

   public :: Factor, 		&
	     IsFactorable, 	&
	     PrintFactor,	&
	     ProdFactor,	&
	     MatchFactor

! !PRIVATE MEMBER FUNCTIONS:

   private :: map,    		&
	      PeanoM, 		&
	      Hilbert, 		&
	      Cinco,  		&
              GenCurve

   private :: FirstFactor, 	&
 	      FindandMark

   integer(int_kind), dimension(:,:), allocatable ::  &
	dir,      &! direction to move along each level
        ordered,  &! the ordering 
	traversal  ! index increments to traverse domain

   integer(int_kind), dimension(:), allocatable ::  &
	pos        ! position along each of the axes
   
   integer(int_kind) ::  &
	maxdim,	  &! dimensionality of entire space
	vcnt       ! visitation count

   logical           :: verbose=.FALSE. 
   
   type (factor_t),  public,save :: fact  ! stores the factorization


!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 function create_distribution(dist_type, nprocs, work_per_block)

! !DESCRIPTION:
!  This routine determines the distribution of blocks across processors
!  by call the appropriate subroutine based on distribution type
!  requested.  Currently only two distributions are supported:
!  2-d Cartesian distribution (cartesian) and a load-balanced
!  distribution (balanced) based on an input amount of work per
!  block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      dist_type             ! method for distributing blocks
                            !  either cartesian or balanced

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distribution   ! resulting structure describing
                            !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('cartesian')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('balanced')

      create_distribution = create_distrb_balanced(nprocs, &
                                                   work_per_block)

   ! ljm tuning
   case('cart_cpuonly')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('cart_micmixed')

      create_distribution = create_distrb_mic(nprocs, work_per_block)
   ! ported from cesm1_2_2, by junmin, Feb 1, 2016
   case('spacecurve')

      create_distribution = create_distrb_spacecurve(nprocs, &
						   work_per_block)

   case default

      call exit_POP(sigAbort,'distribution: unknown distribution type')

   end select

!-----------------------------------------------------------------------
!EOC

 end function create_distribution

!***********************************************************************
!BOP
! !IROUTINE: init_mic_proc_flag
! !INTERFACE:

 subroutine init_mic_proc_flag()

! !DESCRIPTION:
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

   lmic_proc = (my_task >= nproc_cpu_pn * nnode)
   !lmic_proc = (my_task == nproc_cpu_pn * nnode)
   !lmic_proc = .false.
   lmic_trace = .false.
   ! Fine-grained opt, by Junmin, Apr 28,2016
   if (nproc_cpu_pn /= 0) then
      loc_nproc_cpu_pn = nproc_cpu_pn
      loc_nproc_mic_pn = nproc_mic_pn
   else
      loc_nproc_cpu_pn = nproc_mic_pn
      loc_nproc_mic_pn = nproc_cpu_pn
   endif
!EOC

 end subroutine init_mic_proc_flag


!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids
! !INTERFACE:

 subroutine create_local_block_ids(block_ids, distribution)

! !DESCRIPTION:
!  This routine determines which blocks in an input distribution are
!  located on the local processor and creates an array of block ids
!  for all local blocks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which local
                             !  blocks required

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), pointer :: &
      block_ids              ! array of block ids for every block
                             ! that resides on the local processor

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bid, bcount        ! dummy counters

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%proc)
      if (distribution%proc(n) == my_task+1) bcount = bcount + 1
   end do

   !ljm tuning
   write(6,*) 'bcount ', my_task, bcount
   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      do n=1,size(distribution%proc)
         if (distribution%proc(n) == my_task+1) then
            block_ids(distribution%local_block(n)) = n
         endif
      end do
   endif

!EOC

 end subroutine create_local_block_ids

!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids_tiny
! !INTERFACE:

 subroutine create_local_block_ids_tiny(block_ids, distribution)

! !DESCRIPTION:
!  This routine determines which blocks in an input distribution are
!  located on the local processor and creates an array of block ids
!  for all local blocks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which local
                             !  blocks required

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), pointer :: &
      block_ids              ! array of block ids for every block
                             ! that resides on the local processor

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bid, bcount        ! dummy counters

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = nblocks_tot_tiny
   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      do n=1,bcount
         block_ids(n) = n + nblocks_tot ! global id
      end do
   endif

!EOC

 end subroutine create_local_block_ids_tiny


!***********************************************************************
!BOP
! !IROUTINE: create_distrb_cart
! !INTERFACE:

 function create_distrb_cart(nprocs, work_per_block)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs     ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_cart  ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n               ,&! dummy loop indices
      iblock, jblock, nblck ,&!
      is, ie, js, je        ,&! start, end block indices for each proc
      local_block           ,&! block location on this processor
      nprocs_x              ,&! num of procs in x for global domain
      nprocs_y              ,&! num of procs in y for global domain
      nblocks_x_loc         ,&! num of blocks per processor in x
      nblocks_y_loc           ! num of blocks per processor in y

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

   ! ljm tuning
   if (my_task == nproc_cpu_pn*nnode) then ! first MIC process
   !if (my_task == master_task) then
    write(6,*) 'dist_cart: nprocs,nprocs_x,nprocs_y',nprocs,nprocs_x,nprocs_y
    write(6,*) '         : nblocks_tot,nblocks_x,nblocks_y',nblocks_tot,nblocks_x,nblocks_y
   endif
!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   nblocks_x_loc = (nblocks_x-1)/nprocs_x + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_y + 1

   do j=1,nprocs_y
   do i=1,nprocs_x
      n = (j-1)*nprocs_x + i

      is = (i-1)*nblocks_x_loc + 1
      ie =  i   *nblocks_x_loc
      if (ie > nblocks_x) ie = nblocks_x
      js = (j-1)*nblocks_y_loc + 1
      je =  j   *nblocks_y_loc
      if (je > nblocks_y) je = nblocks_y

      local_block = 0
      do jblock = js,je
      do iblock = is,ie
         nblck = (jblock - 1)*nblocks_x + iblock
         if (work_per_block(nblck) /= 0) then
            local_block = local_block + 1
            dist%proc(nblck) = n
            dist%local_block(nblck) = local_block
         else
            dist%proc(nblck) = 0
            dist%local_block(nblck) = 0
         endif
      end do
      end do
   end do
   end do

!----------------------------------------------------------------------

   create_distrb_cart = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_cart

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_spacecurve
! !INTERFACE:

 function create_distrb_spacecurve(nprocs,work_per_block)

! !Description:
!  This function distributes blocks across processors in a 
!  load-balanced manner using space-filling curves
!
! !REVISION HISTORY:
!  added by J. Dennis 3/10/06 

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_spacecurve  ! resulting structure describing
                                ! load-balanced distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,n              ,&! dummy loop indices
      pid                  ,&! dummy for processor id
      local_block          ,&! local block position on processor
      max_work             ,&! max amount of work in any block
      nprocs_x             ,&! num of procs in x for global domain
      nprocs_y               ! num of procs in y for global domain

   integer (int_kind), dimension(:),allocatable :: &
        idxT_i,idxT_j

   integer (int_kind), dimension(:,:),allocatable :: Mesh, Mesh2, Mesh3
   integer (int_kind) :: nblocksL,nblocks,ii,extra,i2,j2,tmp1,s1,ig

   integer (int_kind) :: ierr
   logical, parameter :: Debug = .FALSE.

   integer (int_kind), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      work_tmp           ,&! work per row or column for rake algrthm
      proc_tmp           ,&! temp processor id for rake algrthm
      block_count          ! counter to determine local block indx

   type (distrb) :: dist  ! temp hold distribution

   type (factor_t) :: xdim,ydim
   integer (int_kind) :: it,jj
   integer (int_kind) :: curveSize,sb_x,sb_y,itmp,numfac
   integer (int_kind) :: subNum, sfcNum 
   logical            :: foundx

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!  retain the Cartesian distribution if nblocks_tot = nprocs
!  to avoid processors with no work
!
!----------------------------------------------------------------------
   !------------------------------------------------------
   ! Space filling curves only work if:
   ! 
   ! 	nblocks_x = 2^m1 3^n1 5^o1 where m1,n1,o1 are integers
   ! 	nblocks_y = 2^m2 3^n2 5^o2 where m2,n2,o2 are integers
   !------------------------------------------------------
   if((.not. IsFactorable(nblocks_y)) .or. (.not. IsFactorable(nblocks_x))) then 
     if (my_task == master_task) then
       write(6,*) "Factor:",nblocks_y,IsFactorable(nblocks_y),nblocks_x,IsFactorable(nblocks_x)
     endif
     ! changed from "cart" to "mic", by Junmin, Feb 1,2016.
     create_distrb_spacecurve = create_distrb_mic(nprocs, work_per_block)
     return
   endif

   !-----------------------------------------------
   ! Factor the numbers of blocks in each dimension
   !-----------------------------------------------
   xdim = Factor(nblocks_x)
   ydim = Factor(nblocks_y)
   numfac = xdim%numfact

   !---------------------------------------------
   ! Match the common factors to create SFC curve
   !---------------------------------------------
   curveSize=1
   do it=1,numfac
      call MatchFactor(xdim,ydim,itmp,foundX)
      curveSize = itmp*curveSize
   enddo
   !--------------------------------------
   ! determine the size of the sub-blocks 
   ! within the space-filling curve 
   !--------------------------------------
   sb_x = ProdFactor(xdim)
   sb_y = ProdFactor(ydim)

   call create_communicator(dist%communicator, nprocs)

   dist%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))
   dist%proc=0
   dist%local_block=0


!----------------------------------------------------------------------
!  Create the array to hold the SFC
!----------------------------------------------------------------------
   allocate(Mesh(curveSize,curveSize))
   allocate(Mesh2(nblocks_x,nblocks_y),Mesh3(nblocks_x,nblocks_y))
   Mesh  = 0
   Mesh2 = 0
   Mesh3 = 0

   allocate(idxT_i(nblocks_tot),idxT_j(nblocks_tot))


!----------------------------------------------------------------------
!  Generate the space-filling curve
!----------------------------------------------------------------------
   call GenSpaceCurve(Mesh)
   Mesh = Mesh + 1 ! make it 1-based indexing
   if(Debug) then 
     if(my_task ==0) call PrintCurve(Mesh)
   endif 
   !-----------------------------------------------
   ! Reindex the SFC to address internal sub-blocks  
   !-----------------------------------------------
   do j=1,curveSize
   do i=1,curveSize
      sfcNum = (Mesh(i,j) - 1)*(sb_x*sb_y) + 1
      do jj=1,sb_y
      do ii=1,sb_x
         subNum = (jj-1)*sb_x + (ii-1)
         i2 = (i-1)*sb_x + ii
         j2 = (j-1)*sb_y + jj
         Mesh2(i2,j2) = sfcNum + subNum
      enddo
      enddo
   enddo
   enddo
   !------------------------------------------------
   ! create a linear array of i,j coordinates of SFC
   !------------------------------------------------
   idxT_i=0;idxT_j=0
   do j=1,nblocks_y
     do i=1,nblocks_x
        n = (j-1)*nblocks_x + i
        ig = Mesh2(i,j)
        if(work_per_block(n) /= 0) then 
            idxT_i(ig)=i;idxT_j(ig)=j
        endif
     enddo
   enddo
   !-----------------------------
   ! Compress out the land blocks 
   !-----------------------------
   ii=0
   do i=1,nblocks_tot
      if(IdxT_i(i) .gt. 0) then 
         ii=ii+1
         Mesh3(idxT_i(i),idxT_j(i)) = ii
      endif
   enddo
   if(Debug) then 
     if(my_task==0) call PrintCurve(Mesh3)
   endif
   
   nblocks=ii  
   nblocksL = nblocks/nprocs
   ! every cpu gets nblocksL blocks, but the first 'extra' get nblocksL+1
   extra = mod(nblocks,nprocs)
   s1 = extra*(nblocksL+1)
   ! split curve into two curves:
   ! 1 ... s1  s2 ... nblocks
   !
   !  s1 = extra*(nblocksL+1)         (count be 0)
   !  s2 = s1+1
   !
   ! First region gets nblocksL+1 blocks per partition
   ! Second region gets nblocksL blocks per partition
   if(Debug) print *,'nprocs,extra,nblocks,nblocksL,s1: ', &
                nprocs,extra,nblocks,nblocksL,s1

   do j=1,nblocks_y
   do i=1,nblocks_x
      n = (j-1)*nblocks_x + i
!      i2 = idxT_i(n)
!      j2 = idxT_j(n)
      ii = Mesh3(i,j)
      if(ii>0) then 
        !DBG if(my_task ==0) print *,'i,j,ii:= ',i,j,ii
        if(ii<=s1) then 
           ii=ii-1
           tmp1 = ii/(nblocksL+1)
           dist%proc(n) = tmp1+1
        else
           ii=ii-s1-1
           tmp1 = ii/nblocksL
           dist%proc(n) = extra + tmp1 + 1
        endif
      endif
   enddo
   enddo

!----------------------------------------------------------------------
!  Reset the dist data structure 
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   do n=1,nblocks_tot
      pid = dist%proc(n)
      if(pid>0) then 
        proc_tmp(pid) = proc_tmp(pid) + 1
        dist%local_block(n) = proc_tmp(pid)
      endif
   enddo

   if(Debug) then 
      if(my_task==0) print *,'dist%proc:= ',dist%proc
      print *,'IAM: ',my_task,' SpaceCurve: Number of blocks {total,local} :=', &
                nblocks_tot,nblocks,proc_tmp(my_task+1)
   endif

   deallocate(proc_tmp)
   ierr=1

   deallocate(Mesh,Mesh2,Mesh3)
   deallocate(idxT_i,idxT_j)
!----------------------------------------------------------------------
   create_distrb_spacecurve = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_spacecurve


!**********************************************************************
!BOP
! !IROUTINE: create_distribution_tiny
! !INTERFACE:

 subroutine create_distribution_tiny(dist)

! !DESCRIPTION:
!  This subroutine append the tiny blocks that derived from the normal blocks
!  distributed to this process.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (distrb), intent(inout) :: &
      dist  ! resulting structure describing Cartesian

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      n               ! loop indices

   integer (int_kind),dimension(:),allocatable :: &
      local_proc,local_block  ! tmp buffer for reallocate

   allocate(local_proc(nblocks_tot), &
            local_block(nblocks_tot))

   local_proc(1:nblocks_tot) = dist%proc(1:nblocks_tot)
   local_block(1:nblocks_tot) = dist%local_block(1:nblocks_tot)
   deallocate(dist%proc, dist%local_block)
   allocate(dist%proc(nblocks_tot+nblocks_tot_tiny), &
            dist%local_block(nblocks_tot+nblocks_tot_tiny))
   dist%proc(1:nblocks_tot) = local_proc(1:nblocks_tot)
   dist%local_block(1:nblocks_tot) = local_block(1:nblocks_tot)
   deallocate(local_proc, local_block)
   
!----------------------------------------------------------------------
!  Reaplace the normal(macro) blocks with tiny blocks ? 
!  Or keep the original normal blocks info somewhere?
   do n=1,nblocks_tot_tiny
      dist%proc(n+nblocks_tot) = my_task + 1
      dist%local_block(n+nblocks_tot) = n
   end do

!----------------------------------------------------------------------
!EOC

 end subroutine create_distribution_tiny


!***********************************************************************
!BOP
! !IROUTINE: create_distrb_mic
! !INTERFACE:

 function create_distrb_mic(nprocs, work_per_block)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
!  v2 by Junmin, 2015.10.28.

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs     ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_mic  ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      k, i, j              ,&! dummy loop indices
      n                    ,&! dummy loop indices
      proc_id              ,&! dummy loop indices
      nblocks_eff          ,&! num of non-null blocks
      nblocks_tripole      ,&! num of tripole non-null blocks
      ceilval        ,&! temp ceiling for even distribution
      ceilpos          ! max position where ceiling value be taken

   integer (int_kind),dimension(nnode) :: &
      max_blknum_nodes        !
   integer (int_kind),dimension((nproc_cpu_pn+nproc_mic_pn)*nnode) :: &
      max_blknum_procs        !

   integer (int_kind),dimension(:,:),allocatable :: &
      dist_proc               ! block location on this processor

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))
   dist%proc = 0
   dist%local_block = 0

   max_blknum_nodes = 0

!----------------------------------------------------------------------
!
!  collect the histogram of block distribution in advance
!----------------------------------------------------------------------
   nblocks_eff = count(work_per_block /= 0)
   nblocks_tripole = count(work_per_block(nblocks_tot-nblocks_x+1:nblocks_tot) /= 0)
   if (my_task == master_task) then
    write(6,"(A,5I)") 'dist_cart: nblocks_tot,,nblocks_x,nblocks_y,n_eff,n_tripole:',nblocks_tot,nblocks_x,nblocks_y,nblocks_eff,nblocks_tripole
   endif

   ceilval = (nblocks_eff - 1)/nnode + 1
   ceilpos = nblocks_eff - (ceilval-1)*nnode
   max_blknum_nodes(1:ceilpos) = ceilval
   if (ceilpos < nnode) &
      max_blknum_nodes(ceilpos+1:nnode) = ceilval - 1

   call distrb_assign_by_node(work_per_block, nblocks_eff, max_blknum_nodes, dist%proc)

   call distrb_assign_by_miccenter(work_per_block, nblocks_eff, max_blknum_nodes, dist%proc)

   call distrb_assign_by_cpusurround(work_per_block, nblocks_eff, max_blknum_nodes, dist%proc)

!----------------------------------------------------------------------
!
!  set local_block info based on new distribution
!
!----------------------------------------------------------------------
   max_blknum_procs = 0
   do i=1,nblocks_x
   do j=1,nblocks_y
      n = i + (j-1)*nblocks_x
      proc_id = dist%proc(n)
      if (proc_id > 0) then
         max_blknum_procs(proc_id) = max_blknum_procs(proc_id) + 1
         dist%local_block(n) = max_blknum_procs(proc_id) 
      endif
   enddo
   enddo

   if (my_task == master_task) then
      allocate(dist_proc(nblocks_x,nblocks_y))
      dist_proc = reshape(dist%proc, (/nblocks_x,nblocks_y/))
      do j = 1,nblocks_y
         print *,'y ',j, (dist_proc(i,j), i = 1,nblocks_x)
      end do
      deallocate(dist_proc)
   endif
   ! tuning debug
   !call exit_POP(sigAbort,'dist_map completed.')

!----------------------------------------------------------------------

   create_distrb_mic = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_mic

subroutine distrb_assign_by_cpusurround(work_per_block, nblocks_eff, max_blknum_nodes, proc)

! !DESCRIPTION:
!  This subroutine attempts to allocate continous blocks in square shape to CPU
!  processes. 
!  
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ,&! amount of work per block
      max_blknum_nodes        ! specified block number per node

   integer (int_kind), intent(in) :: &
      nblocks_eff        ! non-land block number in total

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      k,i,j                ,&! dummy loop indices
      ceilval, ceilpos     ,&! max block number to each proc
      h, n                 ,&! dummy loop indices
      pid                    ! sidelen of block distrb for CPU proc per node

   integer (int_kind), dimension(nnode) :: &
      min_x,max_x        ,&! distributed range in x axis for each node
      min_y,max_y          ! distributed range in y axis for each node

   integer (int_kind),dimension(nproc_cpu_pn+nproc_mic_pn) :: &
      max_blknum_procs        !

   real (r8) :: &
      min_score, score        ! evaluate mean ratio of circumference to area 

   ! stats on distrb. range of every node
   min_x = nblocks_x+1
   max_x = 0
   min_y = nblocks_y+1
   max_y = 0
   do j=1,nblocks_y
   do i=1,nblocks_x
      n = i + (j-1) * nblocks_x
      pid = proc(n) ! must in negative value
      if (pid < 0) then 
         pid = - pid
         if (i<min_x(pid)) min_x(pid) = i
         if (i>max_x(pid)) max_x(pid) = i
         if (j<min_y(pid)) min_y(pid) = j
         if (j>max_y(pid)) max_y(pid) = j
      endif
   enddo
   enddo

   do k=1,nnode
      ! by cpu proc here, assign procid value
      ceilval = (max_blknum_nodes(k) - nblock_mic_pp*loc_nproc_mic_pn - 1)/loc_nproc_cpu_pn + 1
      ceilpos = max_blknum_nodes(k) - nblock_mic_pp*loc_nproc_mic_pn - (ceilval-1)*loc_nproc_cpu_pn
      max_blknum_procs(1:ceilpos) = ceilval
      if (ceilpos < loc_nproc_cpu_pn) &
         max_blknum_procs(ceilpos+1:loc_nproc_cpu_pn) = ceilval - 1

      ! assign blocks in the center to MIC processes
      call distrb_assign_number_surround(k,ceilval,max_blknum_procs,&
                      min_x(k),max_x(k),min_y(k),max_y(k),proc)
   enddo

!----------------------------------------------------------------------
!EOC

 end subroutine distrb_assign_by_cpusurround

subroutine distrb_assign_number_surround(k,ceilblk,max_blknum_procs,&
                           min_x, max_x, min_y, max_y, proc)

! !DESCRIPTION:
!  This subroutine attempts to give number to allocated blocks for CPU 
!  processes.
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k            ,&! node number
      ceilblk      ,&! maximum block number to allocate for each proc
      min_x,max_x  ,&!
      min_y,max_y 

   integer (int_kind), dimension(nproc_cpu_pn+nproc_mic_pn), intent(in) :: &
      max_blknum_procs     ! node number

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      mlen              ,&! sidelen of distrb shape
      nalloc,nproc      ,&! counter for allocated blocks and processes
      i,j,jj            ,&! dummy loop indices
      rgnkey            ,&! 
      bi,ei,stepi       ,&! 
      bjj,ejj,stepjj       ,&! 
      pid, n                   ! dummy loop indices

   mlen = ceiling(sqrt(real(ceilblk)))
   rgnkey = -k
   pid = (k-1) * loc_nproc_cpu_pn + 1
   nalloc = 0
   nproc = 1
   stepi = -1
   stepjj = -1
   loop_scanline: do j=min_y,max_y,mlen
     if (stepi == -1) then
       bi = min_x
       ei = max_x
       stepi = 1
     else
       bi = max_x
       ei = min_x
       stepi = -1
     endif
     do i=bi,ei,stepi
       if (stepjj == -1) then
         !.or. (j+mlen-1 >= max_y .and. max_y == nblocks_y) &! heavy tripole blocks
         bjj = j
         ejj = min(j+mlen-1, max_y)
         stepjj = 1
       else
         bjj = min(j+mlen-1, max_y)
         ejj = j
         stepjj = -1
       endif
       do jj=bjj,ejj,stepjj
         n = i + (jj-1) * nblocks_x
         if (proc(n) == rgnkey) then
           proc(n) = pid
           nalloc = nalloc + 1
           if (nalloc >= max_blknum_procs(nproc)) then
              nalloc = 0
              nproc = nproc + 1
              pid = pid + 1
              if (nproc > loc_nproc_cpu_pn) exit loop_scanline
           endif
         endif
       enddo
     enddo
   end do loop_scanline
!----------------------------------------------------------------------
!EOC

 end subroutine distrb_assign_number_surround

subroutine distrb_assign_by_miccenter(work_per_block, nblocks_eff, max_blknum_nodes, proc)

! !DESCRIPTION:
!  This subroutine attempts to allocate continous blocks in square shape to MIC
!  processes. subjected to:
!  in the center of each node range, minimizing the mean ratio of circumference 
!  to area. 
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ,&! amount of work per block
      max_blknum_nodes        ! specified block number per node

   integer (int_kind), intent(in) :: &
      nblocks_eff        ! non-land block number in total

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      k,i,j                ,&! dummy loop indices
      h, n                 ,&! dummy loop indices
      pid, nlen            ,&! sidelen of block distrb for MIC proc per node
      mlen                   ! sidelen of block distrb for each MIC proc

   integer (int_kind), dimension(nnode) :: &
      center_x,center_y  ,&! coordinate of center point for each node
      min_x,max_x        ,&! distributed range in x axis for each node
      min_y,max_y          ! distributed range in y axis for each node

   real (r8) :: &
      min_score, score        ! evaluate mean ratio of circumference to area 

   nlen = ceiling(sqrt(real(nblock_mic_pp)*loc_nproc_mic_pn))
   if (nlen == 0) then
      return
   endif

   ! stats on distrb. range of every node
   min_x = nblocks_x+1
   max_x = 0
   min_y = nblocks_y+1
   max_y = 0
   center_x = 0
   center_y = 0
   do j=1,nblocks_y
   do i=1,nblocks_x
      n = i + (j-1) * nblocks_x
      pid = -proc(n) ! must in negative value
      if (pid /= 0) then 
         if (i<min_x(pid)) min_x(pid) = i
         if (i>max_x(pid)) max_x(pid) = i
         if (j<min_y(pid)) min_y(pid) = j
         if (j>max_y(pid)) max_y(pid) = j
         center_x(pid) = center_x(pid) + i
         center_y(pid) = center_y(pid) + j
      endif
   enddo
   enddo
   do k=1,nnode
!      center_x(k) = (min_x(k) + max_x(k))/2 
!      center_y(k) = (min_y(k) + max_y(k))/2 
      ! try gravity point
      center_x(k) = center_x(k) / max_blknum_nodes(k)
      center_y(k) = center_y(k) / max_blknum_nodes(k)
      ! assign blocks in the center to MIC processes
      call distrb_assign_in_center(k,nlen,center_x(k),center_y(k),&
                      min_x(k),max_x(k),min_y(k),max_y(k),proc)
   enddo

   ! allocate blocks to MIC proc per node
   mlen = ceiling(sqrt(real(nblock_mic_pp)))
   do k=1,nnode
      call distrb_number_in_center(k,mlen,&
                      min_x(k),max_x(k),min_y(k),max_y(k),proc)
   enddo
!----------------------------------------------------------------------
!EOC

 end subroutine distrb_assign_by_miccenter

subroutine distrb_assign_in_center(k,squarelen,center_x, center_y, &
                           min_x, max_x, min_y, max_y, proc)

! !DESCRIPTION:
!  This subroutine attempts to allocate continous blocks by node in global
!  range, minimizing the mean ratio of circumference to area. The method is 
!  similar to space-filling curve.
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k            ,&! node number
      squarelen    ,&! sidelen of square area
      center_x     ,&! 
      center_y     ,&!
      min_x,max_x  ,&!
      min_y,max_y 

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      nalloc            ,&! counter for allocated blocks
      i,j               ,&! dummy loop indices
      pin_x,pin_y       ,&! 
      lefti,righti,step_j,&! 
      cur_x,cur_y,cur_n ,&! 
      assignkey         ,&! 
      h, n                   ! dummy loop indices

   logical (log_kind) :: &
      bexted              ! 

   integer (int_kind), dimension(4) :: &
      dir_x,dir_y            ,&! 4 search directions
      exts,bounds              ! extented range and bounded range

      nalloc = 0
      bexted = .true.
      dir_x = (/1,0,-1,0/)
      dir_y = (/0,1,0,-1/)
      exts = 0
      bounds = (/max_x-1,max_y-1,min_x+1,min_y+1/)
      assignkey = nnode * 2 * loc_nproc_cpu_pn + k

      ! allocate center point
      cur_n = center_x + (center_y - 1) * nblocks_x
      if ((proc(cur_n) == -k .and. &
          .not. distrb_on_boundary(proc,center_x,center_y,-k))) then
         proc(cur_n) = assignkey
         nalloc = nalloc + 1
      endif

      ! dir_x(4),dir_y(4),exts(4)=0,bounds(4)=(min+1,max-1)
      loop_rr_search: do while (nalloc < nblock_mic_pp * loc_nproc_mic_pn .and. bexted)
         bexted = .false.
         loop_four_dir: do i=1,4
            if (exts(i)<bounds(i)) then
               exts(i) = exts(i) + 1
               bexted = .true.
            else
               cycle loop_four_dir
            endif
            ! scanline
            pin_x = center_x + dir_x(i) * exts(i)
            pin_y = center_y + dir_y(i) * exts(i)
            righti = mod(i,4)+1
            lefti = mod(i+2,4)+1

            step_j = (dir_x(lefti) + dir_y(lefti))
            do j=0 * step_j, exts(lefti) * step_j, step_j
               cur_x = pin_x + dir_x(lefti) * abs(j)
               cur_y = pin_y + dir_y(lefti) * abs(j)
               cur_n = cur_x + (cur_y - 1) * nblocks_x
               if ((proc(cur_n) == -k .and. &
                   .not. distrb_on_boundary(proc,cur_x,cur_y,-k))) then
                  proc(cur_n) = assignkey
                  nalloc = nalloc + 1
                  if (nalloc >= nblock_mic_pp * loc_nproc_mic_pn) then
                     exit loop_rr_search
                  endif
               endif
            enddo
            step_j = (dir_x(righti) + dir_y(righti))
            do j=1 * step_j,exts(righti) * step_j, step_j
               cur_x = pin_x + dir_x(righti) * abs(j)
               cur_y = pin_y + dir_y(righti) * abs(j)
               cur_n = cur_x + (cur_y - 1) * nblocks_x
               if ((proc(cur_n) == -k .and. &
                   .not. distrb_on_boundary(proc,cur_x,cur_y,-k))) then
                  proc(cur_n) = assignkey
                  nalloc = nalloc + 1
                  if (nalloc >= nblock_mic_pp * loc_nproc_mic_pn) then
                     exit loop_rr_search
                  endif
               endif
            enddo
         end do loop_four_dir
      end do loop_rr_search
!----------------------------------------------------------------------
!EOC

 end subroutine distrb_assign_in_center

subroutine distrb_number_in_center(k,mlen,&
                           min_x, max_x, min_y, max_y, proc)

! !DESCRIPTION:
!  This subroutine attempts to give number to allocated blocks for MIC 
!  processes.
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k            ,&! node number
      mlen         ,&! sidelen of square area
      min_x,max_x  ,&!
      min_y,max_y 

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      nalloc,nproc      ,&! counter for allocated blocks and processes
      i,j,jj            ,&! dummy loop indices
      rgnkey            ,&! 
      bi,ei,stepi       ,&! 
      bjj,ejj,stepjj       ,&! 
      pid, n                   ! dummy loop indices

   rgnkey = nnode * 2 * loc_nproc_cpu_pn + k
   pid = nnode * loc_nproc_cpu_pn + (k-1) * loc_nproc_mic_pn + 1
   nalloc = 0
   nproc = 0
   stepi = -1
   stepjj = -1
   loop_scanline: do j=min_y,max_y,mlen
     if (stepi == -1) then
       bi = min_x
       ei = max_x
       stepi = 1
     else
       bi = max_x
       ei = min_x
       stepi = -1
     endif
     do i=bi,ei,stepi
       if (stepjj == -1) then
         bjj = j
         ejj = min(j+mlen-1, max_y)
         stepjj = 1
       else
         bjj = min(j+mlen-1, max_y)
         ejj = j
         stepjj = -1
       endif
       do jj=bjj,ejj,stepjj
         n = i + (jj-1) * nblocks_x
         if (proc(n) == rgnkey) then
           proc(n) = pid
           nalloc = nalloc + 1
           if (nalloc >= nblock_mic_pp) then
              nalloc = 0
              nproc = nproc + 1
              pid = pid + 1
              if (nproc > loc_nproc_mic_pn) exit loop_scanline
           endif
         endif
       enddo
     enddo
   end do loop_scanline
!----------------------------------------------------------------------
!EOC

 end subroutine distrb_number_in_center

subroutine distrb_assign_by_node(work_per_block, nblocks_eff, max_blknum_nodes, proc)

! !DESCRIPTION:
!  This subroutine attempts to allocate continous blocks by node in global
!  range, minimizing the mean ratio of circumference to area. The method is 
!  similar to space-filling curve.
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ,&! amount of work per block
      max_blknum_nodes        ! specified block number per node

   integer (int_kind), intent(in) :: &
      nblocks_eff        ! non-land block number in total

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      k, i                 ,&! dummy loop indices
      h, n                 ,&! dummy loop indices
      nbk_lb,nbk_ub        ,&!
      nbkwidth,min_h       ,&! 
      anchor_x,anchor_y    ,&!
      nstrip,step_x      !

   real (r8) :: &
      min_score, score        ! evaluate mean ratio of circumference to area 

   ! define allocation direction and pave width
   nbk_ub = ceiling(sqrt(real(nblocks_tot)/nnode))
   nbk_lb = ceiling(sqrt(real(nblocks_eff)/nnode))
   nbkwidth = ceiling(real(nbk_ub+nbk_lb)/2)
   min_h = 0 
   min_score = 4.0 * nblocks_tot
   do h=max(nbkwidth-1,1),nbkwidth+1,1
     ! initialize proc
     do i=1,nblocks_tot
        if (work_per_block(i) /= 0) then
           proc(i) = availf
        else
           proc(i) = 0
        endif
     enddo
     ! default allocation direction is +-x, S-shaped in incremental y axis
     anchor_x = 1
     anchor_y = 1
     nstrip = 1
     step_x = 1
     do k=1,nnode
        ! distribute with specified num and flag '-k', within blocks having 'availf' values
        ! by node here, assign negative k value
        call distrb_assign_in_group(availf, max_blknum_nodes(k), -k, h, proc,&
                        anchor_x, anchor_y, nstrip,step_x)
     enddo
     call eval_overall_circum_area_ratio(proc,score)
     if (my_task == master_task) then
     write(6,*) 'score=', score
     endif
     if (score < min_score) then
        min_score = score
        min_h = h
     endif
   enddo
   ! allocate with min_h
     ! initialize proc
     do i=1,nblocks_tot
        if (work_per_block(i) /= 0) then
           proc(i) = availf
        else
           proc(i) = 0
        endif
     enddo
     ! default allocation direction is +-x, S-shaped in incremental y axis
     anchor_x = 1
     anchor_y = 1
     nstrip = 1
     step_x = 1
     do k=1,nnode
        ! distribute with specified num and flag '-k', within blocks having 'availf' values
        ! by node here, assign negative k value
        call distrb_assign_in_group(availf, max_blknum_nodes(k), -k, min_h, proc,&
                        anchor_x, anchor_y, nstrip,step_x)
     enddo
!----------------------------------------------------------------------
!EOC

 end subroutine distrb_assign_by_node


subroutine distrb_assign_in_group(rangekey, blocknum, newkey, bkwidth,&
                proc,anchor_x,anchor_y,nstrip,step_x)

! !DESCRIPTION:
!  This subroutine attempts to allocate continous blocks within specified
!  range, minimizing the ratio of circumference to area. The method is similar
!  to space-filling curve.
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      rangekey                     ,&! available blocks having the key value
      blocknum                     ,&! num of blocks to be allocated
      newkey                       ,&! the new value to set, once allocated
      bkwidth                       ! block#/step in the allocation direction

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

   integer (int_kind), intent(inout) :: &
      anchor_x,anchor_y             ,&! start point
      nstrip               ,&! strip order in y axis
      step_x                 ! 1 or -1, 

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, jj, n                ,&! 
      last_x,last_y              ,&!
      lb_y,ub_y              ,&!
      nalloc                     ,&!
      area1, area2

   logical (log_kind) :: &
      multiline
   integer (int_kind) :: &
      pid_neighbor,expand_x,expand_y,shrink_x,shrink_y,&!
      min_order,cur_order                  ,&!
      m                                    ,&!
      ncol_line1,ncol_line2                  ! 

   real (r8) :: &
      ratio1, ratio2                ,&!
      slice_ratio, square_ratio

   i = anchor_x
   j = anchor_y
   lb_y = (nstrip-1)*bkwidth + 1
   ub_y = (nstrip)*bkwidth
   multiline = .false.
   if (ub_y + bkwidth -1 > nblocks_y) then ! extend ub_y
      ub_y = nblocks_y
   endif
   nalloc = 0
   last_x = -1
   last_y = -1
   alloc_loop_1: do while (nalloc < blocknum)
      n = i + (j-1)*nblocks_x
      if (n > nblocks_tot) then
         exit alloc_loop_1
      endif
      if (proc(n) == rangekey) then
         proc(n) = newkey
         nalloc = nalloc + 1
         last_x = i
         last_y = j
      endif
      if (nalloc < blocknum) then
      ! move forward
      if (j < ub_y) then
         j = j+1
      else
         j = lb_y
         i = i+step_x
         if (i<=0 .or. i>nblocks_x) then
            ncol_line1 = abs(i-anchor_x)
            multiline = .true.
            j = j+bkwidth
            if (i> nblocks_x) then
               i = nblocks_x
               step_x = -1
            else
               i = 1
               step_x = 1
            endif
            lb_y = lb_y+bkwidth
            ub_y = ub_y+bkwidth
            if (ub_y + bkwidth -1 > nblocks_y) then ! extend ub_y
               ub_y = nblocks_y
            endif
            nstrip = nstrip + 1
         endif
      endif
      endif ! (nalloc < blocknum) then move forward
   enddo alloc_loop_1
   ! update anchor point
   anchor_x = i
   anchor_y = j
   ! fix extreme long shape in y axis
   ! assume taking up two lines at most
   if (multiline) then
      if (step_x > 0) then
         ncol_line2 = last_x
      else
         ncol_line2 = nblocks_x - last_x + 1
      endif
      ! shorten the narrower line
      if (ncol_line1 < ncol_line2) then
         expand_x = last_x - step_x * (ncol_line2-ncol_line1)
         expand_y = lb_y - 1
         shrink_x = last_x - step_x * (ncol_line2-ncol_line1)
         shrink_y = lb_y - bkwidth
         pid_neighbor = sign(abs(newkey) - 1, newkey)
         do while (shrink_y < expand_y)
            n = expand_x + (expand_y-1)*nblocks_x
            loop_lookup_neigh: do while (proc(n) /= pid_neighbor)
               expand_x = expand_x + step_x
               if (step_x > 0) then
                  if (expand_x > last_x) then ! (step_x*expand_x > step_x*last_x)
                     expand_x = last_x - step_x *(ncol_line2-ncol_line1)
                     expand_y = expand_y - 1
                  endif
               else
                  if (expand_x < last_x) then ! (step_x*expand_x > step_x*last_x)
                     expand_x = last_x - step_x *(ncol_line2-ncol_line1)
                     expand_y = expand_y - 1
                  endif
               endif
               if (expand_y <= shrink_y) exit loop_lookup_neigh
               n = expand_x + (expand_y-1)*nblocks_x
            end do loop_lookup_neigh
            m = shrink_x + (shrink_y-1)*nblocks_x
            loop_lookup_self: do while (proc(m) /= newkey)
               shrink_x = shrink_x - step_x
               if (step_x > 0) then
                  if (shrink_x < 1) then
                     shrink_x = last_x - step_x *(ncol_line2-ncol_line1)
                     shrink_y = shrink_y + 1
                  endif
               else
                  if (shrink_x > nblocks_x) then
                     shrink_x = last_x - step_x *(ncol_line2-ncol_line1)
                     shrink_y = shrink_y + 1
                  endif
               endif
               if (shrink_y >= expand_y) exit loop_lookup_self
                m = shrink_x + (shrink_y-1)*nblocks_x
            end do loop_lookup_self
            if (proc(n) == pid_neighbor .and. proc(m) == newkey) then
               proc(n) = newkey
               proc(m) = pid_neighbor
            endif
         enddo
      else
       if (ncol_line1 > ncol_line2) then
         expand_x = last_x + step_x
         expand_y = lb_y
         shrink_x = last_x + step_x
         shrink_y = ub_y
         min_order = nblocks_tot+1
         ! expand rangekey and newkey region along -y and +y axis, respectively.
         do while (shrink_y > expand_y)
            n = expand_x + (expand_y-1)*nblocks_x
            loop_lookup_avail: do while (proc(n) /= rangekey)
               expand_x = expand_x + step_x
               if (step_x > 0) then
                  if (expand_x > ncol_line1) then
                     expand_x = last_x
                     expand_y = expand_y + 1
                  endif
               else
                  if (expand_x < nblocks_x-ncol_line1+1) then
                     expand_x = last_x
                     expand_y = expand_y + 1
                  endif
               endif
               if (expand_y >= shrink_y) exit loop_lookup_avail
               n = expand_x + (expand_y-1)*nblocks_x
            end do loop_lookup_avail
            m = shrink_x + (shrink_y-1)*nblocks_x
            loop_lookup_self2: do while (proc(m) /= newkey)
               shrink_x = shrink_x - step_x
               if (step_x > 0) then
                  if (shrink_x < 1) then
                     shrink_x = last_x
                     shrink_y = shrink_y - 1
                  endif
               else
                  if (shrink_x > nblocks_x) then
                     shrink_x = last_x
                     shrink_y = shrink_y - 1
                  endif
               endif
               if (shrink_y <= expand_y) exit loop_lookup_self2
                m = shrink_x + (shrink_y-1)*nblocks_x
            end do loop_lookup_self2
            if (proc(n) == rangekey .and. proc(m) == newkey) then
               proc(n) = newkey
               proc(m) = rangekey
               ! update anchor point
               ! find the shrink point with minimium order.
               if (step_x > 0) then
               cur_order = (shrink_x-1)*(ub_y-lb_y) + (shrink_y-lb_y+1)
               else
               cur_order = (nblocks_x-shrink_x)*(ub_y-lb_y) + (shrink_y-lb_y+1)
               endif
               if (cur_order < min_order) then
                  min_order = cur_order
                  anchor_x = shrink_x
                  anchor_y = shrink_y
               endif
            endif
         enddo ! do while, for reshaping
       endif ! Top wide/bottom narrow
      endif ! which is narrower
   endif ! if (multiline)
   call eval_circum_area_ratio(newkey, proc, ratio1, area1)
   call eval_circum_area_ratio(rangekey, proc, ratio2, area2)
   slice_ratio = ratio1*area2 + ratio2*area1
!   write(6,*) 'eval_slice:loc:', newkey, ratio1, area1
!   write(6,*) 'eval_slice:rest:', rangekey, ratio2, area2
!   write(6,*) 'eval_slice:slice:', slice_ratio

!----------------------------------------------------------------------
!EOC

 end subroutine distrb_assign_in_group

 subroutine eval_overall_circum_area_ratio(proc, ratio)

! !DESCRIPTION:
!  This subroutine evaluate the circumference to area ratio for 
!  global grid.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(in) :: &
      proc           ! number of procs in each dimension

! !OUTPUT PARAMETERS:
   real (r8), intent(out) :: &
      ratio                       ! available blocks having the key value

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind), dimension(-nnode:((nproc_cpu_pn+nproc_mic_pn)*nnode)) :: &
      counts     ,&! 
      circums

   integer (int_kind) :: &
      i, j, n , k                ,&! 
      ii, jj, nn                 ,&! 
      pid,ncount,ncircum

   integer (int_kind), dimension(4) :: &
      dirs_x, dirs_y                  !

   dirs_x = (/-1,-1, 1,1/)
   dirs_y = (/-1, 1,-1,1/)

   counts = 0
   circums = 0
   do j=1,nblocks_y
   do i=1,nblocks_x
      n = i + (j-1)*nblocks_x
      pid = proc(n)
      if (pid < -nnode .or. pid >((nproc_cpu_pn+nproc_mic_pn)*nnode) ) then
         if (my_task == master_task) then
         write(6,*) 'error with (x,y):',i,j
         endif
      else
        if (pid /= 0) then
         counts(pid) = counts(pid) + 1
         do k=1,4
           ii = i + dirs_x(k)
           jj = j + dirs_y(k)
           if (ii>=1 .and. ii<=nblocks_x .and. &
               jj>=1 .and. jj<=nblocks_y) then
             nn = ii + (jj-1)*nblocks_x
             if (proc(nn) /= 0 .and. proc(nn) /= pid) then
               circums(pid) = circums(pid) + 1
             endif
           endif
         enddo
        endif
      endif
   enddo
   enddo

   ncount = 0
   ncircum = 0
   do i=-nnode,-1,1
      if (counts(i) /= 0) then
         ncount = ncount + counts(i) 
         ncircum = ncircum + circums(i) 
      else
         if (my_task == master_task) then
         !write(6,*) 'warning: pid empty task.',i
         endif
      endif
   enddo
   do i=1,((nproc_cpu_pn+nproc_mic_pn)*nnode),1
      if (counts(i) /= 0) then
         ncount = ncount + counts(i) 
         ncircum = ncircum + circums(i) 
      else
         if (my_task == master_task) then
         !write(6,*) 'warning: pid empty task.',i
         endif
      endif
   enddo
   ratio = real(ncircum)/ncount
   !write(6,*) 'ratio:', rangekey, ratio, area

!----------------------------------------------------------------------
!EOC

 end subroutine eval_overall_circum_area_ratio


 subroutine eval_circum_area_ratio(rangekey, proc, ratio, area)

! !DESCRIPTION:
!  This subroutine evaluate the circumference to area ratio for 
!  specified region.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      rangekey                       ! available blocks having the key value

   integer (int_kind), dimension(:), intent(in) :: &
      proc           ! number of procs in each dimension

! !OUTPUT PARAMETERS:
   real (r8), intent(out) :: &
      ratio                       ! available blocks having the key value
   integer (int_kind), intent(out) :: &
      area                       ! available blocks having the key value

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n , k                ,&! 
      ii, jj, nn                 ,&! 
      ncount     ,&! 
      ncircum

   integer (int_kind), dimension(4) :: &
      dirs_x, dirs_y                  !

   dirs_x = (/-1,-1, 1,1/)
   dirs_y = (/-1, 1,-1,1/)

   ncount = 0
   ncircum = 0
   do j=1,nblocks_y
   do i=1,nblocks_x
      n = i + (j-1)*nblocks_x
      if (proc(n) == rangekey) then
         ncount = ncount + 1
         do k=1,4
           ii = i + dirs_x(k)
           jj = j + dirs_y(k)
           if (ii>=1 .and. ii<=nblocks_x .and. &
               jj>=1 .and. jj<=nblocks_y) then
             nn = ii + (jj-1)*nblocks_x
             if (proc(nn) /= 0 .and. proc(nn) /= rangekey) then
               ncircum = ncircum + 1
             endif
           endif
         enddo
      endif
   enddo
   enddo
   if (ncount /= 0) then
      ratio = ncircum * 1.0 / ncount
   else
      ratio = 0.0
   endif
   area = ncount
   !write(6,*) 'ratio:', rangekey, ratio, area

!----------------------------------------------------------------------
!EOC

 end subroutine eval_circum_area_ratio

!DIR$ ATTRIBUTES FORCEINLINE :: distrb_fill_by_slice
 subroutine distrb_fill_by_slice(proc, blocknum, rangekey, newkey,&
                    weight_x, weight_y, st_x,st_y, &
                    start_x,end_x,step_x,&
                    start_y,end_y,step_y)
                  

! !DESCRIPTION:
!  This subroutine attempts to allocate continous blocks in a slice-filling
!  way in the specified region.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      blocknum                     ,&! num of blocks to be allocated
      rangekey                     ,&! available blocks having the key value
      newkey                       ,&! assign to blocks if alllcated
      st_x,st_y                    ,&
      start_x,end_x,step_x,         &
      start_y,end_y,step_y

   real (r8), intent(in) :: &
      weight_x, weight_y           ! weight point

! !OUTPUT PARAMETERS:
   integer (int_kind), dimension(:), intent(inout) :: &
      proc           ! number of procs in each dimension

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n                  ,&! 
      nalloc                     ! 

   nalloc = 0
   if (abs(weight_y - st_y) >= abs(weight_x - st_x)) then
   ! allocate in x first
   try_loop1_x: do j=start_y,end_y,step_y
   do i=start_x,end_x,step_x
      n = i + (j-1)*nblocks_x
      if (proc(n) == rangekey) then
         proc(n) = newkey
         nalloc = nalloc + 1
      endif
      if (nalloc >= blocknum) then
         exit try_loop1_x
      endif
   enddo
   enddo try_loop1_x
   else
   ! allocate in y first
   try_loop1_y: do i=start_x,end_x,step_x
   do j=start_y,end_y,step_y
      n = i + (j-1)*nblocks_x
      if (proc(n) == rangekey) then
         proc(n) = newkey
         nalloc = nalloc + 1
      endif
      if (nalloc >= blocknum) then
         exit try_loop1_y
      endif
   enddo
   enddo try_loop1_y

   endif
!----------------------------------------------------------------------
!EOC

 end subroutine distrb_fill_by_slice

!DIR$ ATTRIBUTES FORCEINLINE :: distrb_on_boundary
 function distrb_on_boundary(proc, iblock,jblock, rangekey)

! !DESCRIPTION:
!  This fuction judge if the specified block is located on the boundary 
!  of the region, which means there exists one of its non-land neighbour 
!  belonging to any other node.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(in) :: &
      proc           ! number of procs in each dimension

   integer (int_kind), intent(in) :: &
      iblock                     ,&! num of blocks to be allocated
      jblock                     ,&! available blocks having the key value
      rangekey                     ! assign to blocks if alllcated

! !OUTPUT PARAMETERS:
   logical (log_kind) :: &
      distrb_on_boundary           ! result

! use the global nblocks_x, nblocks_y
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      k, ii, jj, nn              ,&! 
      min_proc_range,max_proc_range,&!
      min_proc_range2,max_proc_range2

   integer (int_kind), dimension(4) :: &
      dirs_x, dirs_y                  !

   dirs_x = (/1,0,-1,0/)
   dirs_y = (/0,1,0,-1/)

   if (iblock == 1 .or. iblock == nblocks_x .or. &
       jblock == 1 .or. jblock == nblocks_y) then
      distrb_on_boundary = .true.
   else
      min_proc_range = (-rangekey-1)*loc_nproc_cpu_pn+1
      max_proc_range = (-rangekey)*loc_nproc_cpu_pn
      ! intermediate value for centralized MIC blocks
      min_proc_range2 = nnode * 2 * loc_nproc_cpu_pn + (-rangekey)
      max_proc_range2 = nnode * 2 * loc_nproc_cpu_pn + (-rangekey)
      do k=1,4
        ii = iblock + dirs_x(k)
        jj = jblock + dirs_y(k)
        !if (ii>=1 .and. ii<=nblocks_x .and. &
        !    jj>=1 .and. jj<=nblocks_y) then
          nn = ii + (jj-1)*nblocks_x
          if ((proc(nn) < 0 .and. proc(nn) /= rangekey) .or.&
              (proc(nn) > 0 .and. .not. &
               ((proc(nn)>=min_proc_range .and. proc(nn)<=max_proc_range) .or.&
                (proc(nn) ==  min_proc_range2)))) then
            distrb_on_boundary = .true.
            return
          endif
        !endif
      enddo
      distrb_on_boundary = .false.
    endif
 
!----------------------------------------------------------------------
!EOC

 end function distrb_on_boundary

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_balanced
! !INTERFACE:

 function create_distrb_balanced(nprocs, work_per_block)

! !DESCRIPTION:
!  This  function distributes blocks across processors in a
!  load-balanced manner based on the amount of work per block.
!  A rake algorithm is used in which the blocks are first distributed
!  in a Cartesian distribution and then a rake is applied in each
!  Cartesian direction.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_balanced  ! resulting structure describing
                              ! load-balanced distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,n              ,&! dummy loop indices
      pid                  ,&! dummy for processor id
      local_block          ,&! local block position on processor
      max_work             ,&! max amount of work in any block
      nprocs_x             ,&! num of procs in x for global domain
      nprocs_y               ! num of procs in y for global domain

   integer (int_kind), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      work_tmp           ,&! work per row or column for rake algrthm
      proc_tmp           ,&! temp processor id for rake algrthm
      block_count          ! counter to determine local block indx

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!  retain the Cartesian distribution if nblocks_tot = nprocs
!  to avoid processors with no work
!
!----------------------------------------------------------------------

   dist = create_distrb_cart(nprocs, work_per_block)
   if (nblocks_tot == nprocs) then
      create_distrb_balanced = dist  ! return the result
      return
   endif

!----------------------------------------------------------------------
!
!  now re-distribute blocks using a rake in each direction
!
!----------------------------------------------------------------------

   max_work = maxval(work_per_block)

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

!----------------------------------------------------------------------
!
!  load-balance using a rake algorithm in the x-direction first
!
!----------------------------------------------------------------------

   allocate(priority(nblocks_tot))

   !*** set highest priority such that eastern-most blocks
   !*** and blocks with the least amount of work are
   !*** moved first

   do j=1,nblocks_y
   do i=1,nblocks_x
      n=(j-1)*nblocks_x + i
      if (work_per_block(n) > 0) then
         priority(n) = (max_work + 1)*(nblocks_x + i) - &
                       work_per_block(n)
      else
         priority(n) = 0
      endif
   end do
   end do

   allocate(work_tmp(nprocs_x), &
            proc_tmp(nprocs_x))

   do j=1,nprocs_y

      work_tmp(:) = 0
      do i=1,nprocs_x
         pid = (j-1)*nprocs_x + i
         proc_tmp(i) = pid
         do n=1,nblocks_tot
            if (dist%proc(n) == pid) then
               work_tmp(i) = work_tmp(i) + work_per_block(n)
            endif
         end do
      end do

      call rake (work_tmp, proc_tmp, work_per_block, priority, dist)

   end do

   deallocate(work_tmp, proc_tmp)

!----------------------------------------------------------------------
!
!  use a rake algorithm in the y-direction now
!
!----------------------------------------------------------------------

   !*** set highest priority for northern-most blocks

   do j=1,nblocks_y
   do i=1,nblocks_x
      n=(j-1)*nblocks_x + i
      if (work_per_block(n) > 0) then
         priority(n) = (max_work + 1)*(nblocks_y + j) - &
                       work_per_block(n)
      else
         priority(n) = 0
      endif
   end do
   end do

   allocate(work_tmp(nprocs_y), &
            proc_tmp(nprocs_y))

   do i=1,nprocs_x

      work_tmp(:) = 0
      do j=1,nprocs_y
         pid = (j-1)*nprocs_x + i
         proc_tmp(j) = pid
         do n=1,nblocks_tot
            if (dist%proc(n) == pid) then
               work_tmp(j) = work_tmp(j) + work_per_block(n)
            endif
         end do
      end do

      call rake (work_tmp, proc_tmp, work_per_block, priority, dist)

   end do

   deallocate(work_tmp, proc_tmp)
   deallocate(priority)

!----------------------------------------------------------------------
!
!  reset local_block info based on new distribution
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   do pid=1,nprocs
      local_block = 0
      do n=1,nblocks_tot
         if (dist%proc(n) == pid) then
            local_block = local_block + 1
            dist%local_block(n) = local_block
            proc_tmp(pid) = proc_tmp(pid) + 1
         endif
      end do
   end do

   if (minval(proc_tmp) < 1) then
      call exit_POP(sigAbort,'Load-balanced distribution failed')
   endif

   deallocate(proc_tmp)

!----------------------------------------------------------------------

   create_distrb_balanced = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_balanced

!**********************************************************************
!BOP
! !IROUTINE: proc_decomposition
! !INTERFACE:

 subroutine proc_decomposition(nprocs, nprocs_x, nprocs_y)

! !DESCRIPTION:
!  This subroutine attempts to find an optimal (nearly square)
!  2d processor decomposition for a given number of processors.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs                       ! total number or processors

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      nprocs_x, nprocs_y           ! number of procs in each dimension

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iguess, jguess               ! guesses for nproc_x,y

   real (r4) :: &
      square                       ! square root of nprocs

!----------------------------------------------------------------------
!
!  start with an initial guess that is closest to square decomp
!
!----------------------------------------------------------------------

   square = sqrt(real(nprocs))
   nprocs_x = 0
   nprocs_y = 0

   iguess = nint(square)

!----------------------------------------------------------------------
!
!  try various decompositions to find the best
!
!----------------------------------------------------------------------

   proc_loop: do
      jguess = nprocs/iguess

      if (iguess*jguess == nprocs) then ! valid decomp

         !***
         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         !***

         if (mod(nblocks_x,iguess) == 0 .and. &
             mod(nblocks_y,jguess) == 0) then
            nprocs_x = iguess
            nprocs_y = jguess
            exit proc_loop

         !***
         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         !***

         else if (mod(nblocks_x,jguess) == 0 .and. &
                mod(nblocks_y,iguess) == 0) then
            nprocs_x = jguess
            nprocs_y = iguess
            exit proc_loop

         !***
         !*** A valid decomposition, but keep searching for
         !***  a better one
         !***

         else
            if (nprocs_x == 0) then
               nprocs_x = iguess
               nprocs_y = jguess
            endif
            iguess = iguess - 1
            if (iguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         endif

      else ! invalid decomp - keep trying

         iguess = iguess - 1
         if (iguess == 0) then
            exit proc_loop
         else
            cycle proc_loop
         endif
      endif
   end do proc_loop

   if (nprocs_x == 0) then
      call exit_POP(sigAbort,'Unable to find 2d processor config')
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine proc_decomposition

!**********************************************************************
!BOP
! !IROUTINE: rake
! !INTERFACE:

 subroutine rake (proc_work, proc_id, block_work, priority, dist)

! !DESCRIPTION:
!  This subroutine performs a rake algorithm to distribute the work
!  along a vector of processors.  In the rake algorithm, a work
!  threshold is first set.  Then, moving from left to right, work
!  above that threshold is raked to the next processor in line.
!  The process continues until the end of the vector is reached
!  and then the threshold is reduced by one for a second rake pass.
!  In this implementation, a priority for moving blocks is defined
!  such that the rake algorithm chooses the highest priority
!  block to be moved to the next processor.  This can be used
!  for example to always choose the eastern-most block or to
!  ensure a block does not stray too far from its neighbors.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent(inout), dimension(:) :: &
      proc_work           ,&! amount of work per processor
      priority              ! priority for moving a given block

   integer (int_kind), intent(in), dimension(:) :: &
      block_work          ,&! amount of work per block
      proc_id               ! global processor number

   type (distrb), intent(inout) :: &
      dist                  ! distribution to change

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n, m, np1, &
      iproc, inext, &
      nprocs, nblocks, &
      last_priority, last_loc, &
      residual, &
      work_mean, work_max, work_diff, &
      iter, niters, itransfer, ntransfers, &
      min_priority

!----------------------------------------------------------------------
!
!  initialization
!
!----------------------------------------------------------------------

   nprocs  = size(proc_work)
   nblocks = size(block_work)

   !*** mean work per processor

   work_mean = sum(proc_work)/nprocs + 1
   work_max  = maxval(proc_work)
   residual = mod(work_mean,nprocs)

   min_priority = 1000000
   do n=1,nprocs
      iproc = proc_id(n)
      do i=1,nblocks
         if (dist%proc(i) == iproc) then
            min_priority = min(min_priority,priority(i))
         endif
      end do
   end do

!----------------------------------------------------------------------
!
!  do two sets of transfers
!
!----------------------------------------------------------------------

   transfer_loop: do

!----------------------------------------------------------------------
!
!     do rake across the processors
!
!----------------------------------------------------------------------

      ntransfers = 0
       do n=1,nprocs
          if (n < nprocs) then
             np1   = n+1
          else
             np1   = 1
          endif
         iproc = proc_id(n)
         inext = proc_id(np1)

         if (proc_work(n) > work_mean) then !*** pass work to next
            work_diff = proc_work(n) - work_mean

            rake1: do while (work_diff > 1)

               !*** attempt to find a block with the required
               !*** amount of work and with the highest priority
               !*** for moving (eg boundary blocks first)

               last_priority = 0
               last_loc = 0
               do i=1,nblocks
                  if (dist%proc(i) == iproc) then
                     if (priority(i) > last_priority ) then
                        last_priority = priority(i)
                        last_loc = i
                     endif
                  endif
               end do
               if (last_loc == 0) exit rake1 ! could not shift work

               ntransfers = ntransfers + 1
               dist%proc(last_loc) = inext
               if (np1 == 1) priority(last_loc) = min_priority
               work_diff = work_diff - block_work(last_loc)

               proc_work(n  ) = proc_work(n  )-block_work(last_loc)
               proc_work(np1) = proc_work(np1)+block_work(last_loc)
            end do rake1
         endif

      end do

!----------------------------------------------------------------------
!
!     increment work_mean by one and repeat
!
!----------------------------------------------------------------------

      work_mean = work_mean + 1
      if (ntransfers == 0 .or. work_mean > work_max) exit transfer_loop

   end do transfer_loop

!----------------------------------------------------------------------
!EOC

end subroutine rake

!***********************************************************************
! Begin of porting from POP_SpaceCurveMod
!BOP
! !IROUTINE: Cinco
! !INTERFACE:

   recursive function Cinco(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This subroutine implements a Cinco space-filling curve.
!  Cinco curves connect a Nb x Nb block of points where 
!  		
!        Nb = 5^p 
!
! !REVISION HISTORY:
!  same as module
!


! !INPUT PARAMETERS 

   integer(int_kind), intent(in) ::  &
	l, 	& ! level of the space-filling curve 
        type,   & ! type of SFC curve
	ma,     & ! Major axis [0,1]
	md,  	& ! direction of major axis [-1,1]
	ja,	& ! joiner axis [0,1]
	jd	  ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) :: &
	lma,		&! local major axis (next level)
	lmd,		&! local major direction (next level)
	lja,		&! local joiner axis (next level)
	ljd,		&! local joiner direction (next level)
	ltype,          &! type of SFC on next level 
        ll		 ! next level down 

   logical     :: debug = .FALSE.

!-----------------------------------------------------------------------
     ll = l
     if(ll .gt. 1) ltype = fact%factors(ll-1) ! Set the next type of space curve

     !--------------------------------------------------------------
     !  Position [0,0]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'Cinco: After Position [0,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,0]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [1,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [2,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [2,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = ma
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [2,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [1,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [1,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = MOD(ma+1,maxdim)
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [0,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [0,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,30) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [0,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,31) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [0,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = MOD(ma+1,maxdim)
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,32) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [1,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = ma
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,33) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [1,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,34) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [2,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,35) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [2,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,36) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [3,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = MOD(ma+1,maxdim)
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,37) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [4,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,3]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,38) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [4,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,39) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [3,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = ma
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,40) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [3,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,2]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = MOD(ma+1,maxdim)
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,41) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [4,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,42) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [4,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,43) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [3,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = ma
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,44) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [3,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,0]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = ja
     ljd       = jd

     if(ll .gt. 1) then
        if(debug) write(*,45) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'After Position [4,0] ',pos
     endif

 21   format('Call Cinco Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 22   format('Call Cinco Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 23   format('Call Cinco Pos [2,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 24   format('Call Cinco Pos [2,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 25   format('Call Cinco Pos [2,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 26   format('Call Cinco Pos [1,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 27   format('Call Cinco Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 28   format('Call Cinco Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 29   format('Call Cinco Pos [0,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 30   format('Call Cinco Pos [0,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 31   format('Call Cinco Pos [0,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 32   format('Call Cinco Pos [1,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 33   format('Call Cinco Pos [1,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 34   format('Call Cinco Pos [2,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 35   format('Call Cinco Pos [2,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 36   format('Call Cinco Pos [3,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 37   format('Call Cinco Pos [4,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 38   format('Call Cinco Pos [4,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 39   format('Call Cinco Pos [3,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 40   format('Call Cinco Pos [3,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 41   format('Call Cinco Pos [4,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 42   format('Call Cinco Pos [4,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 43   format('Call Cinco Pos [3,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 44   format('Call Cinco Pos [3,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 45   format('Call Cinco Pos [4,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

!EOC
!-----------------------------------------------------------------------

   end function Cinco

!***********************************************************************
!BOP
! !IROUTINE: PeanoM
! !INTERFACE:

   recursive function PeanoM(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This function implements a meandering Peano 
!  space-filling curve. A meandering Peano curve 
!  connects a Nb x Nb block of points where
!
!        Nb = 3^p
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS

   integer(int_kind), intent(in) ::  &
        l,      & ! level of the space-filling curve
        type,   & ! type of SFC curve
        ma,     & ! Major axis [0,1]
        md,     & ! direction of major axis [-1,1]
        ja,     & ! joiner axis [0,1]
        jd        ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


   integer(int_kind) :: &
        lma,            &! local major axis (next level)
        lmd,            &! local major direction (next level)
        lja,            &! local joiner axis (next level)
        ljd,            &! local joiner direction (next level)
        ltype,          &! type of SFC on next level
        ll               ! next level down

   logical     :: debug = .FALSE.

!-----------------------------------------------------------------------

     ll = l
     if(ll .gt. 1) ltype = fact%factors(ll-1) ! Set the next type of space curve
     !--------------------------------------------------------------
     !  Position [0,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [0,0] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [0,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [0,1] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [0,2]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [0,2] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [1,2]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [1,2] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [2,2]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = md
     lja        = MOD(lma+1,maxdim)
     ljd        = -lmd

     if(ll .gt. 1) then
        if(debug) write(*,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [2,2] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [2,1]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = -md
     lja        = lma
     ljd        = lmd

     if(ll .gt. 1) then
        if(debug) write(*,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [2,1] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [1,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja        = lma
     ljd        = lmd

     if(ll .gt. 1) then
        if(debug) write(*,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [1,1] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [1,0]
     !--------------------------------------------------------------
     lma        = MOD(ma+1,maxdim)
     lmd        = -md
     lja        = MOD(lma+1,maxdim)
     ljd        = -lmd

     if(ll .gt. 1) then
        if(debug) write(*,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [1,0] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [2,0]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = md
     lja        = ja
     ljd        = jd

     if(ll .gt. 1) then
        if(debug) write(*,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'PeanoM: After Position [2,0] ',pos
     endif

 21   format('Call PeanoM Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 22   format('Call PeanoM Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 23   format('Call PeanoM Pos [0,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 24   format('Call PeanoM Pos [1,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 25   format('Call PeanoM Pos [2,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 26   format('Call PeanoM Pos [2,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 27   format('Call PeanoM Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 28   format('Call PeanoM Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 29   format('Call PeanoM Pos [2,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

!EOC
!-----------------------------------------------------------------------

   end function PeanoM

!***********************************************************************
!BOP
! !IROUTINE: Hilbert
! !INTERFACE:

   recursive function Hilbert(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This function implements a Hilbert space-filling curve.
!  A Hilbert curve connect a Nb x Nb block of points where
!
!        Nb = 2^p
!
! !REVISION HISTORY:
!  same as module
!


! !INPUT PARAMETERS

   integer(int_kind), intent(in) ::  &
        l,      & ! level of the space-filling curve
        type,   & ! type of SFC curve
        ma,     & ! Major axis [0,1]
        md,     & ! direction of major axis [-1,1]
        ja,     & ! joiner axis [0,1]
        jd        ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


   integer(int_kind) :: &
        lma,            &! local major axis (next level)
        lmd,            &! local major direction (next level)
        lja,            &! local joiner axis (next level)
        ljd,            &! local joiner direction (next level)
        ltype,          &! type of SFC on next level
        ll               ! next level down

   logical     :: debug = .FALSE.

!-----------------------------------------------------------------------
     ll = l
     if(ll .gt. 1) ltype = fact%factors(ll-1) ! Set the next type of space curve
     !--------------------------------------------------------------
     !  Position [0,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'Hilbert: After Position [0,0] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [0,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'Hilbert: After Position [0,1] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [1,1]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = md
     lja        = MOD(ma+1,maxdim)
     ljd        = -md

     if(ll .gt. 1) then
        if(debug) write(*,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'Hilbert: After Position [1,1] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [1,0]
     !--------------------------------------------------------------
     lma        = MOD(ma+1,maxdim)
     lmd        = -md
     lja        = ja
     ljd        = jd

     if(ll .gt. 1) then
        if(debug) write(*,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) write(*,*) 'Hilbert: After Position [1,0] ',pos
     endif

 21   format('Call Hilbert Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 22   format('Call Hilbert Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 23   format('Call Hilbert Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 24   format('Call Hilbert Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

!EOC
!-----------------------------------------------------------------------

   end function hilbert

!***********************************************************************
!BOP
! !IROUTINE: IncrementCurve
! !INTERFACE:

   function IncrementCurve(ja,jd) result(ierr)

! !DESCRIPTION:
!   This function creates the curve which is store in the 
!   the ordered array.  The curve is implemented by 
!   incrementing the curve in the direction [jd] of axis [ja].
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS:
     integer(int_kind)  :: &
	ja, 	&! axis to increment
	jd	 ! direction along axis

! !OUTPUT PARAMETERS:
     integer(int_kind) :: ierr ! error return code

     !-----------------------------
     ! mark the newly visited point
     !-----------------------------
     ordered(pos(0)+1,pos(1)+1) = vcnt
	
     !------------------------------------
     ! increment curve and update position
     !------------------------------------
     vcnt  = vcnt + 1
     pos(ja) = pos(ja) + jd

     ierr = 0
!EOC
!-----------------------------------------------------------------------

   end function IncrementCurve

!***********************************************************************
!BOP
! !IROUTINE: log2
! !INTERFACE:

   function log2( n)

! !DESCRIPTION:
!  This function calculates the log2 of its integer 
!  input.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: n  ! integer value to find the log2
   
! !OUTPUT PARAMETERS: 

   integer(int_kind) :: log2

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) ::  tmp

   !-------------------------------
   !  Find the log2 of input value
   !-------------------------------
   if (n .le. 1) then
      log2 = 0
   else
      log2 = 1
      tmp =n
      do while (tmp/2 .ne. 1) 
         tmp=tmp/2
         log2=log2+1
      enddo 
   endif

!EOP
!-----------------------------------------------------------------------

   end function log2

!***********************************************************************
!BOP
! !IROUTINE: IsLoadBalanced
! !INTERFACE:

   function  IsLoadBalanced(nelem,npart)
   
! !DESCRIPTION:
!  This function determines if we can create 
!  a perfectly load-balanced partitioning.
!
! !REVISION HISTORY:
!  same as module

! !INTPUT PARAMETERS:

   integer(int_kind), intent(in) ::  &
	nelem,		&  ! number of blocks/elements to partition
	npart              ! size of partition

! !OUTPUT PARAMETERS:
   logical        :: IsLoadBalanced   ! .TRUE. if a perfectly load balanced 
				      ! partition is possible
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
	
   integer(int_kind)   :: tmp1 ! temporary int

!-----------------------------------------------------------------------
   tmp1 = nelem/npart

   if(npart*tmp1 == nelem ) then 
	IsLoadBalanced=.TRUE.
   else
        IsLoadBalanced=.FALSE.
   endif

!EOP
!-----------------------------------------------------------------------

   end function IsLoadBalanced

!***********************************************************************
!BOP
! !IROUTINE: GenCurve
! !INTERFACE:

   function GenCurve(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This subroutine generates the next level down
!  space-filling curve
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS

   integer(int_kind), intent(in) ::  &
        l,      & ! level of the space-filling curve
        type,   & ! type of SFC curve
        ma,     & ! Major axis [0,1]
        md,     & ! direction of major axis [-1,1]
        ja,     & ! joiner axis [0,1]
        jd        ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------

   !-------------------------------------------------
   ! create the space-filling curve on the next level  
   !-------------------------------------------------
   if(type == 2) then
      ierr = Hilbert(l,type,ma,md,ja,jd)
   elseif ( type == 3) then
      ierr = PeanoM(l,type,ma,md,ja,jd)
   elseif ( type == 5) then 
      ierr = Cinco(l,type,ma,md,ja,jd)
   endif

!EOP
!-----------------------------------------------------------------------

   end function GenCurve

    function FirstFactor(fac) result(res)
       type (factor_t) :: fac
       integer :: res
       logical :: found
       integer (int_kind) :: i

       found = .false.
       i=1
       do while (i<=fac%numfact .and. (.not. found))
          if(fac%used(i) == 0) then
                res = fac%factors(i)
                found = .true.
          endif
          i=i+1
        enddo

    end function FirstFactor

    function FindandMark(fac,val,f2) result(found)
       type (factor_t) :: fac
       integer :: val
       logical :: found
       logical :: f2
       integer (int_kind) :: i 

       found = .false.
       i=1
       do while (i<=fac%numfact .and. (.not. found))
          if(fac%used(i) == 0) then
                if(fac%factors(i) .eq. val) then
                   if(f2)  then
                      fac%used(i) = 1
                      found = .true.
                   else if( .not. f2) then
                      fac%used(i) = -1
                      found = .false.
                   endif
                endif
          endif
          i=i+1
        enddo

    end function FindandMark


   subroutine MatchFactor(fac1,fac2,val,found)
      type (factor_t) :: fac1
      type (factor_t) :: fac2
      integer :: val
      integer :: val1
      logical :: found
      logical :: tmp

      found = .false.

      val1 = FirstFactor(fac1)
!pw      print *,'Matchfactor: found value: ',val
      found = FindandMark(fac2,val1,.true.)
      tmp = FindandMark(fac1,val1,found)
      if (found) then 
        val = val1
      else
        val = 1
      endif

   end subroutine MatchFactor

   function ProdFactor(fac) result(res)

   type (factor_t) :: fac
   integer :: res
   integer (int_kind) :: i 

     res = 1
     do i=1,fac%numfact
        if(fac%used(i) <= 0) then
          res = res * fac%factors(i)
        endif
     enddo

   end function ProdFactor


   subroutine PrintFactor(msg,fac)

   
      character(len=*) :: msg
      type (factor_t) :: fac
      integer (int_kind) :: i 

      write(*,*) ' '
      write(*,*) 'PrintFactor: ',msg
      write(*,*) (fac%factors(i),i=1,fac%numfact)
      write(*,*) (fac%used(i),i=1,fac%numfact)


   end subroutine PrintFactor

!***********************************************************************
!BOP
! !IROUTINE: Factor
! !INTERFACE:

   function Factor(num) result(res)

! !DESCRIPTION:
!  This function factors the input value num into a 
!  product of 2,3, and 5.
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS:

   integer(int_kind), intent(in)  :: num  ! number to factor

! !OUTPUT PARAMETERS:

   type (factor_t)     :: res

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind)   ::  &
	tmp,tmp2,tmp3,tmp5   ! tempories for the factorization algorithm
   integer(int_kind)   :: i,n    ! loop tempories
   logical             :: found  ! logical temporary

   ! --------------------------------------
   ! Allocate allocate for max # of factors
   ! --------------------------------------
   tmp = num
   tmp2 = log2(num)
   allocate(res%factors(tmp2))
   allocate(res%used(tmp2))

   res%used = 0
   n=0

   !-----------------------
   !  Look for factors of 5
   !-----------------------
   found=.TRUE.
   do while (found)
      found = .FALSE.
      tmp5 = tmp/5
      if( tmp5*5 == tmp ) then
        n = n + 1
        res%factors(n) = 5
        found = .TRUE.
        tmp = tmp5
      endif
   enddo

   !-----------------------
   !  Look for factors of 3
   !-----------------------
   found=.TRUE.
   do while (found)
      found = .FALSE.
      tmp3 = tmp/3
      if( tmp3*3 == tmp ) then
        n = n + 1
        res%factors(n) = 3
        found = .TRUE.
        tmp = tmp3
      endif
   enddo

   !-----------------------
   !  Look for factors of 2
   !-----------------------
   found=.TRUE.
   do while (found)
      found = .FALSE.
      tmp2 = tmp/2
      if( tmp2*2 == tmp ) then
        n = n + 1
        res%factors(n) = 2
        found = .TRUE.
        tmp = tmp2
      endif
   enddo

   !------------------------------------
   ! make sure that the input value 
   ! only contains factors of 2,3,and 5  
   !------------------------------------
   tmp=1
   do i=1,n
     tmp = tmp * res%factors(i)
   enddo
   if(tmp == num) then
     res%numfact = n
   else
     res%numfact = -1
   endif

!EOP
!---------------------------------------------------------
   end function Factor

!***********************************************************************
!BOP
! !IROUTINE: IsFactorable
! !INTERFACE:

   function IsFactorable(n)
   
! !DESCRIPTION:
!  This function determines if we can factor
!   n into 2,3,and 5.  
!
! !REVISION HISTORY:
!  same as module


! !INTPUT PARAMETERS:

   integer(int_kind), intent(in)  :: n  ! number to factor

! !OUTPUT PARAMETERS:
   logical  :: IsFactorable  ! .TRUE. if it is factorable

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (factor_t)     :: fact  ! data structure to store factor information

   fact = Factor(n)
   if(fact%numfact .ne. -1) then
     IsFactorable = .TRUE.
   else
     IsFactorable = .FALSE.
   endif

!EOP
!-----------------------------------------------------------------------

   end function IsFactorable

!***********************************************************************
!BOP
! !IROUTINE: map
! !INTERFACE:

   subroutine map(l)

! !DESCRIPTION:
!   Interface routine between internal subroutines and public 
!   subroutines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer(int_kind)  :: l   ! level of space-filling curve


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind)  :: &
	d, 		 & ! dimension of curve only 2D is supported
	type,		 & ! type of space-filling curve to start off
        ierr   		   ! error return code

   d = SIZE(pos)

   pos=0
   maxdim=d
   vcnt=0

   type = fact%factors(l)
   ierr = GenCurve(l,type,0,1,0,1)


!EOP
!-----------------------------------------------------------------------

   end subroutine map

!***********************************************************************
!BOP
! !IROUTINE: PrintCurve
! !INTERFACE:

   subroutine PrintCurve(Mesh)


! !DESCRIPTION:
!  This subroutine prints the several low order 
!  space-filling curves in an easy to read format
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

     integer(int_kind), intent(in), target ::  Mesh(:,:) ! SFC to be printed

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
     integer(int_kind) ::  &
        gridsize,	   &! order of space-filling curve
        i		    ! loop temporary

!-----------------------------------------------------------------------

     gridsize = SIZE(Mesh,dim=1)

     if(gridsize == 2) then
        write (*,*) "A Level 1 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,2) Mesh(1,i),Mesh(2,i)
        enddo
     else if(gridsize == 3) then
        write (*,*) "A Level 1 Peano Meandering Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,3) Mesh(1,i),Mesh(2,i),Mesh(3,i)
        enddo
     else if(gridsize == 4) then
        write (*,*) "A Level 2 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,4) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i)
        enddo
     else if(gridsize == 5) then
        write (*,*) "A Level 1 Cinco Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,5) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i),Mesh(5,i)
        enddo
     else if(gridsize == 6) then
        write (*,*) "A Level 1 Hilbert and Level 1 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,6) Mesh(1,i),Mesh(2,i),Mesh(3,i), &
	    	      Mesh(4,i),Mesh(5,i),Mesh(6,i)
        enddo
     else if(gridsize == 8) then
        write (*,*) "A Level 3 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,8) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                      Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i)
         enddo
     else if(gridsize == 9) then
        write (*,*) "A Level 2 Peano Meandering Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,9) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                      Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                      Mesh(9,i)
         enddo
     else if(gridsize == 10) then
        write (*,*) "A Level 1 Hilbert and Level 1 Cinco Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,10) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                      Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                      Mesh(9,i),Mesh(10,i)
         enddo
     else if(gridsize == 12) then
        write (*,*) "A Level 2 Hilbert and Level 1 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,12) Mesh(1,i),Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i),Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i)
        enddo
     else if(gridsize == 15) then
        write (*,*) "A Level 1 Peano and Level 1 Cinco Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,15) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                       Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                       Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i)
        enddo
     else if(gridsize == 16) then
        write (*,*) "A Level 4 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,16) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                       Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                       Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i)
        enddo
     else if(gridsize == 18) then
        write (*,*) "A Level 1 Hilbert and Level 2 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,18) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i)
        enddo
     else if(gridsize == 20) then
        write (*,*) "A Level 2 Hilbert and Level 1 Cinco Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,20) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i)
        enddo
     else if(gridsize == 24) then
        write (*,*) "A Level 3 Hilbert and Level 1 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,24) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i)
        enddo
     else if(gridsize == 25) then
        write (*,*) "A Level 2 Cinco Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,25) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
		       Mesh(25,i)
        enddo
     else if(gridsize == 27) then
        write (*,*) "A Level 3 Peano Meandering Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,27) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                       Mesh(25,i),Mesh(26,i),Mesh(27,i)
        enddo
     else if(gridsize == 32) then
        write (*,*) "A Level 5 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,32) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i),  &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i),  &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                       Mesh(25,i),Mesh(26,i),Mesh(27,i),Mesh(28,i), &
                       Mesh(29,i),Mesh(30,i),Mesh(31,i),Mesh(32,i)
        enddo
     endif
 2 format('|',2(i2,'|'))
 3 format('|',3(i2,'|'))
 4 format('|',4(i2,'|'))
 5 format('|',5(i2,'|'))
 6 format('|',6(i2,'|'))
 8 format('|',8(i2,'|'))
 9 format('|',9(i2,'|'))
10 format('|',10(i2,'|'))
12 format('|',12(i3,'|'))
15 format('|',15(i3,'|'))
16 format('|',16(i3,'|'))
18 format('|',18(i3,'|'))
20 format('|',20(i3,'|'))
24 format('|',24(i3,'|'))
25 format('|',25(i3,'|'))
27 format('|',27(i3,'|'))
32 format('|',32(i4,'|'))

!EOC
!-----------------------------------------------------------------------

   end subroutine PrintCurve

!***********************************************************************
!BOP
! !IROUTINE: GenSpaceCurve
! !INTERFACE:

  subroutine  GenSpaceCurve(Mesh)

! !DESCRIPTION:
!  This subroutine is the public interface into the 
!  space-filling curve functionality
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT/OUTPUT PARAMETERS:
   integer(int_kind), target,intent(inout) :: &
	Mesh(:,:)		! The SFC ordering in 2D array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) ::  &
	level,   &! Level of space-filling curve		
	dim       ! dimension of SFC... currently limited to 2D

   integer(int_kind) :: gridsize   ! number of points on a side
   
!-----------------------------------------------------------------------

   !-----------------------------------------
   !  Setup the size of the grid to traverse
   !-----------------------------------------
   dim = 2
   gridsize = SIZE(Mesh,dim=1)
   fact     = factor(gridsize)
   level    = fact%numfact

   if(verbose) write(*,*) 'GenSpacecurve: level is ',level
   allocate(ordered(gridsize,gridsize))

   !--------------------------------------------
   ! Setup the working arrays for the traversal
   !--------------------------------------------
   allocate(pos(0:dim-1))
   
   !-----------------------------------------------------
   !  The array ordered will contain the visitation order
   !-----------------------------------------------------
   ordered(:,:) = 0

   call map(level) 

   Mesh(:,:) = ordered(:,:)

   deallocate(pos,ordered)

!EOP
!-----------------------------------------------------------------------

  end subroutine GenSpaceCurve 
! End of porting from POP_SpaceCurveMod
!***********************************************************************

end module distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
