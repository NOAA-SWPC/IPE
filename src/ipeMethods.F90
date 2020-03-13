module ipeMethods

  use ESMF
  use NUOPC
  use IPE_Wrapper
  use IPE_Constants_Dictionary

  implicit none

  ! -- debug
  real(ESMF_KIND_R8), parameter :: BAD_VALUE = -999._ESMF_KIND_R8

  type IPE_Model_Type
    integer,          pointer :: nodeToIndexMap(:,:)
    type (IPE_Model), pointer :: ipe
  end type

  type IPE_InternalState_Type
    type (IPE_Model_type), pointer :: model
  end type

  private

! public :: iHemiStart, iHemiEnd
! public :: numLineNodes, numLocalNodes
! public :: jmin, jSouth
  public :: IPE_InternalState_Type
  public :: IPE_Model_Type
  public :: prec

  public :: IPEFieldDiagnostics
  public :: IPEMeshCreate

contains

  !---------------------------------------------------------------------------
  ! SUBROUTINE: IPEMeshCreate
  !
  ! DESCRIPTION:
  !> \brief
  !> Build a 3D mesh from magnetic QD coordinates distributed in rectangular
  !! tiles in the (magnetic longitude, magnetic latitude) index space (mp,lp)
  !
  !> \details
  !> This subroutine selects all available heights between
  !! \c mesh_height_min and \c mesh_height_max, read as km from the
  !! \c &ipecap namelist in IPE.inp, then internally converted to m.
  !! The mesh's nodes and elements are created and indexed sequentially from
  !! the Northern to the Southern hemisphere, with nodes on the magnetic
  !! equator belonging only to the Northern hemisphere.
  !! The hemispheres are connected as hexahedral cells are built along the
  !! magnetic equator.
  !! If optional 'fill' argument is set to .true., the empty polar and
  !! equatorial regions are filled by adding geometric polar nodes at each
  !! included height and further cells to occupy the space where no magnetic
  !! field lines are available.
  !! Information about the horizontal domain decomposition is retrieved
  !! from the IPE module:
  !! \code
  !!   module_input_parameters
  !! \endcode
  !! as bounds of the local tile in magnetic index space (lps:lpe,mps:mpe).
  !! For debugging purposes, each hemisphere can be built independently by
  !! setting the module \c iHemiStart and \c iHemiEnd parameters to
  !! 0 (Northern hemisphere) or 1 (Southern hemisphere).
  !---------------------------------------------------------------------------

  subroutine IPEMeshCreate(gcomp, mesh, fill, rc)

    ! --- input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_Mesh)                :: mesh
    logical, optional, intent(in)  :: fill
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: verbosity, diagnostic
    integer :: localPet, petCount, pet, poleOwnerPet
    integer :: lps, lpe, mps, mpe
    integer :: lp, lpp, lpu, mp, mpp, mpu, kp, kpp, kpe
    integer :: lpStart, lpOffset, lHalo, mpStart, mpOffset, mHalo
    integer :: kpStart, kpEnd, kpStep, kpOffset
    integer :: numNodes, numHemiNodes, numLocalNodes
    integer :: numElems, numLongElems
    integer :: nCount
    integer :: id, in, is, im, i, j
    integer :: iHemi, jHemi
    logical :: fillFlag, fillEquator, fillPoles

    integer, dimension(1)                :: jtop
    integer, dimension(:),   allocatable :: jmin, jSouth
    integer, dimension(:),   allocatable :: numLongNodes
    integer, dimension(:),   allocatable :: numLineElems
    integer, dimension(:),   allocatable :: nodeIds, nodeOwners
    integer, dimension(:),   allocatable :: elemIds, elemTypes, elemConn
    integer, dimension(:,:), allocatable :: idNodeOffset
    integer, dimension(:,:), allocatable :: numLineNodes
    integer, dimension(:,:), allocatable :: petMap

    integer(ESMF_KIND_I4), dimension(:),       allocatable :: localBounds, globalBounds
    real   (ESMF_KIND_R8), dimension(:),       allocatable :: nodeCoords
    real   (ESMF_KIND_R8), dimension(:),       allocatable :: sendData, recvData
    real   (ESMF_KIND_R8), dimension(:,:),     allocatable :: poleCoords
    real   (ESMF_KIND_R4), dimension(:,:,:,:), allocatable :: local_grid_3d

    character(len=ESMF_MAXSTR) :: name

    type(ESMF_VM) :: vm

    type(IPE_InternalState_Type)  :: istate
    type(IPE_Model_Type), pointer :: this
    type(IPE_Model),      pointer :: ipe

    ! -- local parameters
    integer,               parameter :: iHemiStart  = 0  ! Use Northern (0) and
    integer,               parameter :: iHemiEnd    = 1  ! Southern (1) hemisphere (for debug only)
    integer, dimension(8), parameter :: iConn       = (/ 0,0,0,1,1,1,1,0 /)  ! connection array for hexahedron face (i,i+1) -> (i+2,i+3)...
    real(ESMF_KIND_R8),    parameter :: rad2deg     = 57.29577951308232087721_ESMF_KIND_R8

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- get verbosity and diagnostic levels
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! -- get internal state
    nullify(this, ipe)
    call ESMF_GridCompGetInternalState(gcomp, istate, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    this => istate % model
    ipe  => istate % model % ipe

    ! -- set local (lp,mp) bounds
    lps = ipe % mpi_layer % lp_low
    lpe = ipe % mpi_layer % lp_high
    mps = ipe % mpi_layer % mp_low
    mpe = ipe % mpi_layer % mp_high

    ! -- set local fill flags for equatorial and polar regions if this PET
    ! -- builds adjacent mesh portions
    fillFlag = .false.
    if (present(fill)) fillFlag = fill
    fillEquator = fillFlag .and. (lpe == ipe%grid%NLP)
    fillPoles   = fillFlag .and. (lps ==   1)

    ! -- check that input parameters are acceptable
    ! -- verify that min and max mesh height are provided in the right order
    if (ipe % parameters % mesh_height_min >= ipe % parameters % mesh_height_max) then
      call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
      msg="Minimim height must be smaller than maximum height", &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)
      return
    end if

    ! -- identify neighboring PETs
    ! -- get VM for current gridded component
    call ESMF_GridCompGet(gcomp, vm=vm, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- retrieve PET information
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- gather information about domain decomposition
    allocate(localBounds(4), globalBounds(4*petCount), stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    localBounds = (/ mps, mpe, lps, lpe /)
    globalBounds = 0
    call ESMF_VMAllGather(vm, localBounds, globalBounds, 4, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    deallocate(localBounds, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- debug
    if (btest(verbosity,1)) then
      if (localPet == 0) then
        write(6,'("====== IPEMeshCreate: tile bounds ======"/&
                 &"localPET     mps     mpe     lps     lpe"/&
                 &40("-"))')
        do i = 1, size(globalBounds), 4
          write(6,'(5i8)') i/4, globalBounds(i:i+3)
        end do
        write(6,'(40("="))')
        if (fillFlag) write(6,'("Filling equatorial and polar gaps ...")')
      end if
    end if

    ! -- get halo size for local DE
    call IPEGetHalo(lUBound=lpu, lHaloSize=lHalo, mUBound=mpu, mHaloSize=mHalo, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- retrieve 2D coordinates (including halo regions)
    allocate(local_grid_3d(ipe % grid % nFluxTube, lps:lpu, mps:mpu, 2), &
      stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return
    call IPEGetGridCoord(vm, local_grid_3d, globalBounds, &
      ipe % grid % nFluxTube, lps, lpu, lHalo, mps, mpu, mHalo, &
      diagnostic=diagnostic, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- determine PETs that own neighboring DEs
    allocate(petMap(0:lpu-lpe,0:mpu-mpe), stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- initialize PET map array
    petMap = -1

    do mpp = mpe, mpu
      mp = mod(mpp-1, ipe % grid % NMP) + 1    ! returns periodic longitude index
      do lp = lpe, lpu
        petMap(lp-lpe,mpp-mpe) = IPEGetOwnerPet(lp, mp, globalBounds)
      end do
    end do

    ! -- check if PET map is complete
    if (any(petMap < 0)) then
      call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
      msg="Failure to determine tile to PET mapping", &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)
      return
    end if

    ! -- select grid points between min and max height
    allocate(jmin(ipe % grid % NLP), &
      jSouth(ipe % grid % NLP), &
      numLineNodes(ipe % grid % NLP,0:1), &
      numLongNodes(0:1), &
      stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    jmin = 0
    jSouth = 0
    numLineNodes = 0
    numLongNodes = 0
    do lp = 1, ipe % grid % NLP
      in = 1
      is = ipe % grid % flux_tube_max(lp)
      im = (in + is)/2    ! highest point (apex) along field line
      if (ipe % grid % altitude(im,lp) >= ipe % parameters % mesh_height_min) then
        jtop = minloc(ipe % grid % altitude(in:im,lp), ipe % grid % altitude(in:im,lp) >= ipe % parameters % mesh_height_min)
        jmin(lp) = jtop(1) + in - 1
        is = is - jtop(1) + 1
        in = jmin(lp)
        if (ipe % grid % altitude(im,lp) > ipe % parameters % mesh_height_max) then
          jtop = maxloc(ipe % grid % altitude(in:im,lp), ipe % grid % altitude(in:im,lp) <= ipe % parameters % mesh_height_max)
          jSouth(lp) = is - jtop(1) + 1
          numLineNodes(lp,0) = jtop(1)
        else
          jSouth(lp) = im + 1
          numLineNodes(lp,0) = im - in + 1
        end if
        numLineNodes(lp,1) = is - jSouth(lp) + 1
      end if
      numLongNodes = numLongNodes + numLineNodes(lp,:)
    end do

    numNodes      = (mpu-mps+1) * sum(numLineNodes(lps:lpu,iHemiStart:iHemiEnd))
    numLocalNodes = (mpe-mps+1) * sum(numLineNodes(lps:lpe,iHemiStart:iHemiEnd))

    ! -- include supplemental polar nodes, if filling
    if (fillPoles) then
      nCount = sum(numLineNodes(lps,iHemiStart:iHemiEnd))
      numNodes      = numNodes      + nCount
      ! -- polar nodes are assigned to PET holding (lp=1, mp=1) point
      poleOwnerPet = IPEGetOwnerPet(1, 1, globalBounds)
      if (localPet == poleOwnerPet) numLocalNodes = numLocalNodes + nCount
    end if

    ! -- free up memory
    deallocate(globalBounds, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- build mesh nodes
    allocate(nodeIds(numNodes), &
      nodeOwners(numNodes), &
      nodeCoords(3*numNodes), &
      stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    nodeIds = 0
    nodeOwners = 0
    nodeCoords = 0._ESMF_KIND_R8

    i = 0
    nCount = 0
    do iHemi = iHemiStart, iHemiEnd
      mpStart = (iHemi-iHemiStart)*ipe % grid % NMP*numLongNodes(0)
      lpStart = sum(numLineNodes(:lps-1,iHemi))
      do mpp = mps, mpu
        mp = mod(mpp-1, ipe % grid % NMP) + 1
        j = max(mpp-mpe, 0)
        mpOffset = (mp-1) * numLongNodes(iHemi) + mpStart
        lpOffset = lpStart
        do lp = lps, lpu
          kpStep   = numLineNodes(lp,iHemi) - 1
          kpStart  = iHemi * kpStep + 1
          kpEnd    = (1 - iHemi) * kpStep + 1
          kpStep   = 1 - 2 * iHemi
          kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
          pet = petMap(max(lp-lpe,0),j)
          id = 0
          do kpp = kpStart, kpEnd, kpStep
            id = id + 1
            nCount = nCount + 1
            nodeIds(nCount) = id + lpOffset + mpOffset
            nodeOwners(nCount) = pet
            kp = kpp + kpOffset
            nodeCoords(i + 1) = rad2deg * local_grid_3d(kp,lp,mpp,1)
            nodeCoords(i + 2) = 90._ESMF_KIND_R8 - rad2deg * local_grid_3d(kp,lp,mpp,2)
            nodeCoords(i + 3) = 1.e-03_ESMF_KIND_R8 * ipe % grid % altitude(kp,lp)
            i = i + 3
          end do
          lpOffset = lpOffset + numLineNodes(lp,iHemi)
        end do
      end do
    end do

    ! -- compute geometric pole coordinates, if needed
    if (fillFlag) then
      ! -- retrieve pole coordinates on all PETs for simplicity
      allocate(poleCoords(ipe%grid%nFluxTube,2), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
      call IPEGetPoleCoord(vm, poleCoords, ipe%grid%nFluxTube, lps, mps, mpe, ipe%grid%nmp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

    ! -- add supplemental pole coordinates to mesh only on PETs with adjacent regions
    if (fillPoles) then
      lp = 1 ! select adjacent field lines
      id = ipe % grid % NMP * sum(numLongNodes(iHemiStart:iHemiEnd))
      do iHemi = iHemiStart, iHemiEnd
        kpStep   = numLineNodes(lp,iHemi) - 1
        kpStart  = iHemi * kpStep + 1
        kpEnd    = (1 - iHemi) * kpStep + 1
        kpStep   = 1 - 2 * iHemi
        kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
        do kpp = kpStart, kpEnd, kpStep
          id = id + 1
          nCount = nCount + 1
          nodeIds(nCount) = id
          nodeOwners(nCount) = poleOwnerPet
          kp = kpp + kpOffset
          nodeCoords(i + 1) = rad2deg * poleCoords(kp,1)
          nodeCoords(i + 2) = 90._ESMF_KIND_R8 - rad2deg * poleCoords(kp,2)
          nodeCoords(i + 3) = 1.e-03_ESMF_KIND_R8 * ipe % grid % altitude(kp,lp)
          i = i + 3
        end do
      end do
    end if

    ! -- free up memory
    if (fillFlag) then
      deallocate(poleCoords, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

    deallocate(local_grid_3d, petMap, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- create mesh object
    mesh = ESMF_MeshCreate(parametricDim=3, spatialDim=3,  &
      coordSys=ESMF_COORDSYS_SPH_DEG, &
      rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- add nodes to mesh
    call ESMF_MeshAddNodes(mesh, nodeIds=nodeIds, &
      nodeCoords=nodeCoords, nodeOwners=nodeOwners, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- free up memory
    deallocate(nodeIds, nodeCoords, nodeOwners, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- create mesh elements
    ! -- allocate work arrays
    allocate(numLineElems(ipe % grid % NLP-1), idNodeOffset(lps:lpu,iHemiStart:iHemiEnd), stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    numLineElems = 0
    numLongElems = 0
    do lp = 1, ipe % grid % NLP - 1
      numLineElems(lp) = max(maxval(numLineNodes(lp:lp+1,0)) - 1, 0)
      numLongElems = numLongElems + numLineElems(lp)
    end do

    numElems = (iHemiEnd-iHemiStart+1)* (mpe-mps+1) * sum(numLineElems(lps:lpu-1))

    ! -- add supplemental elements if filling empty regions
    if (fillEquator) numElems = numElems + (mpe-mps+1) * (numLineNodes(ipe%grid%NLP,iHemiStart) - 1)
    if (fillPoles)   numElems = numElems + (mpe-mps+1) * (sum(numLineNodes(1,iHemiStart:iHemiEnd)) - 2)

    ! -- allocate storage for mesh element definition
    allocate(elemIds(numElems), &
      elemTypes(numElems), &
      elemConn(8*numElems), &
      stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    elemIds   = 0
    elemConn  = 0
    elemTypes = ESMF_MESHELEMTYPE_HEX

    ! -- recompute for local nodes
    numLongNodes = sum(numLineNodes(lps:lpu,:), dim=1)
    numHemiNodes = (mpu-mps+1) * numLongNodes(0)

    idNodeOffset = 0
    do lp = lps, lpu - 1
      idNodeOffset(lp + 1,:) = idNodeOffset(lp,:) + numLineNodes(lp,:)
    end do

    ! -- define three-dimensional mesh elements
    j = 0
    nCount = 0
    lpStart = sum(numLineElems(:lps-1))
    do iHemi = iHemiStart, iHemiEnd
      mpOffset = (mps-1) * numLongElems + (iHemi-iHemiStart) * ipe % grid % NMP * numLongElems
      do mp = mps, mpe
        im = mp-mps+1
        im = im * min(1,ipe % grid % NMP-im)
        lpOffset = lpStart
        do lp = lps, lpu - 1
          kpOffset = lpOffset + mpOffset
          do kp = 1, numLineElems(lp)
            nCount = nCount + 1
            elemIds(nCount) = kp + kpOffset
            do i = 1, 7, 2
              j = j + 1
              kpp = kp + iConn(i)
              lpp = lp + iConn(i+1)
              kpe = min(kpp, numLineNodes(lpp,0))
              jHemi = iHemi * max(min(numLineNodes(lpp,iHemi)-kpe,0)+1,0) !  select hemisphere for connected nodes
              id = min(kpp,numLineNodes(lpp,jHemi)) + idNodeOffset(lpp,jHemi) + jHemi*(iHemi-iHemiStart)*numHemiNodes
              elemConn(j  ) = id + (mp-mps)*numLongNodes(jHemi)
              elemConn(j+4) = id + im * numLongNodes(jHemi)
            end do
            j = j + 4
          end do
          lpOffset = lpOffset + numLineElems(lp)
        end do
        mpOffset = mpOffset + numLongElems
      end do
    end do

    ! -- add supplemental equatorial cells, if requested
    if (fillEquator) then
      lp = ipe%grid%NLP  ! select adjacent field line
      kpEnd = numLineNodes(lp, iHemiStart) - 1
      ! -- compute offset for global element Ids
      kpOffset = (iHemiEnd-iHemiStart+1) * ipe % grid % NMP * numLongElems + (mps-1) * kpEnd
      do mp = mps, mpe
        im = mp-mps+1
        im = im * min(1,ipe % grid % NMP-im)
        do kp = 1, kpEnd
          nCount = nCount + 1
          elemIds(nCount) = kp + kpOffset
          do i = 1, 7, 2
            j = j + 1
            kpp = kp + iConn(i)
            jHemi = iConn(i+1)
            id = min(kpp,numLineNodes(lp,jHemi)) + idNodeOffset(lp,jHemi) + jHemi*numHemiNodes
            elemConn(j  ) = id + (mp-mps)*numLongNodes(jHemi)
            elemConn(j+4) = id + im * numLongNodes(jHemi)
          end do
          j = j + 4
        end do
        kpOffset = kpOffset + kpEnd
      end do
    end if

    ! -- add supplemental polar cells, if requested
    if (fillPoles) then
      lp = 1  ! select adjacent field line
      ! -- compute offset for polar nodes local Ids
      kpe = (mpu-mps+1) * sum(numLineNodes(lps:lpu,iHemiStart:iHemiEnd))
      ! -- compute offset for global element Ids
      kpOffset = (iHemiEnd-iHemiStart+1) * ipe % grid % NMP * numLongElems
      if (fillFlag) kpOffset = kpOffset + ipe % grid % NMP * (numLineNodes(ipe%grid%NLP, iHemiStart)-1)
      do iHemi = iHemiStart, iHemiEnd
        kpOffset = kpOffset + (mps-1) * (numLineNodes(lp, iHemi)-1)
        do mp = mps, mpe
          im = mp-mps+1
          im = im * min(1,ipe % grid % NMP-im)
          do kp = 1, numLineNodes(lp, iHemi) - 1
            nCount = nCount + 1
            elemIds(nCount) = kp + kpOffset
            do i = 1, 7, 2
              j = j + 1
              kpp = kp + iConn(i)
              id = kpp + idNodeOffset(lp,iHemi) + iHemi*numHemiNodes
              elemConn(j  ) = (1-iConn(i+1)) * (kpp + kpe) + iConn(i+1) * (id + (mp-mps) * numLongNodes(iHemi))
              elemConn(j+4) = (1-iConn(i+1)) * (kpp + kpe) + iConn(i+1) * (id +       im * numLongNodes(iHemi))
            end do
            j = j + 4
          end do
          kpOffset = kpOffset + (numLineNodes(lp, iHemi)-1)
        end do
        kpOffset = kpOffset + ipe % grid % NMP * (numLineNodes(lp, iHemiStart)-1)
        kpe = kpe + numLineNodes(lp, iHemiStart)
      end do
    end if

    ! -- free up memory
    deallocate(numLineElems, idNodeOffset, numLongNodes, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- add elements to mesh
    call ESMF_MeshAddElements(mesh, elementIds=elemIds,&
           elementTypes=elemTypes, elementConn=elemConn, &
           rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- free up memory
    deallocate(elemIds, elemTypes, elemConn, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    ! -- store map connecting local nodes to indices of internal arrays
    ! -- Note: internal arrays do not include poles
    numLocalNodes = (mpe-mps+1) * sum(numLineNodes(lps:lpe,iHemiStart:iHemiEnd))
    allocate(this % nodeToIndexMap(numLocalNodes,3), stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    id = 0
    do iHemi = iHemiStart, iHemiEnd
      do mp = mps, mpe
        do lp = lps, lpe
          kpStep   = numLineNodes(lp,iHemi) - 1
          kpStart  = iHemi * kpStep + 1
          kpEnd    = (1 - iHemi) * kpStep + 1
          kpStep   = 1 - 2 * iHemi
          kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
          do kpp = kpStart, kpEnd, kpStep
            kp = kpp + kpOffset
            id = id + 1
            this % nodeToIndexMap(id, 1) = kp
            this % nodeToIndexMap(id, 2) = lp
            this % nodeToIndexMap(id, 3) = mp
          end do
        end do
      end do
    end do

    ! -- free up memory
    deallocate(jmin, jSouth, numLineNodes, stat=localrc)
    if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

  contains

    subroutine IPEGetHalo(lUBound, lHaloSize, mUBound, mHaloSize, rc)
      integer,    optional, intent(out) :: lUBound
      integer,    optional, intent(out) :: lHaloSize
      integer,    optional, intent(out) :: mUBound
      integer,    optional, intent(out) :: mHaloSize
      integer,    optional, intent(out) :: rc
  
      ! -- local variables
      integer :: lpu, lHalo, mpu, mHalo
  
      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS
  
      ! -- set bounds for local DE
      mHalo = min(ipe % grid % NMP-mpe+mps-1,1)   ! magnetic longitude: halo size
      mpu   = mpe + mHalo            ! magnetic longitude: periodic dimension
      lpu   = min(lpe + 1, ipe % grid % NLP)      ! magnetic latitude:  non-periodic dimension
      lHalo = min(ipe % grid % NLP-lpe+lps-1,1)   ! magnetic latituse:  halo size
  
      if (present(lUBound))     lUBound = lpu
      if (present(mUBound))     mUBound = mpu
      if (present(lHaloSize)) lHaloSize = lHalo
      if (present(mHaloSize)) mHaloSize = mHalo
  
    end subroutine IPEGetHalo
  
    subroutine IPEGetGridCoord(vm, coord, globalBounds, MaxFluxTube, &
                               lps, lpu, lHalo, mps, mpu, mHalo, diagnostic, rc)
      type(ESMF_VM),            intent(in)  :: vm
      real(ESMF_KIND_R4),       intent(out) :: coord(MaxFluxTube, lps:lpu, mps:mpu, 2)
      integer, dimension(:),    intent(in)  :: globalBounds
      integer,                  intent(in)  :: MaxFluxTube
      integer,                  intent(in)  :: lps, lpu
      integer,                  intent(in)  :: lHalo
      integer,                  intent(in)  :: mps, mpu
      integer,                  intent(in)  :: mHalo
      integer,        optional, intent(in)  :: diagnostic
      integer,        optional, intent(out) :: rc
  
      ! -- local variables
      integer :: localrc
      integer :: l_diagnostic
      integer :: i, j, region, nCount, localPet
      integer :: lp, lpp, kp
      integer :: mp, mr, ms, mrs, mre, mss, mse, mpp
      integer :: srcPet, dstPet
      logical :: doSend, doRecv, update
      real(ESMF_KIND_R4), dimension(:), allocatable :: sendData, recvData
  
      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS
  
      l_diagnostic = 0
      if (present(diagnostic)) l_diagnostic = diagnostic
  
      ! -- perform halo update for horizontal coordinates
      coord = BAD_VALUE
      coord(:, lps:lpe, mps:mpe, 1) = ipe % grid % longitude(:, lps:lpe, mps:mpe)
      coord(:, lps:lpe, mps:mpe, 2) = ipe % grid % colatitude(:, lps:lpe, mps:mpe)
  
      if (mHalo > 0) then
        ! -- update periodic halo region (East to West)
        mp = mpe + mHalo
        if (mp > ipe % grid % NMP) mp = mp - ipe % grid % NMP
        srcPet = IPEGetOwnerPet(lps, mp, globalBounds)
        if (srcPet < 0) then
          call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
            msg="Could not map tile to PET", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return
        end if
        mp = mps - mHalo
        if (mp < 1) mp = mp + ipe % grid % NMP
        dstPet = IPEGetOwnerPet(lps, mp, globalBounds)
        if (dstPet < 0) then
          call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
            msg="Could not map tile to PET", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return
        end if
        nCount = 2 * MaxFluxTube * (lpe-lps+1) * mHalo
        allocate(sendData(nCount), recvData(nCount), stat=localrc)
        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        sendData = 0._ESMF_KIND_R4
        recvData = 0._ESMF_KIND_R4
        mpp = mps + mHalo - 1
        i = 0
        do j = 1, 2
          do mp = mps, mpp
            do lp = lps, lpe
              do kp = 1, MaxFluxTube
                i = i + 1
                sendData(i) = coord(kp, lp, mp, j)
              end do
            end do
          end do
        end do
        call ESMF_VMSendRecv(vm, sendData, nCount, dstPet, recvData, nCount, srcPet, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        mpp = mpe + 1
        i = 0
        do j = 1, 2
          do mp = mpp, mpu
            do lp = lps, lpe
              do kp = 1, MaxFluxTube
                i = i + 1
                coord(kp, lp, mp, j) = recvData(i)
              end do
            end do
          end do
        end do
        deallocate(sendData, recvData, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
  
      ! -- update S-N (0) and SE-NW (1) halo regions
      region = 0
      update = lHalo > 0
      doSend = (lps - lHalo >  0  )
      doRecv = (lpe + lHalo <= ipe % grid % NLP)
      do while (update)
        select case (region)
        case (0)
          nCount = 2 * MaxFluxTube * lHalo * (mpe-mps+1)
          ms  = mps
          mss = mps
          mse = mpe
          mr  = mps
          mrs = mps
          mre = mpe
        case (1)
          nCount = 2 * MaxFluxTube * lHalo * mHalo
          ms  = mps - mHalo
          if (ms < 1) ms = ms + ipe % grid % NMP
          mss = mps
          mse = mps + mHalo - 1
          mr  = mpe + mHalo
          if (mr > ipe % grid % NMP) mr = mr - ipe % grid % NMP
          mrs = mpe + 1
          mre = mpu
        end select
        ! -- update region
        if (doRecv) then
          lp = lpe + lHalo
          srcPet = IPEGetOwnerPet(lp, mr, globalBounds)
          if (srcPet < 0) then
            call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
              msg="Could not map tile to PET", &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            return
          end if
          allocate(recvData(nCount), stat=localrc)
          recvData = 0._ESMF_KIND_R4
          if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
          call ESMF_VMRecv(vm, recvData, nCount, srcPet, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
          lpp = lpe + lHalo
          i = 0
          do j = 1, 2
            do mp = mrs, mre
              do lp = lpe + 1, lpp
                do kp = 1, MaxFluxTube
                  i = i + 1
                  coord(kp, lp, mp, j) = recvData(i)
                end do
              end do
            end do
          end do
          deallocate(recvData, stat=localrc)
          if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
        end if
        if (doSend) then
          lp = lps - lHalo
          dstPet = IPEGetOwnerPet(lp, ms, globalBounds)
          if (dstPet < 0) then
            call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
              msg="Could not map tile to PET", &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            return
          end if
          allocate(sendData(nCount), stat=localrc)
          sendData = 0._ESMF_KIND_R4
          if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
          lpp = lps + lHalo - 1
          i = 0
          do j = 1, 2
            do mp = mss, mse
              do lp = lps, lpp
                do kp = 1, MaxFluxTube
                  i = i + 1
                  sendData(i) = coord(kp, lp, mp, j)
                end do
              end do
            end do
          end do
          call ESMF_VMSend(vm, sendData, nCount, dstPet, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
          deallocate(sendData, stat=localrc)
          if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
        end if
        update = update .and. (mHalo > 0) .and. (region < 1)
        region = region + 1
      end do
  
      if (btest(l_diagnostic,8)) then
        call ESMF_VMGet(vm, localPet=localPet, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
        write(6,'(" Checking tile on PET ",i0,":"2(" (",i0,":",i0,") x (",i0,":",i0,")"))') &
          localPet, lps, lpe, mps, mpe, lps, lpu, mps, mpu
        nCount = count(abs(coord-BAD_VALUE) < 0.1_ESMF_KIND_R4)
        if (nCount > 0) then
          write(6,'(" - Total number of bad points: ",i0)') nCount
          nCount = count(abs(coord(:,lps:lpe,mps:mpe,:) - BAD_VALUE) < 0.1_ESMF_KIND_R4)
          write(6,'(" - Local number of bad points/flux tube/axis: ",g0)') 0.5 * nCount / MaxFluxTube
          if (mHalo > 0) then
            nCount = count(abs(coord(:,lps:lpe,mpu:mpu,:) - BAD_VALUE) < 0.1_ESMF_KIND_R4)
            write(6,'(" - Number of bad points in Eastern halo region/flux tube/axis: ",g0)') 0.5 * nCount / MaxFluxTube
          end if
          if (lHalo > 0) then
            nCount = count(abs(coord(:,lpu:lpu,mps:mpe,:) - BAD_VALUE) < 0.1_ESMF_KIND_R4)
            write(6,'(" - Number of bad points in Southern halo region/flux tube/axis: ",g0)') 0.5 * nCount / MaxFluxTube
          end if
          if ((lHalo > 0) .and. (mHalo > 0)) then
            nCount = count(abs(coord(:,lpu:lpu,mpu:mpu,:) - BAD_VALUE) < 0.1_ESMF_KIND_R4)
            write(6,'(" - Number of bad points in Southeastern halo region/flux tube/axis: ",g0)') 0.5 * nCount / MaxFluxTube
          end if
        else
          write(6,'(" - No bad point found")')
        end if
      end if
  
    end subroutine IPEGetGridCoord
  
    subroutine IPEGetPoleCoord(vm, coord, MaxFluxTube, lps, mps, mpe, nmp, rc)
  
      type(ESMF_VM),            intent(in)  :: vm
      real(ESMF_KIND_R8),       intent(out) :: coord(MaxFluxTube, 2)
      integer,                  intent(in)  :: MaxFluxTube
      integer,                  intent(in)  :: lps
      integer,                  intent(in)  :: mps
      integer,                  intent(in)  :: mpe
      integer,                  intent(in)  :: nmp
      integer,        optional, intent(out) :: rc
  
      ! -- local variables
      integer :: localrc
      integer :: dataCount
      real(ESMF_KIND_R8), dimension(:), allocatable :: sendData, recvData
  
      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS
  
      ! -- compute buffers size
      dataCount = 2 * MaxFluxTube
  
      ! -- allocate send and receive buffers
      allocate(sendData(dataCount), recvData(dataCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
  
      ! -- use higher floating point precision to reduce impact of partial sum roundoff on final results
      sendData = 0._ESMF_KIND_R8
      recvData = 0._ESMF_KIND_R8
  
      ! -- populate send buffer if on region adjacent to pole
      if (lps == 1) then
        sendData(            1:MaxFluxTube) = sum(real(ipe%grid%longitude (1:MaxFluxTube, lps, mps:mpe), kind=ESMF_KIND_R8), dim=2)
        sendData(MaxFluxTube+1:dataCount  ) = sum(real(ipe%grid%colatitude(1:MaxFluxTube, lps, mps:mpe), kind=ESMF_KIND_R8), dim=2)
      end if
  
      ! -- sum across PETs
      call ESMF_VMAllReduce(vm, sendData, recvData, dataCount, ESMF_REDUCE_SUM, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
  
      ! -- compute 2D pole coordinates as geometric average of surrounding nodes
      coord(1:MaxFluxTube, 1) = recvData(1:MaxFluxTube) / NMP
      coord(1:MaxFluxTube, 2) = recvData(1+MaxFluxTube:dataCount) / NMP
  
      ! -- free up memory
      deallocate(sendData, recvData, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
  
    end subroutine IPEGetPoleCoord
  
    integer function IPEGetOwnerPet(lp, mp, globalBounds)
      integer, intent(in) :: lp, mp
      integer, dimension(:), intent(in) :: globalBounds
  
      integer :: pet, pos
  
      IPEGetOwnerPet = -1
      pos = 0
      do pet = 0, size(globalBounds)/4 - 1
        if (((globalBounds(pos + 1) <= mp) .and. (mp <= globalBounds(pos + 2))) .and.  &
            ((globalBounds(pos + 3) <= lp) .and. (lp <= globalBounds(pos + 4)))) then
          IPEGetOwnerPet = pet
          exit
        end if
        pos = pos + 4
      end do
  
    end function IPEGetOwnerPet

  end subroutine IPEMeshCreate

  !-----------------------------------------------------------------------------

  subroutine IPEFieldDiagnostics(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Field), pointer :: fieldList(:)
    type(ESMF_State)          :: importState
    type(ESMF_Mesh)           :: mesh
    type(ESMF_VM)             :: vm

    type(IPE_InternalState_Type)  :: is
    type(IPE_Model_Type), pointer :: this

    logical :: isConnected
    integer :: item
!   integer :: iHemi
!   integer :: kp, kpp, kpStart, kpEnd, kpStep, kpOffset
!   integer :: lp, mp, mpp
!   integer :: nCount, numOwnedNodes, spatialDim
    integer :: id, js, je, kp, lp, mp
    integer :: numLocalNodes, numOwnedNodes, spatialDim
    integer :: localrc
    integer :: verbosity, diagnostic
    integer :: fieldCount
    character(len=ESMF_MAXSTR) :: errmsg
    character(len=ESMF_MAXSTR) :: msgString
    character(len=ESMF_MAXSTR) :: name
    character(len=ESMF_MAXSTR), pointer :: connectedList(:)
    character(len=ESMF_MAXSTR), pointer :: fieldNames(:)
    real(ESMF_KIND_R8)         :: dataValue
    real(ESMF_KIND_R8), dimension(:),     pointer :: dataPtr, ownedNodeCoords
    real(ESMF_KIND_R8), dimension(:),     pointer :: localMin, localMax, globalMin, globalMax

    ! local parameters
    character(len=*), parameter :: rName = "FieldDiagnostics"

    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
      diagnostic=diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- get internal state
    nullify(this)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % model

    if (.not.associated(this % nodeToIndexMap)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg="nodeToIndexMap unavailable", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    numLocalNodes = size(this % nodeToIndexMap, 1)

    ! query the Component for its VM and importState
    call ESMF_GridCompGet(gcomp, vm=vm, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_GetStateMemberLists(importState, StandardNameList=fieldNames, &
      ConnectedList=connectedList, fieldList=fieldList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    isConnected = .false.
    if (associated(fieldNames)) isConnected = any(connectedList == "true")

    nullify(fieldList, fieldNames, connectedList, dataPtr, ownedNodeCoords)
    nullify(localMin, localMax, globalMin, globalMax)

    ! -- check values of imported fields, if requested
    if (isConnected) then

      fieldCount = size(fieldList)

      allocate(localMin(fieldCount), localMax(fieldCount), &
        globalMin(fieldCount), globalMax(fieldCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      localMin = huge(0._ESMF_KIND_R8)
      localMax = -localMin

      do item = 1, fieldCount
        if (connectedList(item) == "true") then
          ! --- get field data
          nullify(dataPtr)
          call ESMF_FieldGet(fieldList(item), farrayPtr=dataPtr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
  
          if (.not.associated(ownedNodeCoords)) then
            ! -- retrieve node coordinates
            call ESMF_FieldGet(fieldList(item), mesh=mesh, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__)) &
              return  ! bail out
  
            call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedNodes=numOwnedNodes, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__)) &
              return  ! bail out
  
            allocate(ownedNodeCoords(spatialDim * numOwnedNodes), stat=localrc)
            if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)) return
  
            call ESMF_MeshGet(mesh, ownedNodeCoords=ownedNodeCoords, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__)) &
              return  ! bail out
  
          end if
  
          js = 1
          je = spatialDim
          do id = 1, numLocalNodes
            kp = this % nodeToIndexMap(id, 1)
            lp = this % nodeToIndexMap(id, 2)
            mp = this % nodeToIndexMap(id, 3)
            dataValue = dataPtr(id)
            errmsg = ""
            if (dataValue == BAD_VALUE) then
               ! -- field array contains unfilled points
               errmsg = "BAD_VALUE"
            else if (abs(dataValue - BAD_VALUE) < 0.1_ESMF_KIND_R8) then
               ! -- field array contains approximate unfilled points
               errmsg = "approx BAD_VALUE"
            else if (.not.(dataValue == dataValue)) then
               ! -- field array contains NaN
               errmsg = "NaN"
            end if
            if (len_trim(errmsg) > 0) then
              write(msgString,'("ERROR: ",a,": ",a," found (",g20.8,") at (mp,lp,kp) = (",2(i0,","),i0,")"' &
                //'," coord: ",2(f12.6,","),g20.8)') &
                trim(fieldNames(item)), trim(errmsg), mp, lp, kp, ownedNodeCoords(js:je)
              call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_WARNING, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__)) &
                return  ! bail out
            endif
            errmsg = ""
            select case (trim(fieldNames(item)))
              case ("temp_neutral", "O_Density", "O2_Density", "N2_Density")
                ! -- these fields must contain positive values
                if (dataValue <= 0._ESMF_KIND_R8) errmsg = "value <= 0.0"
            end select
            if (len_trim(errmsg) > 0) then
              write(msgString,'("ERROR: ",a,": ",a," found (",g20.8,") at (mp,lp,kp) = (",2(i0,","),i0,")"' &
                //'," coord: ",2(f12.6,","),g20.8)') &
                trim(fieldNames(item)), trim(errmsg), dataValue, mp, lp, kp, ownedNodeCoords(js:je)
              call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_WARNING, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__)) &
                return  ! bail out
            endif
            localMin(item) = min(dataValue, localMin(item))
            localMax(item) = max(dataValue, localMax(item))
            js = js + spatialDim
            je = je + spatialDim
          end do
        end if

      end do
      nullify(dataPtr)

      if (associated(ownedNodeCoords)) then
        deallocate(ownedNodeCoords, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if

      call ESMF_VMAllReduce(vm, localMin, globalMin, fieldCount, ESMF_REDUCE_MIN, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_VMAllReduce(vm, localMax, globalMax, fieldCount, ESMF_REDUCE_MAX, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      do item = 1, fieldCount
        if (connectedList(item) == "true") then
          write(msgString,'(a,": ",a,": ",a,"[",i0,"]: min/max =",2g20.8)') &
            trim(name), rName, trim(fieldNames(item)), globalMin(item), globalMax(item)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
        end if
      end do

      deallocate(localMin, localMax, globalMin, globalMax, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

    end if

  end subroutine IPEFieldDiagnostics
 
end module ipeMethods
