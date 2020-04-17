!nm20140919: originally copied from MODEL.F90 but modified for IPE
! DATE: 19 September, 2014
!********************************************
!***      Copyright 2014 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Whole Atmosphere Model(WAM)-Ionosphere Plasmasphere Electrodynamics(IPE) model Interface
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!-------------------------------------------- 
!
module ipeCap

  !-----------------------------------------------------------------------------
  ! IPE Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS           => SetServices,          &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance,        &
    model_label_Finalize       => label_Finalize

  use ipeMethods
  use IPE_Wrapper


  implicit none

  integer, parameter :: importFieldCount = 7
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "temp_neutral          ", &
      "eastward_wind_neutral ", &
      "northward_wind_neutral", &
      "upward_wind_neutral   ", &
      "O_Density             ", &
      "O2_Density            ", &
      "N2_Density            "  &
      /)
  
  integer, parameter :: exportFieldCount = 21
  character(len=*), dimension(exportFieldCount), parameter :: &
    exportFieldNames = (/ &
      "temp_neutral          ", &
      "eastward_wind_neutral ", &
      "northward_wind_neutral", &
      "upward_wind_neutral   ", &
      "O_Density             ", &
      "O2_Density            ", &
      "N2_Density            ", &
      "O_plus_density        ", &
      "H_plus_density        ", &
      "He_plus_density       ", &
      "N_plus_density        ", &
      "NO_plus_density       ", &
      "O2_plus_density       ", &
      "N2_plus_density       ", &
      "O_plus_2D_density     ", &
      "O_plus_2P_density     ", &
      "ion_temperature       ", &
      "electron_temperature  ", &
      "eastward_exb_velocity ", &
      "northward_exb_velocity", &
      "upward_exb_velocity   "  &
      /)
  
  private

  public :: SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! begin
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv01
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p3"/), userRoutine=InitializeP3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, specRoutine=InitializeData, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! -- advance method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- finalize method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
     specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: name

    ! local parameters
    character(len=*), parameter :: rName = "InitializeP0"
    
    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv02p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine InitializeP0

  !-----------------------------------------------------------------------------

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: name

    ! local parameters
    character(len=*), parameter :: rName = "InitializeP1"

    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! import fields from WAM
    call NUOPC_Advertise(importState, StandardNames=importFieldNames, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! export fields from WAM
    call NUOPC_Advertise(exportState, StandardNames=exportFieldNames, &
      SharePolicyField="share", TransferOfferGeomObject="will provide", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP1
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field) :: field
    type(ESMF_Mesh)  :: mesh
    type(ESMF_VM)    :: vm

    type(IPE_InternalState_Type) :: is
    type(IPE_Model), pointer     :: ipe

    integer :: item, stat
    integer :: verbosity
    logical :: isImportConnected, isExportConnected
    character(len=ESMF_MAXSTR) :: name

    ! local parameters
    character(len=*), parameter :: rName = "InitializeP3"
    
    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! check if fields are connected, i.e. if the model is coupled
    ! - import state
    isImportConnected = IPEIsStateConnected(importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    ! - export state
    isExportConnected = IPEIsStateConnected(exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! allocate memory for the internal state and store it into component
    allocate(is % model, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    allocate(is % model % ipe, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    ! initialize internal state
    nullify(is % model % nodeToIndexMap)

    ipe => is % model % ipe

    ! Get local VM
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! initialize IPE
    CALL Initialize_IPE(ipe, clock, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set IPE coupled flag
    ipe % forcing % coupled = isImportConnected

    if (isImportConnected .or. isExportConnected) then
      ! create 3D IPE mesh
      call IPEMeshCreate(gcomp, mesh, fill=(ipe % parameters % mesh_fill > 0), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! write IPE mesh to VTK file if requested
      if (ipe % parameters % mesh_write > 0) then
        call ESMF_MeshWrite(mesh, trim(ipe % parameters % mesh_write_file), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

    ! realize connected Fields in the importState
    do item = 1, importFieldCount
      call NUOPC_Realize(importState, mesh, fieldName=trim(importFieldNames(item)), &
        typekind=ESMF_TYPEKIND_R8, selection="realize_connected_remove_others", &
        dataFillScheme="const", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    ! realize connected Fields in the exportState
    do item = 1, exportFieldCount
      call NUOPC_Realize(exportState, mesh, fieldName=trim(exportFieldNames(item)), &
        typekind=ESMF_TYPEKIND_R8, selection="realize_connected_remove_others", &
        dataFillScheme="const", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP3
  
  !-----------------------------------------------------------------------------

  subroutine InitializeData(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_State)     :: importState

    integer :: verbosity
    character(len=ESMF_MAXSTR) :: name

    ! local parameters
    character(len=*), parameter :: rName = "DataInitialize"

    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
        
    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeData

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)   :: clock
    type(ESMF_State)   :: importState, exportState
    type(ESMF_Field), pointer :: fieldList(:)

    type(IPE_InternalState_Type)  :: is
    type(IPE_Model_Type), pointer :: this
    type(IPE_Model),      pointer :: ipe

    integer :: id, item, stat
    integer :: kps, kpe, lps, lpe, mph, mps, mpe
    integer :: kp, lp, mp
    integer :: numLocalNodes
    integer :: verbosity, diagnostic
    integer(ESMF_KIND_R8) :: advanceCount
    character(len=ESMF_MAXSTR) :: name
    character(len=ESMF_MAXSTR), pointer :: standardNameList(:)
    character(len=ESMF_MAXSTR), pointer :: connectedList(:)
    real(ESMF_KIND_R8), dimension(:),     pointer :: fieldPtr
    real(prec),         dimension(:,:,:), pointer :: modelPtr

    ! local parameters
    character(len=*), parameter :: rName = "Run"

    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, &
      diagnostic=diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! HERE IPE ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.

    ! -- do not import fields at startup since they are not available
    call ESMF_ClockGet(clock, advanceCount=advanceCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- get internal state
    nullify(this, ipe)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % model
    ipe  => is % model % ipe

    if (.not.associated(this % ipe)) then
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg="IPE model unavailable", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    kps = 1
    lps = ipe % mpi_layer % lp_low
    mps = ipe % mpi_layer % mp_low
    mph = mps - ipe % mpi_layer % mp_halo_size

    ! -- import data
    nullify(fieldList, connectedList, standardNameList)

    if (advanceCount > 0) then
      call NUOPC_GetStateMemberLists(importState, &
        StandardNameList=standardNameList, ConnectedList=connectedList, &
        fieldList=fieldList, nestedFlag=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
    end if

    if (associated(fieldList)) then

      if (.not.associated(this % nodeToIndexMap)) then
        call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
          msg="nodeToIndexMap unavailable", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if

      numLocalNodes = size(this % nodeToIndexMap, 1)

      do item = 1, size(fieldList)
        ! --- get field data
        nullify(fieldPtr)
        call ESMF_FieldGet(fieldList(item), farrayPtr=fieldPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out

        ! -- identify IPE neutral array receiving imported field data
        nullify(modelPtr)
        select case (trim(standardNameList(item)))
          case ("temp_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % temperature
          case ("eastward_wind_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % velocity_geographic(1,:,:,:)
          case ("northward_wind_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % velocity_geographic(2,:,:,:)
          case ("upward_wind_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % velocity_geographic(3,:,:,:)
          case ("O_Density")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % oxygen
          case ("O2_Density")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % molecular_oxygen
          case ("N2_Density")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % molecular_nitrogen
          case default
            ! -- unavailable neutrals array, skip it
            cycle
        end select

        do id = 1, numLocalNodes
          kp = this % nodeToIndexMap(id, 1)
          lp = this % nodeToIndexMap(id, 2)
          mp = this % nodeToIndexMap(id, 3)
          modelPtr(kp, lp, mp) = fieldPtr(id)
        end do

      end do

      ! -- check values of imported fields, if requested
      if (btest(diagnostic,17)) then
        call IPEFieldDiagnostics(gcomp, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if

    end if

    if (associated(fieldList)) then
      deallocate(fieldList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if
    if (associated(connectedList)) then
      deallocate(connectedList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if
    if (associated(standardNameList)) then
      deallocate(standardNameList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

    ! -- advance IPE model
    call Update_IPE(ipe, clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- export data
    nullify(fieldList, connectedList, standardNameList)
    call NUOPC_GetStateMemberLists(exportState, &
      StandardNameList=standardNameList, ConnectedList=connectedList, &
      fieldList=fieldList, nestedFlag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (associated(fieldList)) then

      if (.not.associated(this % nodeToIndexMap)) then
        call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
          msg="nodeToIndexMap unavailable", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if

      numLocalNodes = size(this % nodeToIndexMap, 1)

      do item = 1, size(fieldList)
        ! --- get field data
        nullify(fieldPtr)
        call ESMF_FieldGet(fieldList(item), farrayPtr=fieldPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out

        ! -- identify IPE neutral array receiving imported field data
        nullify(modelPtr)
        select case (trim(standardNameList(item)))
          case ("temp_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % temperature
          case ("eastward_wind_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % velocity_geographic(1,:,:,:)
          case ("northward_wind_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % velocity_geographic(2,:,:,:)
          case ("upward_wind_neutral")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % velocity_geographic(3,:,:,:)
          case ("O_Density")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % oxygen
          case ("O2_Density")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % molecular_oxygen
          case ("N2_Density")
            modelPtr(kps:,lps:,mps:) => ipe % neutrals % molecular_nitrogen
          case ("ion_temperature")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_temperature
          case ("electron_temperature")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % electron_temperature
          case ("O_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(1,:,:,:)
          case ("H_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(2,:,:,:)
          case ("He_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(3,:,:,:)
          case ("N_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(4,:,:,:)
          case ("NO_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(5,:,:,:)
          case ("O2_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(6,:,:,:)
          case ("N2_plus_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(7,:,:,:)
          case ("O_plus_2D_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(8,:,:,:)
          case ("O_plus_2P_density")
            modelPtr(kps:,lps:,mph:) => ipe % plasma % ion_densities(9,:,:,:)
          case ("eastward_exb_velocity")
            modelPtr(kps:,lps:,mps:) => ipe % eldyn % v_exb_geographic(1,:,:,:)
          case ("northward_exb_velocity")
            modelPtr(kps:,lps:,mps:) => ipe % eldyn % v_exb_geographic(2,:,:,:)
          case ("upward_exb_velocity")
            modelPtr(kps:,lps:,mps:) => ipe % eldyn % v_exb_geographic(3,:,:,:)
          case default
            ! -- unavailable neutrals array, skip it
            cycle
        end select

        do id = 1, numLocalNodes
          kp = this % nodeToIndexMap(id, 1)
          lp = this % nodeToIndexMap(id, 2)
          mp = this % nodeToIndexMap(id, 3)
          fieldPtr(id) = modelPtr(kp, lp, mp)
        end do
      end do

    end if

    if (associated(fieldList)) then
      deallocate(fieldList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if
    if (associated(connectedList)) then
      deallocate(connectedList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if
    if (associated(standardNameList)) then
      deallocate(standardNameList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelAdvance
 
  !-----------------------------------------------------------------------------

  subroutine Finalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! -- local variables
    integer :: localrc, stat
    integer :: verbosity
    character(len=ESMF_MAXSTR) :: name

    type(IPE_InternalState_Type)  :: is
    type(IPE_Model_Type), pointer :: this

    ! local parameters
    character(len=*), parameter :: rName = "Finalize"

    ! begin
    rc = ESMF_SUCCESS

    ! Get component information
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get internal state
    nullify(this)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % model

    if (associated(this % ipe)) then
      ! finalize IPE model
      call Finalize_IPE(this % ipe, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      deallocate(this % ipe, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      nullify(this % ipe)
    end if

    ! deallocate work memory
    if (associated(this % nodeToIndexMap)) then
      deallocate(this % nodeToIndexMap, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      nullify(this % nodeToIndexMap)
    end if

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine Finalize

end module ipeCap
