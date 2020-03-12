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
#define ESMF_CONTEXT  line=__LINE__,file=__FILE__,method=ESMF_METHOD

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
  use IPE_Constants_Dictionary


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

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::SetServices()"

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! begin
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv01
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p3"/), userRoutine=InitializeP3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, specRoutine=InitializeData, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    
    ! -- advance method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! -- finalize method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
     specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

  end subroutine SetServices
  
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeP0()"

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
      ESMF_CONTEXT)) &
      return  ! bail out

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine InitializeP0

  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeP1()"

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

    rc = ESMF_SUCCESS
    
    ! import fields from WAM
    call NUOPC_Advertise(importState, StandardNames=importFieldNames, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! export fields from WAM
    call NUOPC_Advertise(exportState, StandardNames=exportFieldNames, &
      SharePolicyField="share", TransferOfferGeomObject="will provide", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP1
  
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeP3()"

  subroutine InitializeP3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field) :: field
    type(ESMF_Mesh)  :: mesh

    integer :: item
    integer :: verbosity
    logical :: isConnected
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

    ! check if all required fields are connected
    item = 0
    isConnected = .true.
    do while (isConnected .and. (item < importFieldCount))
      item = item + 1
      isConnected = NUOPC_IsConnected(importState, &
        fieldName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    if (item < importFieldCount) then
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="Not all required fields are connected", &
        ESMF_CONTEXT, &
        rcToReturn=rc)
      return
    end if

    ! initialize IPE
    CALL Initialize_IPE(clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

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

    ! realize connected Fields in the importState
    do item = 1, importFieldCount
      field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
        name=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldFill(field, dataFillScheme="const", &
        const1=0._ESMF_KIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_Realize(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    ! realize connected Fields in the exportState
    do item = 1, exportFieldCount
      isConnected = NUOPC_IsConnected(exportState, &
        fieldName=trim(exportFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (isConnected) then
        field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
          name=trim(exportFieldNames(item)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end if
    end do

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP3
  
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeData()"

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
      ESMF_CONTEXT)) &
      return  ! bail out
        
    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeData

  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::ModelAdvance()"

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)   :: clock
    type(ESMF_State)   :: importState, exportState
    type(ESMF_Field)   :: field, efield
    type(ESMF_Mesh)    :: mesh
    type(ESMF_VM)      :: vm

    integer :: i, item
    integer :: iHemi
    integer :: kp, kpp, kpStart, kpEnd, kpStep, kpOffset
    integer :: lp, mp, mpp
    integer :: nCount
    integer :: localrc, localPet
    integer :: verbosity, diagnostic
    integer(ESMF_KIND_R8) :: advanceCount
    character(len=ESMF_MAXSTR) :: errmsg
    character(len=ESMF_MAXSTR) :: name
    real(ESMF_KIND_R8), dimension(:),     pointer :: dataPtr
    real(prec),         dimension(:,:,:), pointer :: neutralsPtr

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
      ESMF_CONTEXT)) &
      return  ! bail out

    ! HERE IPE ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.

    ! -- do not import fields at startup since they are not available
    call ESMF_ClockGet(clock, advanceCount=advanceCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out


    if (advanceCount /= 0) then

      ! -- import data

      do item = 1, importFieldCount
        ! --- retrieve field
        call ESMF_StateGet(importState, field=field, &
          itemName=trim(importFieldNames(item)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          ESMF_CONTEXT)) &
          return  ! bail out

        ! --- get field data
        nullify(dataPtr)
        call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          ESMF_CONTEXT)) &
          return  ! bail out

        ! -- check if data size is consistent with local mesh size
        if (size(dataPtr) /= numLocalNodes) then
          errmsg = ""
          write(errmsg, '("Field: ",a,": retrieved data size (",i0,&
            &") does not match expected data size (",i0,")")') &
            trim(importFieldNames(item)), size(dataPtr), numLocalNodes
          call ESMF_LogSetError(ESMF_RC_PTR_BAD, &
            msg=errmsg, &
            ESMF_CONTEXT, &
            rcToReturn=rc)
          return  ! bail out
        end if

        ! -- identify IPE neutral array receiving imported field data
        nullify(neutralsPtr)
        select case (trim(importFieldNames(item)))
          case ("temp_neutral")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % temperature
          case ("eastward_wind_neutral")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % velocity_geographic(1,:,:,:)
          case ("northward_wind_neutral")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % velocity_geographic(2,:,:,:)
          case ("upward_wind_neutral")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % velocity_geographic(3,:,:,:)
          case ("O_Density")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % oxygen
          case ("O2_Density")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % molecular_oxygen
          case ("N2_Density")
            neutralsPtr(1:,lps:,mps:) => ipe % neutrals % molecular_nitrogen
          case default
            ! -- unavailable neutrals array, skip it
            cycle
        end select

        nCount = 0
        do iHemi = iHemiStart, iHemiEnd
          do mp = mps, mpe
            do lp = lps, lpe
              kpStep   = numLineNodes(lp,iHemi) - 1
              kpStart  = iHemi * kpStep + 1
              kpEnd    = (1 - iHemi) * kpStep + 1
              kpStep   = 1 - 2 * iHemi
              kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
              do kpp = kpStart, kpEnd, kpStep
                nCount = nCount + 1
                kp = kpp + kpOffset
                neutralsPtr(kp,lp,mp) = dataPtr(nCount)
              end do
            end do
          end do
        end do

      end do

    end if

    ! -- check values of imported fields, if requested
    if (btest(diagnostic,8)) then
      call IPEFieldDiagnostics(gcomp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
    end if

    nullify(dataPtr)

    ! -- advance IPE model
    call Update_IPE(clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! -- copy import neutral field data to export fields for I/O purposes
    do item = 1, importFieldCount
      ! --- retrieve import field
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
      ! --- retrieve export field
      call ESMF_StateGet(exportState, field=efield, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
      ! -- copy content
      call ESMF_FieldCopy(efield, field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
    end do

    ! -- export ion densities
    do item = importFieldCount + 1, exportFieldCount
      ! --- retrieve export field
      call ESMF_StateGet(exportState, field=efield, &
        itemName=trim(exportFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out

      ! --- get field data
      nullify(dataPtr)
      call ESMF_FieldGet(efield, farrayPtr=dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out

      i = item - importFieldCount
      nCount = 0
      do iHemi = iHemiStart, iHemiEnd
        do mp = mps, mpe
          do lp = lps, lpe
            kpStep   = numLineNodes(lp,iHemi) - 1
            kpStart  = iHemi * kpStep + 1
            kpEnd    = (1 - iHemi) * kpStep + 1
            kpStep   = 1 - 2 * iHemi
            kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
            do kpp = kpStart, kpEnd, kpStep
              nCount = nCount + 1
              kp = kpp + kpOffset
              select case (trim(exportFieldNames(item)))
                case ("ion_temperature")
                  dataPtr(nCount) = ipe % plasma % ion_temperature(kp,lp,mp)
                case ("electron_temperature")
                  dataPtr(nCount) = ipe % plasma % electron_temperature(kp,lp,mp)
                case ("eastward_exb_velocity")
                  dataPtr(nCount) = ipe % eldyn % v_exb_geographic(1,kp,lp,mp)
                case ("northward_exb_velocity")
                  dataPtr(nCount) = ipe % eldyn % v_exb_geographic(2,kp,lp,mp)
                case ("upward_exb_velocity")
                  dataPtr(nCount) = ipe % eldyn % v_exb_geographic(3,kp,lp,mp)
                case default
                  dataPtr(nCount) = ipe % plasma % ion_densities(i,kp,lp,mp)
              end select
            end do
          end do
        end do
      end do

    end do

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelAdvance
 
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::Finalize()"

  subroutine Finalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: verbosity
    character(len=ESMF_MAXSTR) :: name

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

    ! -- finalize IPE model
    call Finalize_IPE(rc)
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

  end subroutine Finalize

end module ipeCap
