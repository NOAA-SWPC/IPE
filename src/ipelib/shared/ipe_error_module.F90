module ipe_error_module

#ifdef HAVE_MPI
  use mpi
#endif

  implicit none

  integer, parameter :: IPE_SUCCESS = 0, &
                        IPE_FAILURE = 1

  integer, parameter :: IPE_LOG_MSGLEN  = 1024
  
  integer, parameter :: IPE_LOG_UNIT    = 6

  integer, parameter :: IPE_LOG_ERROR   = 1
  integer, parameter :: IPE_LOG_INFO    = 2
  integer, parameter :: IPE_LOG_WARNING = 3

  character(len=*), parameter :: IPE_ERR_MSG_DEFAULT = "Internal subroutine call returned error"
  character(len=*), parameter :: IPE_ERR_MSG_ALLOC   = "Unable to allocate memory"
  character(len=*), parameter :: IPE_ERR_MSG_DEALLOC = "Unable to free up memory"
  character(len=*), parameter :: IPE_ERR_MSG_IOFILE  = "File I/O error"
  character(len=*), parameter :: IPE_ERR_MSG_IOEOF   = "End of file/end of record occurred"

  private

  public :: IPE_SUCCESS, & 
            IPE_FAILURE
  public :: IPE_LOG_ERROR, &
            IPE_LOG_INFO,  &
            IPE_LOG_WARNING

  public :: ipe_error_check,    &
            ipe_status_check,   &
            ipe_iostatus_check, &
            ipe_alloc_check,    &
            ipe_dealloc_check,  &
            ipe_error_set,      &
            ipe_warning_log
  
contains

  subroutine ipe_log(msg, severity, unit)
    character(len=*),  intent(in) :: msg
    integer, optional, intent(in) :: severity
    integer, optional, intent(in) :: unit

    ! -- local variables
    integer :: loglevel, logunit
    integer :: comm, info, rank, ierr
    logical :: is_mpi_on
    character(len=7) :: logtype

    ! -- begin
    loglevel = IPE_LOG_INFO
    if (present(severity)) loglevel = severity

    logunit  = IPE_LOG_UNIT
    if (present(unit)) logunit = unit

    logtype = ""
    select case (loglevel)
      case (IPE_LOG_ERROR)
        logtype = "ERROR"
      case (IPE_LOG_INFO)
        logtype = "INFO"
      case (IPE_LOG_WARNING)
        logtype = "WARNING"
      case default
        write(logtype, '("LEV",i4.4)') loglevel
    end select

    is_mpi_on = .false.
    rank = 0
#ifdef HAVE_MPI
    call mpi_initialized(is_mpi_on, ierr)
    if (is_mpi_on) then
!     call mpi_comm_get_info(comm, info, ierr)
      ! -- to be fixed
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    end if
#endif
    if (is_mpi_on) then
      write(logunit, '("IPE",i3.3,": ",a," - ",a)') rank, trim(logtype), trim(msg)
    else
      write(logunit, '("IPE: ",a," - ",a)') trim(logtype), trim(msg)
    end if

  end subroutine ipe_log

  subroutine ipe_warning_log(msg, file, line)
    character(len=*), optional, intent(in) :: msg
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: line

    character(len=IPE_LOG_MSGLEN) :: linestr

    linestr = ""
    if (present(line)) write(linestr,'(i0)') line

    if (present(file)) then
      if (present(msg)) then
        call ipe_log(trim(file) // ":" // trim(linestr) &
          // "    " // trim(msg), severity=IPE_LOG_WARNING)
      else
        call ipe_log(trim(file) // ":" // trim(linestr) &
          // "    " // IPE_ERR_MSG_DEFAULT, severity=IPE_LOG_WARNING)
      end if
    else
      if (present(msg)) then
        call ipe_log(trim(linestr) // "    " // trim(msg), &
          severity=IPE_LOG_WARNING)
      else
        call ipe_log(trim(linestr) // "    " // IPE_ERR_MSG_DEFAULT, &
          severity=IPE_LOG_WARNING)
      end if
    end if

  end subroutine ipe_warning_log

  subroutine ipe_error_log(msg, file, line)
    character(len=*), optional, intent(in) :: msg
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: line

    character(len=IPE_LOG_MSGLEN) :: linestr

    linestr = ""
    if (present(line)) write(linestr,'(i0)') line

    if (present(file)) then
      if (present(msg)) then
        call ipe_log(trim(file) // ":" // trim(linestr) &
          // "    " // trim(msg), severity=IPE_LOG_ERROR)
      else
        call ipe_log(trim(file) // ":" // trim(linestr) &
          // "    " // IPE_ERR_MSG_DEFAULT, severity=IPE_LOG_ERROR)
      end if
    else
      if (present(msg)) then
        call ipe_log(trim(linestr) // "    " // trim(msg), &
          severity=IPE_LOG_ERROR)
      else
        call ipe_log(trim(linestr) // "    " // IPE_ERR_MSG_DEFAULT, &
          severity=IPE_LOG_ERROR)
      end if
    end if

  end subroutine ipe_error_log

  logical function ipe_error_check(irc, msg, file, line, rc)
    integer,                    intent(in)  :: irc
    character(len=*), optional, intent(in)  :: msg
    character(len=*), optional, intent(in)  :: file
    integer,          optional, intent(in)  :: line
    integer,          optional, intent(out) :: rc

    if (present(rc)) rc = irc

    ipe_error_check = (irc /= IPE_SUCCESS)

    if (ipe_error_check) &
      call ipe_error_log(msg=msg, file=file, line=line)

  end function ipe_error_check

  logical function ipe_status_check(status, msg, file, line, rc)
    logical,                    intent(in)  :: status
    character(len=*), optional, intent(in)  :: msg
    character(len=*), optional, intent(in)  :: file
    integer,          optional, intent(in)  :: line
    integer,          optional, intent(out) :: rc

    integer :: localrc

    localrc = IPE_FAILURE
    if (status) localrc = IPE_SUCCESS

    ipe_status_check = ipe_error_check(localrc, &
      msg=msg, file=file, line=line)

    if (present(rc)) rc = localrc

  end function ipe_status_check

  logical function ipe_alloc_check(status, msg, file, line, rc)
    integer,                    intent(in)  :: status
    character(len=*), optional, intent(in)  :: msg
    character(len=*), optional, intent(in)  :: file
    integer,          optional, intent(in)  :: line
    integer,          optional, intent(out) :: rc

    if (present(msg)) then
      ipe_alloc_check = &
        ipe_status_check(status == 0, msg=msg, &
        file=file, line=line, rc=rc)
    else
      ipe_alloc_check = &
        ipe_status_check(status == 0, msg=IPE_ERR_MSG_ALLOC, &
        file=file, line=line, rc=rc)
    end if

  end function ipe_alloc_check

  logical function ipe_dealloc_check(status, msg, file, line, rc)
    integer,                    intent(in)  :: status
    character(len=*), optional, intent(in)  :: msg
    character(len=*), optional, intent(in)  :: file
    integer,          optional, intent(in)  :: line
    integer,          optional, intent(out) :: rc

    if (present(msg)) then
      ipe_dealloc_check = &
        ipe_status_check(status == 0, msg=msg, &
        file=file, line=line, rc=rc)
    else
      ipe_dealloc_check = &
        ipe_status_check(status == 0, msg=IPE_ERR_MSG_DEALLOC, &
        file=file, line=line, rc=rc)
    end if

  end function ipe_dealloc_check

  logical function ipe_iostatus_check(status, msg, file, line, rc)
    integer,                    intent(in)  :: status
    character(len=*), optional, intent(in)  :: msg
    character(len=*), optional, intent(in)  :: file
    integer,          optional, intent(in)  :: line
    integer,          optional, intent(out) :: rc

    if (present(msg)) then
      ipe_iostatus_check = &
        ipe_status_check(status == 0, msg=msg, &
        file=file, line=line, rc=rc)
    else
      ipe_iostatus_check = .false.
      if (status < 0) then
        ipe_iostatus_check = &
          ipe_status_check(.true., msg=IPE_ERR_MSG_IOEOF, &
          file=file, line=line, rc=rc)
      else if (status > 0) then
        ipe_iostatus_check = &
          ipe_status_check(.true., msg=IPE_ERR_MSG_IOFILE, &
          file=file, line=line, rc=rc)
      end if
    end if

  end function ipe_iostatus_check

  subroutine ipe_error_set(msg, file, line, rc)
    character(len=*), optional, intent(in)  :: msg
    character(len=*), optional, intent(in)  :: file
    integer,          optional, intent(in)  :: line
    integer,          optional, intent(out) :: rc

    logical :: status

    status = ipe_status_check(.false., msg=msg, &
        file=file, line=line, rc=rc)

  end subroutine ipe_error_set

end module ipe_error_module
