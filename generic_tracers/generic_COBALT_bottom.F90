module generic_COBALT_bottom

  use time_manager_mod,  only: time_type
  use field_manager_mod, only: fm_string_len

  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use g_tracer_utils, only : g_tracer_type, g_tracer_get_common
  use g_tracer_utils, only : g_tracer_set_values, g_tracer_get_values
  use g_tracer_utils, only : g_tracer_get_pointer, g_send_data

implicit none ; private

public generic_COBALT_update_from_bottom_simple_slab
public CBED_update_from_bottom
public allocate_cobalt_btm
public deallocate_cobalt_btm
public generic_COBALT_btm_register_diag
public generic_COBALT_btm_update_from_source

type :: COBALT_btm_type
  real, dimension(:,:), ALLOCATABLE, private :: &
       fcadet_arag_btm, &
       fcadet_calc_btm, & ! Used in L 11445
       ffedet_btm, & ! used in L 11489
       flithdet_btm, &
       fndet_btm, &
       fpdet_btm, &
       fsidet_btm
  integer, private :: &
       id_fcadet_arag_btm = -1, &
       id_fcadet_calc_btm = -1, &
       id_ffedet_btm      = -1, &
       id_flithdet_btm    = -1, &
       id_fndet_btm       = -1, &
       id_fpdet_btm       = -1, &
       id_fsidet_btm      = -1
  contains
    procedure, public :: get_fcadet_calc_btm
    procedure, public :: get_ffedet_btm
end type COBALT_btm_type

type vardesc
  character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
  character(len=fm_string_len) :: longname ! The long name of that variable.
  character(len=1)  :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
  character(len=1)  :: z_grid   ! The vert. grid:  L, i, or 1.
  character(len=1)  :: t_grid   ! The time description: s, a, m, or 1.
  character(len=fm_string_len) :: units    ! The dimensions of the variable.
  character(len=1)  :: mem_size ! The size in memory: d or f.
end type vardesc

type(COBALT_btm_type) :: cobalt_btm

contains

  subroutine get_fcadet_calc_btm(self, fcadet_calc_btm)
    class(COBALT_btm_type), intent(in) :: self
    real, dimension(:,:), pointer, intent(out) :: fcadet_calc_btm

    fcadet_calc_btm = self%fcadet_calc_btm

  end subroutine get_fcadet_calc_btm


  subroutine get_ffedet_btm(self, ffedet_btm)
    class(COBALT_btm_type), intent(in) :: self
    real, dimension(:,:), pointer, intent(out) :: ffedet_btm

    ffedet_btm = self%ffedet_btm

  end subroutine get_ffedet_btm

  subroutine generic_COBALT_btm_register_diag(package_name, missing_value1, axes, init_time)
    character(len=*), intent(in) :: package_name
    real            , intent(in) :: missing_value1
    integer         , intent(in) :: axes(3)
    type(time_type) , intent(in) :: init_time

    type(vardesc) :: vardesc_temp

    vardesc_temp = vardesc("fcadet_arag_btm","CaCO3 sinking flux at bottom",'h','1','s','mol m-2 s-1','f')
    cobalt_btm%id_fcadet_arag_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_btm","CaCO3 sinking flux at bottom",'h','1','s','mol m-2 s-1','f')
    cobalt_btm%id_fcadet_calc_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_btm","fedet sinking flux burial",'h','1','s','mol m-2 s-1','f')
    cobalt_btm%id_ffedet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_btm","Lithogenic detrital sinking flux burial",'h','1','s','g m-2 s-1','f')
    cobalt_btm%id_flithdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_btm","ndet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt_btm%id_fndet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_btm","pdet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt_btm%id_fpdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_btm","sidet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt_btm%id_fsidet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
  end subroutine generic_COBALT_btm_register_diag

  subroutine generic_COBALT_btm_update_from_source(model_time, grid_tmask, isc, jsc, iec, jec)
    type(time_type)                , intent(in) :: model_time
    real, dimension(:,:,:) ,pointer, intent(in) :: grid_tmask
    integer                        , intent(in) :: isc, jsc, iec, jec

    logical :: used

    if (cobalt_btm%id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fcadet_arag_btm,   cobalt_btm%fcadet_arag_btm,      &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt_btm%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fcadet_calc_btm,   cobalt_btm%fcadet_calc_btm,      &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt_btm%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_ffedet_btm,   cobalt_btm%ffedet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt_btm%id_fndet_btm .gt. 0)            &
         used = g_send_data(cobalt_btm%id_fndet_btm,    cobalt_btm%fndet_btm,              &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt_btm%id_fpdet_btm .gt. 0)            &
         used = g_send_data(cobalt_btm%id_fpdet_btm,    cobalt_btm%fpdet_btm,              &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt_btm%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fsidet_btm,   cobalt_btm%fsidet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt_btm%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_flithdet_btm,   cobalt_btm%flithdet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  end subroutine generic_COBALT_btm_update_from_source

  subroutine allocate_cobalt_btm(isd, ied, jsd, jed)
    integer, intent(in) :: isd, ied, jsd, jed

    allocate(cobalt_btm%fcadet_arag_btm(isd:ied, jsd:jed))    ; cobalt_btm%fcadet_arag_btm=0.0
    allocate(cobalt_btm%fcadet_calc_btm(isd:ied, jsd:jed))    ; cobalt_btm%fcadet_calc_btm=0.0
    allocate(cobalt_btm%ffedet_btm(isd:ied, jsd:jed))         ; cobalt_btm%ffedet_btm=0.0
    allocate(cobalt_btm%flithdet_btm(isd:ied, jsd:jed))       ; cobalt_btm%flithdet_btm=0.0
    allocate(cobalt_btm%fndet_btm(isd:ied, jsd:jed))          ; cobalt_btm%fndet_btm=0.0
    allocate(cobalt_btm%fpdet_btm(isd:ied, jsd:jed))          ; cobalt_btm%fpdet_btm=0.0
    allocate(cobalt_btm%fsidet_btm(isd:ied, jsd:jed))         ; cobalt_btm%fsidet_btm=0.0
  end subroutine allocate_cobalt_btm

  subroutine deallocate_cobalt_btm()

    deallocate(cobalt_btm%fcadet_arag_btm)
    deallocate(cobalt_btm%fcadet_calc_btm)
    deallocate(cobalt_btm%ffedet_btm)
    deallocate(cobalt_btm%flithdet_btm)
    deallocate(cobalt_btm%fndet_btm)
    deallocate(cobalt_btm%fpdet_btm)
    deallocate(cobalt_btm%fsidet_btm)
  end subroutine deallocate_cobalt_btm

  subroutine generic_COBALT_update_from_bottom_simple_slab(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    !
    ! The bottom reservoirs of aragonite and calcite are immediately redistributed to the
    ! water column as a bottom flux (btf) where they impact the alkalinity and DIC
    !
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',cobalt_btm%fcadet_arag_btm,isd,jsd)
    cobalt_btm%fcadet_arag_btm = cobalt_btm%fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (cobalt_btm%id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fcadet_arag_btm,cobalt_btm%fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',cobalt_btm%fcadet_calc_btm,isd,jsd)
    cobalt_btm%fcadet_calc_btm = cobalt_btm%fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (cobalt_btm%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fcadet_calc_btm, cobalt_btm%fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',cobalt_btm%ffedet_btm,isd,jsd)
    cobalt_btm%ffedet_btm = cobalt_btm%ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = cobalt_btm%ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (cobalt_btm%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_ffedet_btm, cobalt_btm%ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',cobalt_btm%flithdet_btm,isd,jsd)
    cobalt_btm%flithdet_btm = cobalt_btm%flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (cobalt_btm%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_flithdet_btm, cobalt_btm%flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',cobalt_btm%fndet_btm,isd,jsd)
    cobalt_btm%fndet_btm = cobalt_btm%fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (cobalt_btm%id_fndet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fndet_btm,cobalt_btm%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',cobalt_btm%fpdet_btm,isd,jsd)
    cobalt_btm%fpdet_btm = cobalt_btm%fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (cobalt_btm%id_fpdet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fpdet_btm,cobalt_btm%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',cobalt_btm%fsidet_btm,isd,jsd)
    cobalt_btm%fsidet_btm = cobalt_btm%fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (cobalt_btm%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fsidet_btm,    cobalt_btm%fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine generic_COBALT_update_from_bottom_simple_slab

  subroutine CBED_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    !
    ! The bottom reservoirs of aragonite and calcite are immediately redistributed to the
    ! water column as a bottom flux (btf) where they impact the alkalinity and DIC
    !
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',cobalt_btm%fcadet_arag_btm,isd,jsd)
    cobalt_btm%fcadet_arag_btm = cobalt_btm%fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (cobalt_btm%id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fcadet_arag_btm,cobalt_btm%fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',cobalt_btm%fcadet_calc_btm,isd,jsd)
    cobalt_btm%fcadet_calc_btm = cobalt_btm%fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (cobalt_btm%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fcadet_calc_btm, cobalt_btm%fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',cobalt_btm%ffedet_btm,isd,jsd)
    cobalt_btm%ffedet_btm = cobalt_btm%ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = cobalt_btm%ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (cobalt_btm%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_ffedet_btm, cobalt_btm%ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',cobalt_btm%flithdet_btm,isd,jsd)
    cobalt_btm%flithdet_btm = cobalt_btm%flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (cobalt_btm%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_flithdet_btm, cobalt_btm%flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',cobalt_btm%fndet_btm,isd,jsd)
    cobalt_btm%fndet_btm = cobalt_btm%fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (cobalt_btm%id_fndet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fndet_btm,cobalt_btm%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',cobalt_btm%fpdet_btm,isd,jsd)
    cobalt_btm%fpdet_btm = cobalt_btm%fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (cobalt_btm%id_fpdet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fpdet_btm,cobalt_btm%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',cobalt_btm%fsidet_btm,isd,jsd)
    cobalt_btm%fsidet_btm = cobalt_btm%fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt_btm%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (cobalt_btm%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt_btm%id_fsidet_btm,    cobalt_btm%fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine CBED_update_from_bottom

end module generic_COBALT_bottom
