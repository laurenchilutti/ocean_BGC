module generic_COBALT_update_from_bottom_mod

  use time_manager_mod,  only: time_type

  use g_tracer_utils, only : g_tracer_type, g_tracer_get_common
  use g_tracer_utils, only : g_tracer_set_values, g_tracer_get_values
  use g_tracer_utils, only : g_tracer_get_pointer, g_send_data

  use generic_COBALT, only : generic_COBALT_type

implicit none ; private

public generic_COBALT_update_from_bottom_simple_slab
public CBED_update_from_bottom

contains

  ! <SUBROUTINE NAME="generic_COBALT_update_from_bottom_simple_slab">
  ! 
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !
  !   This routine calculates bottom fluxes for tracers with bottom reservoirs.
  !   It is called near the end of the time step, meaning that the fluxes 
  !   calculated pertain to the next time step.  
  ! 
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_bottom_simple_slab(cobalt,tracer_list,dt, tau, model_time) 
  !  </TEMPLATE>
  !
  !  <IN NAME="cobalt" TYPE="type(generic_COBALT_type)">
  !  </IN>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment 
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  !
  ! </SUBROUTINE>
  subroutine generic_COBALT_update_from_bottom_simple_slab(tracer_list, dt, tau, model_time, fcadet_arag_btm, id_fcadet_arag_btm, fcadet_calc_btm, id_fcadet_calc_btm, ffedet_btm, id_ffedet_btm, flithdet_btm, id_flithdet_btm, fndet_btm, id_fndet_btm, fpdet_btm, id_fpdet_btm, fsidet_btm, id_fsidet_btm)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time
    real, dimension(:,:), intent(inout):: fcadet_arag_btm
    integer,              intent(in) :: id_fcadet_arag_btm
    real, dimension(:,:), intent(inout) :: fcadet_calc_btm
    integer,              intent(in) :: id_fcadet_calc_btm
    real, dimension(:,:), intent(inout) :: ffedet_btm
    integer,              intent(in) :: id_ffedet_btm
    real, dimension(:,:), intent(inout) :: flithdet_bt
    integer,              intent(in) :: id_flithdet_btm
    real, dimension(:,:), intent(inout) :: fndet_btm
    integer,              intent(in) :: id_fndet_btm
    real, dimension(:,:), intent(inout) :: fpdet_btm
    integer,              intent(in) :: id_fpdet_btm
    real, dimension(:,:), intent(inout) :: fsidet_btm
    integer,              intent(in) :: id_fsidet_btm

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    !
    ! The bottom reservoirs of aragonite and calcite are immediately redistributed to the
    ! water column as a bottom flux (btf) where they impact the alkalinity and DIC
    !
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',fcadet_arag_btm,isd,jsd)
    fcadet_arag_btm = fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(id_fcadet_arag_btm,fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',fcadet_calc_btm,isd,jsd)
    fcadet_calc_btm = fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(id_fcadet_calc_btm, fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',ffedet_btm,isd,jsd)
    ffedet_btm = ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (id_ffedet_btm .gt. 0)           &
         used = g_send_data(id_ffedet_btm, ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',flithdet_btm,isd,jsd)
    flithdet_btm = flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (id_flithdet_btm .gt. 0)           &
         used = g_send_data(id_flithdet_btm, flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of 
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',fndet_btm,isd,jsd)
    fndet_btm = fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (id_fndet_btm .gt. 0)           &
         used = g_send_data(id_fndet_btm,fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',fpdet_btm,isd,jsd)
    fpdet_btm = fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (id_fpdet_btm .gt. 0)           &
         used = g_send_data(id_fpdet_btm,fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',fsidet_btm,isd,jsd)
    fsidet_btm = fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (id_fsidet_btm .gt. 0)           &
         used = g_send_data(id_fsidet_btm,    fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine generic_COBALT_update_from_bottom_simple_slab

  ! <SUBROUTINE NAME="CBED_update_from_bottom">
  ! 
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !
  !   This routine calculates bottom fluxes for tracers with bottom reservoirs.
  !   It is called near the end of the time step, meaning that the fluxes 
  !   calculated pertain to the next time step.  
  ! 
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call CBED_update_from_bottom(cobalt,tracer_list,dt, tau, model_time) 
  !  </TEMPLATE>
  !
  !  <IN NAME="cobalt" TYPE="type(generic_COBALT_type)">
  !  </IN>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment 
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  !
  ! </SUBROUTINE>
  subroutine CBED_update_from_bottom(cobalt, tracer_list, dt, tau, model_time, fcadet_arag_btm, id_fcadet_arag_btm, fcadet_calc_btm, id_fcadet_calc_btm, ffedet_btm, id_ffedet_btm, flithdet_btm, id_flithdet_btm, fndet_btm, id_fndet_btm, fpdet_btm, id_fpdet_btm, fsidet_btm, id_fsidet_btm)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time
    real, dimension(:,:), intent(inout):: fcadet_arag_btm
    integer,              intent(in) :: id_fcadet_arag_btm
    real, dimension(:,:), intent(inout) :: fcadet_calc_btm
    integer,              intent(in) :: id_fcadet_calc_btm
    real, dimension(:,:), intent(inout) :: ffedet_btm
    integer,              intent(in) :: id_ffedet_btm
    real, dimension(:,:), intent(inout) :: flithdet_bt
    integer,              intent(in) :: id_flithdet_btm
    real, dimension(:,:), intent(inout) :: fndet_btm
    integer,              intent(in) :: id_fndet_btm
    real, dimension(:,:), intent(inout) :: fpdet_btm
    integer,              intent(in) :: id_fpdet_btm
    real, dimension(:,:), intent(inout) :: fsidet_btm
    integer,              intent(in) :: id_fsidet_btm

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    !
    ! The bottom reservoirs of aragonite and calcite are immediately redistributed to the
    ! water column as a bottom flux (btf) where they impact the alkalinity and DIC
    !
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',fcadet_arag_btm,isd,jsd)
    fcadet_arag_btm = fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(id_fcadet_arag_btm,fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',fcadet_calc_btm,isd,jsd)
    fcadet_calc_btm = fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(id_fcadet_calc_btm, fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',ffedet_btm,isd,jsd)
    ffedet_btm = ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (id_ffedet_btm .gt. 0)           &
         used = g_send_data(id_ffedet_btm, ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',flithdet_btm,isd,jsd)
    flithdet_btm = flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (id_flithdet_btm .gt. 0)           &
         used = g_send_data(id_flithdet_btm, flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of 
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',fndet_btm,isd,jsd)
    fndet_btm = fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (id_fndet_btm .gt. 0)           &
         used = g_send_data(id_fndet_btm,fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',fpdet_btm,isd,jsd)
    fpdet_btm = fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (id_fpdet_btm .gt. 0)           &
         used = g_send_data(id_fpdet_btm,fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',fsidet_btm,isd,jsd)
    fsidet_btm = fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (id_fsidet_btm .gt. 0)           &
         used = g_send_data(id_fsidet_btm,    fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine CBED_update_from_bottom

end module generic_COBALT_update_from_bottom_mod
