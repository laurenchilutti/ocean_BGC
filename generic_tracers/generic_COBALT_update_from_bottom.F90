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
  subroutine generic_COBALT_update_from_bottom_simple_slab(cobalt, tracer_list, dt, tau, model_time)
    type(generic_COBALT_type), intent(inout) :: cobalt
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
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',cobalt%fcadet_arag_btm,isd,jsd)
    cobalt%fcadet_arag_btm = cobalt%fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (cobalt%id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_arag_btm,cobalt%fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',cobalt%fcadet_calc_btm,isd,jsd)
    cobalt%fcadet_calc_btm = cobalt%fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (cobalt%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_calc_btm, cobalt%fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',cobalt%ffedet_btm,isd,jsd)
    cobalt%ffedet_btm = cobalt%ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = cobalt%ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (cobalt%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_ffedet_btm, cobalt%ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',cobalt%flithdet_btm,isd,jsd)
    cobalt%flithdet_btm = cobalt%flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (cobalt%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_flithdet_btm, cobalt%flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of 
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',cobalt%fndet_btm,isd,jsd)
    cobalt%fndet_btm = cobalt%fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (cobalt%id_fndet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fndet_btm,cobalt%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',cobalt%fpdet_btm,isd,jsd)
    cobalt%fpdet_btm = cobalt%fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (cobalt%id_fpdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fpdet_btm,cobalt%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',cobalt%fsidet_btm,isd,jsd)
    cobalt%fsidet_btm = cobalt%fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (cobalt%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fsidet_btm,    cobalt%fsidet_btm,          &
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
  subroutine CBED_update_from_bottom(cobalt, tracer_list, dt, tau, model_time)
    type(generic_COBALT_type), intent(inout) :: cobalt
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
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',cobalt%fcadet_arag_btm,isd,jsd)
    cobalt%fcadet_arag_btm = cobalt%fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (cobalt%id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_arag_btm,cobalt%fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',cobalt%fcadet_calc_btm,isd,jsd)
    cobalt%fcadet_calc_btm = cobalt%fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (cobalt%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_calc_btm, cobalt%fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',cobalt%ffedet_btm,isd,jsd)
    cobalt%ffedet_btm = cobalt%ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = cobalt%ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (cobalt%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_ffedet_btm, cobalt%ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',cobalt%flithdet_btm,isd,jsd)
    cobalt%flithdet_btm = cobalt%flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (cobalt%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_flithdet_btm, cobalt%flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of 
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',cobalt%fndet_btm,isd,jsd)
    cobalt%fndet_btm = cobalt%fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (cobalt%id_fndet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fndet_btm,cobalt%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',cobalt%fpdet_btm,isd,jsd)
    cobalt%fpdet_btm = cobalt%fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (cobalt%id_fpdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fpdet_btm,cobalt%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',cobalt%fsidet_btm,isd,jsd)
    cobalt%fsidet_btm = cobalt%fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (cobalt%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fsidet_btm,    cobalt%fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine CBED_updiate_from_bottom

end module generic_COBALT_update_from_bottom_mod
