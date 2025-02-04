module generic_CBED

  use time_manager_mod,  only: time_type

  use g_tracer_utils, only : g_tracer_type, g_tracer_get_common
  use g_tracer_utils, only : g_tracer_set_values, g_tracer_get_values
  use g_tracer_utils, only : g_tracer_get_pointer, g_send_data

implicit none ; private

public generic_CBED_sediments_update_from_bottom
public generic_CBED_sediments_update_from_source

real, parameter :: sperd = 24.0 * 3600.0
real, parameter :: spery = 365.25 * sperd
real, parameter :: epsln=1.0e-30

contains

  ! subroutine generic_CBED_sediments_update_from_bottom is intended to be modified to become the CBED bottom layer model
  subroutine generic_CBED_sediments_update_from_bottom(fcadet_arag_btm, id_fcadet_arag_btm, &
           fcadet_calc_btm, id_fcadet_calc_btm, ffedet_btm, id_ffedet_btm, &
           flithdet_btm, id_flithdet_btm, fndet_btm, id_fndet_btm, &
           fpdet_btm, id_fpdet_btm, fsidet_btm, id_fsidet_btm, &
		   tracer_list, dt, tau, model_time)
    real, dimension(:,:), intent(inout) :: fcadet_arag_btm, fcadet_calc_btm, ffedet_btm, &
                                           flithdet_btm, fndet_btm, fpdet_btm, fsidet_btm
	integer, intent(inout)              :: id_fcadet_arag_btm, id_fcadet_calc_btm, id_ffedet_btm, &
                                           id_flithdet_btm, id_fndet_btm, id_fpdet_btm, id_fsidet_btm
    type(g_tracer_type), pointer        :: tracer_list
    real,               intent(in)      :: dt
    integer,            intent(in)      :: tau
    type(time_type),    intent(in)      :: model_time

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

  end subroutine generic_CBED_sediments_update_from_bottom

  ! Subroutine generic_CBED_sediments_update_from_source is intended to be modified to become the CBED bottom layer model
  subroutine generic_CBED_sediments_update_from_source(jfe_coast, fe_coast, f_ndet_btf, &
        c_2_n, frac_burial, zt, z_burial, fndet_burial, fpdet_burial, &
        f_pdet_btf, fno3denit_sed, f_no3, Rho_0, n_2_n_denit, &
        k_no3_denit, f_o2, o2_min, fnoxic_sed, o2_2_nh4, fnfeso4red_sed, &
        ffe_sed, fe_2_n_sed, ffe_sed_max, ffe_geotherm, ffe_iceberg, &
        ffe_iceberg_ratio, jprod_fed, fcased_redis_surfresp, f_cadet_calc_btf, &
        phi_surfresp_cased, cased_redis_coef, gamma_cased*max, omega_calc, &
        phi_deepresp_cased, alpha_cased, cased_redis_delz, f_lithdet_btf, &
        beta_cased, fcased_redis, f_cased, cased_steady, fcased_burial, &
        Co_cased, z_se, b_alk, f_cadet_arag_btf, alk_2_n_denit, b_dic, &
        b_fed, f_fedet_btf, b_nh4, b_no3, b_o2, b_po4, b_sio4, &
        f_sidet_btf, mask_coast, grid_tmask, grid_dat, grid_kmt, isc,iec, jsc,jec, &
        isd, jsd, nk, r_dt, internal_heat, dt, frunoff, rho_dzt, id_clock_ballast_loops)
    real, intent(inout) :: jfe_coast, fe_coast, f_ndet_btf, &
        c_2_n, frac_burial, zt, z_burial, fndet_burial, fpdet_burial, &
        f_pdet_btf, fno3denit_sed, f_no3, Rho_0, n_2_n_denit, &
        k_no3_denit, f_o2, o2_min, fnoxic_sed, o2_2_nh4, fnfeso4red_sed, &
        ffe_sed, fe_2_n_sed, ffe_sed_max, ffe_geotherm, ffe_iceberg, &
        ffe_iceberg_ratio, jprod_fed, fcased_redis_surfresp, f_cadet_calc_btf, &
        phi_surfresp_cased, cased_redis_coef, gamma_cased*max, omega_calc, &
        phi_deepresp_cased, alpha_cased, cased_redis_delz, f_lithdet_btf, &
        beta_cased, fcased_redis, f_cased, cased_steady, fcased_burial, &
        Co_cased, z_se, b_alk, f_cadet_arag_btf, alk_2_n_denit, b_dic, &
        b_fed, f_fedet_btf, b_nh4, b_no3, b_o2, b_po4, b_sio4, &
        f_sidet_btf
	real, dimension(ilb:,jlb:), intent(in) :: grid_dat
    real, dimension(:,:,:), intent(in)     :: grid_tmask
    integer, dimension(:,:), intent(in)    :: mask_coast, grid_kmt
    integer, intent(in)                    :: isc,iec, jsc,jec, isd, jsd, nk
    real, intent(in)                       :: r_dt, dt
    real, dimension(:,:), intent(in)       :: internal_heat
    real, dimension(:,:), intent(in)       :: frunoff
    real, dimension(:,:,:), intent(in)     :: rho_dzt
    integer, intent(in)                    :: id_clock_ballast_loops

    integer :: i, j, k
    real :: fpoc_btm, log_fpoc_btm
	
      !
    ! Coastal iron input (default is 0)
    !
    !do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
    !   jfe_coast(i,j,k) = fe_coast * mask_coast(i,j) * grid_tmask(i,j,k) / &
    !        sqrt(grid_dat(i,j))
    !     ! uncomment if running "no mass change" test
    !     !jfe_coast(i,j,k) = 0.0
    !enddo; enddo; enddo  !} i,j,k

    do j = jsc, jec; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{
          !
          ! Nitrogen flux from the sediments
          ! 
          if (f_ndet_btf(i,j,1) .gt. 0.0) then !{
             ! fpoc_bottom in mmoles C m-2 day-1 for burial relationship
             fpoc_btm = f_ndet_btf(i,j,1)*c_2_n*sperd*1000.0
             !frac_burial(i,j) = (0.013 + 0.53*fpoc_btm**2.0)/((7.0+fpoc_btm)**2.0)
             frac_burial(i,j) = 0.013 + 0.53*fpoc_btm**2.0/((7.0+fpoc_btm)**2.0) * &
                  zt(i,j,k) / (z_burial + zt(i,j,k))
             ! uncomment for "no mass change" test
             !frac_burial(i,j) = 0.0
             fndet_burial(i,j) = frac_burial(i,j)*f_ndet_btf(i,j,1)
             fpdet_burial(i,j) = frac_burial(i,j)*f_pdet_btf(i,j,1)
             ! fpoc_bottom in micromoles C cm-2 day-1 for denitrification relationship, cap at 43
             ! to prevent anomalous extrapolation of the relationship
             log_fpoc_btm = log(min(43.0,0.1*fpoc_btm))
             fno3denit_sed(i,j) = min(f_no3(i,j,k)*Rho_0*r_dt,  &      
                  min((f_ndet_btf(i,j,1)-fndet_burial(i,j))*n_2_n_denit, & 
                  10.0**(-0.9543+0.7662*log_fpoc_btm - 0.235*log_fpoc_btm**2.0)/(c_2_n*sperd*100.0)* &
                  n_2_n_denit*f_no3(i,j,k)/(k_no3_denit + f_no3(i,j,k)))) * &
                  zt(i,j,k) / (z_burial + zt(i,j,k))
             ! uncomment "no mass change" test 
             !fno3denit_sed(i,j) = 0.0             
             if (f_o2(i,j,k) .gt. o2_min) then  !{
                fnoxic_sed(i,j) = max(0.0, min(f_o2(i,j,k)*Rho_0*r_dt*(1.0/o2_2_nh4), &
                                         f_ndet_btf(i,j,1) - fndet_burial(i,j) - &
                                         fno3denit_sed(i,j)/n_2_n_denit))
             else
                fnoxic_sed(i,j) = 0.0
             endif !}
             fno3denit_sed(i,j) = fno3denit_sed(i,j) + &
                                         min(f_no3(i,j,k)*Rho_0*r_dt-fno3denit_sed(i,j), &
                                         (f_ndet_btf(i,j,1)-fnoxic_sed(i,j)-fndet_burial(i,j) - &
                                         fno3denit_sed(i,j)/n_2_n_denit)*n_2_n_denit)
             fnfeso4red_sed(i,j) = max(0.0, f_ndet_btf(i,j,1)-fnoxic_sed(i,j)- &
                                          fndet_burial(i,j)-fno3denit_sed(i,j)/n_2_n_denit)
          else
             fnfeso4red_sed(i,j) = 0.0
             fno3denit_sed(i,j) = 0.0
             fnoxic_sed(i,j) = 0.0
          endif !}

          ! iron from sediment (Elrod) 
          !ffe_sed(i,j) = fe_2_n_sed * f_ndet_btf(i,j,1)
          ! iron from sediment (Dale)
          ffe_sed(i,j) = ffe_sed_max * tanh( (f_ndet_btf(i,j,1)*c_2_n*sperd*1.0e3)/ &
                                max(f_o2(i,j,k)*1.0e6,epsln) )

          ffe_geotherm(i,j) = ffe_geotherm_ratio*internal_heat(i,j)*4184.0/dt
          ! default for icebergs: 40 nanomoles fe dissolved per kg of icemelt
          ! sediments: Raiswell et al., 2008: 0.5 kg sed per m-3 of iceberg; 0.1% mean Fe, 5-10% soluble
          ! ~500 nanomoles Fe per kg-1 icemelt 
          ffe_iceberg(i,j) = ffe_iceberg_ratio*max(frunoff(i,j),0.0)
          jprod_fed(i,j,1) = jprod_fed(i,j,1) + ffe_iceberg(i,j)/rho_dzt(i,j,1) 

          !
          ! Calcium carbonate flux and burial
          ! 2015/11/18 JGJ: fix from JPD to cap the absolute cased dissolution rate to 10 mmol m-2 d-1
          ! Calcite cycling in the sediments is based on the model of   Dunne et al., 2012.
          !
          ! phi_surfresp_cased = 0.14307   ! const for enhanced diss., surf sed respiration (dimensionless)
          ! phi_deepresp_cased = 4.1228    ! const for enhanced diss., deep sed respiration (dimensionless)
          ! alpha_cased = 2.7488 ! exponent controlling non-linearity of deep dissolution           
          ! beta_cased = -2.2185 ! exponent controlling non-linearity of effective thickness
          ! gamma_cased = 0.03607/spery   ! dissolution rate constant
          ! Co_cased = 8.1e3        ! moles CaCo3 m-3 for pure calcite sediment with porosity = 0.7
          !
          ! if cased_steady is true, burial is calculated from Dunne's eq. (2) assuming dcased/dt = 0.
          ! This ensures that all the calcite bottom flux is partitioned between burial and redissolution.     
          ! The steady state cased value of cased is calculated to reflect the changing bottom conditions.
          ! This influences the the partitioning of burial and redissolution over time, but there are
          ! no alkalinity changes/drifts associated with the long-term evolution of cased
          ! 
          ! If cased_steady is false, calcite is partitioned between dissolution, burial and evolving
          ! cased as described in Dunne et al. (2012).  The multi-century scale evolution of cased
          ! impacts alkalinity, but care must to ensure that cased starts in equilibrium with the 
          ! mean ocean state to avoid unrealistic drifts. 

          ! Enhanced dissolution by fast respiration near the sediment surface, proportional 
          ! to organic flux, moles Ca m-2 s-1, limited to a max 1/2 the instantaneous calcite flux
          fcased_redis_surfresp(i,j)=min(0.5*f_cadet_calc_btf(i,j,1), &
            phi_surfresp_cased*f_ndet_btf(i,j,1)*c_2_n)
          ! Ca-specific dissolution coeficient, depends on calcite saturation state and is enhanced by
          ! respiration deep in the sediment (s-1), non-linearity controlled by alpha_cased
          cased_redis_coef(i,j) = gamma_cased*max(0.0,1.0-omega_calc(i,j,k)+ &
            phi_deepresp_cased*f_ndet_btf(i,j,1)*c_2_n*spery)**alpha_cased
          ! Effective thickness term that enhances burial of calcite when total sediment accumulation is high
          ! dimensionless value between 0 and 1
          cased_redis_delz(i,j) = max(1.0, &
            f_lithdet_btf(i,j,1)*spery+f_cadet_calc_btf(i,j,1)*100.0*spery)**beta_cased  
          ! calculate the sediment redissolution rate (moles Ca m-2 sec-1). This calculation is subject to
          ! three limiters: a) a maximum of 1/2 of the total cased over one time step; b) a maximum of 0.01
          ! moles Ca per day; and c) a minimum of 0.0
          fcased_redis(i,j) = max(0.0, min(0.01/sperd, min(0.5*f_cased(i,j,1)*r_dt,  &
            fcased_redis_surfresp(i,j)+cased_redis_coef(i,j)*cased_redis_delz(i,j)*f_cased(i,j,1))) ) 
          !
          ! Old expression
          !
          !fcased_redis(i,j) = max(0.0, min(0.01/sperd,min(0.5 * f_cased(i,j,1) * r_dt, min(0.5 *       &                          
          !   f_cadet_calc_btf(i,j,1), 0.14307 * f_ndet_btf(i,j,1) * c_2_n) +        &
          !   0.03607 / spery * max(0.0, 1.0 - omega_calc(i,j,k) +   &
          !   4.1228 * f_ndet_btf(i,j,1) * c_2_n * spery)**(2.7488) *                        &
          !   max(1.0, f_lithdet_btf(i,j,1) * spery + f_cadet_calc_btf(i,j,1) * 100.0 *  &
          !   spery)**(-2.2185) * f_cased(i,j,1))))*grid_tmask(i,j,k)

          if (cased_steady) then
            fcased_burial(i,j) = f_cadet_calc_btf(i,j,1) - fcased_redis(i,j)
            f_cased(i,j,1) = fcased_burial(i,j)*Co_cased/f_cadet_calc_btf(i,j,1)
          else
            fcased_burial(i,j) = max(0.0, f_cadet_calc_btf(i,j,1) * f_cased(i,j,1) / &
              Co_cased)
            f_cased(i,j,1) = f_cased(i,j,1) + (f_cadet_calc_btf(i,j,1) -            &
              fcased_redis(i,j) - fcased_burial(i,j)) / z_sed * dt *                &
              grid_tmask(i,j,k)
          endif

          ! uncomment for "no mass change" test (next 3 lines)
          !fcased_redis(i,j) = f_cadet_calc_btf(i,j,1)
          !fcased_burial(i,j) = 0.0
          !f_cased(i,j,1) = f_cased(i,j,1)
          !
          ! Bottom flux boundaries passed to the vertical mixing routine 
          !

          b_alk(i,j) = - 2.0*(fcased_redis(i,j)+f_cadet_arag_btf(i,j,1)) -    &
             fnoxic_sed(i,j) - fno3denit_sed(i,j)*alk_2_n_denit
          b_dic(i,j) =  - fcased_redis(i,j) - f_cadet_arag_btf(i,j,1) -       &
             (f_ndet_btf(i,j,1) - fndet_burial(i,j)) * c_2_n
          ! uncomment for "no mass change" test (next 2 lines)
          !b_dic(i,j) =  - f_cadet_calc_btf(i,j,1)  - f_cadet_arag_btf(i,j,1) -            &
          !   (f_ndet_btf(i,j,1) - fndet_burial(i,j)) * c_2_n 
          b_fed(i,j) = - ffe_sed(i,j) - ffe_geotherm(i,j)
          ! uncomment for "no mass change" test (next line)
          !b_fed(i,j) = - f_fedet_btf(i,j,1)
          b_nh4(i,j) = - f_ndet_btf(i,j,1) + fndet_burial(i,j)
          b_no3(i,j) = fno3denit_sed(i,j)
          b_o2(i,j)  = o2_2_nh4 * (fnoxic_sed(i,j) + fnfeso4red_sed(i,j))
          b_po4(i,j) = - f_pdet_btf(i,j,1) + fpdet_burial(i,j)
          b_sio4(i,j)= - f_sidet_btf(i,j,1)

       endif !}
    enddo; enddo  !} i, j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       f_cased(i,j,k) = 0.0
    enddo; enddo ; enddo  !} i,j,k

    call mpp_clock_end(id_clock_ballast_loops)

    call g_tracer_set_values(tracer_list,'alk',  'btf', b_alk ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic',  'btf', b_dic ,isd,jsd)
    call g_tracer_set_values(tracer_list,'fed',  'btf', b_fed ,isd,jsd)
    call g_tracer_set_values(tracer_list,'nh4',  'btf', b_nh4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'no3',  'btf', b_no3 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2',   'btf', b_o2  ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4',  'btf', b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'sio4', 'btf', b_sio4,isd,jsd)
!
    call mpp_clock_begin(id_clock_source_sink_loop1)
!

  end subroutine generic_CBED_sediments_update_from_source

end module generic_CBED