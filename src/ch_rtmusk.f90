      subroutine ch_rtmusk
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes a daily flow through a reach using the
!!    Muskingum method

!!    code provided by Dr. Valentina Krysanova, Pottsdam Institute for
!!    Climate Impact Research, Germany
!!    Modified by Balaji Narasimhan
!!    Spatial Sciences Laboratory, Texas A&M University

      use basin_module
      use channel_data_module
      use channel_module
      use hydrograph_module !, only : ob, icmd, jrch, isdch, fp_stor, ch_stor, wet
      use time_module
      use channel_velocity_module
      use sd_channel_module
      use sd_channel_icejam_module
      use climate_module
      use reservoir_module
      use reservoir_data_module
      use water_body_module
      use conditional_module
      use gwflow_module
      
      implicit none
      
      external :: rcurv_interp_flo
      real, external :: qman, theta
      
      integer :: ii = 0     !none              |current day of simulation
      integer :: ihru = 0
      integer :: iihru = 0
      integer :: icha = 0
      integer :: irtstep = 0
      integer :: isubstep = 0
      
      real :: ch_stor_init = 0. !m3             |storage in channel at beginning of day
      real :: fp_stor_init = 0. !m3             |storage in flood plain above wetlands emergency spillway at beginning of day
      real :: wet_stor_init = 0.  !m3           |storage in flood plain wetlands at beginning of day
      real :: tot_stor_init = 0.
      real :: inout = 0.        !m3             |inflow - outflow for day
      real :: del_stor = 0.     !m3             |change in storage of channel + flood plain + wetlands
      real :: topw = 0.         !m              |top width of main channel
      real :: qinday = 0.       !units          |description 
      real :: qoutday = 0.      !units          |description   
      real :: inflo = 0.        !m^3            |inflow water volume
      real :: inflo_rate = 0.   !m^3/s          |inflow rate
      real :: outflo = 0.       !m^3            |outflow water volume
      real :: trans_loss = 0.   !m^3            |transmission losses during day
      real :: evap = 0.         !m^3            |evaporation losses during day
      real :: rto = 0.
      real :: outflo_rate = 0.
      real :: dts = 0.               !seconds    |time step interval for substep
      real :: dthr = 0.
      real :: scoef = 0.
      real :: sum_inflo = 0.
      real :: sum_outflo = 0.
      real :: wet_evol = 0. 
      real :: ratio = 0.

      ! Icejam dynamic Muskingum override variables.  Wedge water-balance
      ! capture/release is owned by sd_channel_icejam; this routine only routes
      ! the ice-adjusted inflow and applies phase-dependent K/X modifiers.
      type(icejam_param_type), save :: ice_prm
      logical, save :: ice_prm_initialized = .false.
      real :: k_norm_hr = 0.
      real :: k_cur_hr = 0.
      real :: x_cur = 0.
      real :: denom_msk = 0.
      real :: c1_ice = 0.
      real :: c2_ice = 0.
      real :: c3_ice = 0.
      real :: k_lower = 0.

      if (.not. ice_prm_initialized) then
        call icejam_default_params(ice_prm)
        call icejam_validate_params(ice_prm)
        ice_prm_initialized = .true.
      endif

      jrch = isdch
      jhyd = sd_dat(jrch)%hyd
      
      qinday = 0
      qoutday = 0
      ht2 = hz
      ob(icmd)%hyd_flo = 0.
      hyd_rad = 0.
      trav_time = 0.
      flo_dep = 0.
      trans_loss = 0.
      evap = 0.
      ch_wat_d(jrch)%evap = 0.
      ch_wat_d(jrch)%seep = 0.
      
      sum_inflo = sum (ob(icmd)%tsin)
        
      !! total wetland volume at start of day
      wet_stor(jrch) = hz
      wet_evol = 0.
      do ihru = 1, sd_ch(jrch)%fp%hru_tot
        iihru = sd_ch(jrch)%fp%hru(ihru)
        wet_stor(jrch) = wet_stor(jrch) + wet(iihru)
        wet_evol = wet_evol + wet_ob(iihru)%evol
      end do
      wet_stor_init = wet_stor(jrch)%flo
      ch_stor_init = ch_stor(jrch)%flo
      fp_stor_init = fp_stor(jrch)%flo
      tot_stor(jrch) = ch_stor(jrch) + fp_stor(jrch)
      tot_stor_init = tot_stor(jrch)%flo

      !! Icejam module has already adjusted ht1%flo for wedge capture/release.
      !! ch_rtmusk remains the routing engine.  Optional conduit inflow is added
      !! here as an additional channel inflow and is not captured by icejam in
      !! this first clean-boundary implementation.

      !add groundwater conduit if applied
      if (bsn_cc%gwflow > 0 .and. gw_conduit_flag > 0) then
          ht1 = ht1 + gw_conduit_info(jrch)%output
          !ratio = gw_conduit_info(jrch)%output%flo / ht1%flo
          !if (jrch == 129) write (9003,*) "rtmusk, conduit inflow:", gw_conduit_info(jrch)%output%flo, "ratio:", ratio
      endif

      ! Channel water balance sees the current inflow hydrograph.  Any icejam
      ! wedge capture/release has been applied upstream in sd_channel_icejam.
      sum_inflo = ht1%flo

      !! keep Muskingum substeps computed in sd_hydsed_init
      !! For daily simulations, substeps may be greater than 1 to satisfy
      !! the Muskingum stability criterion. Do not reset them here.
      !! set for daily time step
      !if (time%step == 1) then
      !  sd_ch(jrch)%msk%nsteps = 1
      !  sd_ch(jrch)%msk%substeps = 1
      !end if
      irtstep = 1
      isubstep = 0
      dts = time%dtm / sd_ch(jrch)%msk%substeps * 60.
      dthr = dts / 3600.
      
      !! subdaily time step
      do ii = 1, sd_ch(jrch)%msk%nsteps
        !! water entering reach during time step - substeps for stability
        isubstep = isubstep + 1
        if (isubstep > sd_ch(jrch)%msk%substeps) then
          irtstep = irtstep + 1
          isubstep = 1
        end if
        
        !! inflow for the current Muskingum substep
        !! ob(icmd)%tsin(irtstep) is the inflow volume for the routing step;
        !! divide it by substeps to obtain the volume entering this substep.
        !!inflo = ob(icmd)%tsin(irtstep) / sd_ch(jrch)%msk%substeps
        
        !! add inflow and associated constituents to total storage
        if (ht1%flo > 1.e-6) then
          inflo = ht1%flo / sd_ch(jrch)%msk%substeps
          
          rto = inflo / ht1%flo
          rto = Max(0., rto)
          rto = Min(1., rto)
          tot_stor(jrch) = tot_stor(jrch) + rto * ht1
        else
          inflo = 0.
        end if    ! ht1%flo > 1.e-6
        
        !! interpolate rating curve using inflow rate for this substep
        icha = jrch
        inflo_rate = inflo / dts
        call rcurv_interp_flo (icha, inflo_rate)
        ch_rcurv(jrch)%in2 = rcurv
        
        !! if no water in channel - skip routing and set rating curves to zero
        if (tot_stor(jrch)%flo < 1.e-6) then
          ch_rcurv(jrch)%in1 = rcz
          ch_rcurv(jrch)%out1 = rcz
          sd_ch(jrch)%in1_vol = 0.
          sd_ch(jrch)%out1_vol = 0.
        else
          if (bsn_cc%rte == 1) then
          !! Muskingum flood routing method. Icejam optionally overrides
          !! K and X to represent ice-cover/jam hydraulic resistance only.
            if (bsn_cc%icejam == 1 .and. sd_ch(jrch)%ice_hydro_active == 1) then
              k_norm_hr = max(1.0e-3, ch_rcurv(jrch)%in2%ttime)
              if (sd_ch(jrch)%stor_dis_bf > 1.0e-6) k_norm_hr = max(1.0e-3, sd_ch(jrch)%stor_dis_bf)

              k_cur_hr = k_norm_hr * max(ice_prm%k_min_mult, min(ice_prm%k_max_mult, sd_ch(jrch)%ice_k_mult))
              x_cur = max(0.0, min(0.49, sd_ch(jrch)%ice_x_current))
              ! Active ice jam is reservoir-like; enforce tiny X even if an
              ! upstream state assignment fails to pass x_jam correctly.
              if (sd_ch(jrch)%ice_phase == 3 .and. sd_ch(jrch)%is_jamming) then
                x_cur = min(x_cur, ice_prm%x_jam)
              endif

              !! Stability guards for 2 K X <= dt <= 2 K (1-X).
              if (2.0 * k_cur_hr * x_cur > dthr) then
                x_cur = min(x_cur, 0.49 * dthr / max(k_cur_hr, 1.0e-6))
              endif
              k_lower = dthr / max(2.0 * (1.0 - x_cur), 1.0e-6)
              if (k_cur_hr < k_lower) k_cur_hr = 1.001 * k_lower

              denom_msk = 2.0 * k_cur_hr * (1.0 - x_cur) + dthr
              c1_ice = (dthr - 2.0 * k_cur_hr * x_cur) / denom_msk
              c2_ice = (dthr + 2.0 * k_cur_hr * x_cur) / denom_msk
              c3_ice = (2.0 * k_cur_hr * (1.0 - x_cur) - dthr) / denom_msk

              outflo = c1_ice * inflo + c2_ice * sd_ch(jrch)%in1_vol + c3_ice * sd_ch(jrch)%out1_vol
            else
              outflo = sd_ch(jrch)%msk%c1 * inflo + sd_ch(jrch)%msk%c2 * sd_ch(jrch)%in1_vol +     &
                                                  sd_ch(jrch)%msk%c3 * sd_ch(jrch)%out1_vol
            endif
            outflo = Min (outflo, tot_stor(jrch)%flo)
            outflo = Max (outflo, 0.)
               
            !! save inflow/outflow volumes for next time step (and day) for Muskingum
            sd_ch(jrch)%in1_vol = inflo
            sd_ch(jrch)%out1_vol = outflo
          else

            !! Variable Storage Coefficient method - sc=2*dt/(2*ttime+dt) - ttime=(in2+out1)/2
            scoef = dthr / (ch_rcurv(jrch)%in2%ttime + ch_rcurv(jrch)%out1%ttime + dthr)
            scoef = bsn_prm%scoef * 2. * dthr / (2.* ch_rcurv(jrch)%out1%ttime + dthr)
            scoef = Min (scoef, 1.)
            outflo = scoef * tot_stor(jrch)%flo
          end if
          
          !! compute outflow rating curve for next time step
          outflo_rate = outflo / dts      !convert to cms
          call rcurv_interp_flo (jrch, outflo_rate)
          ch_rcurv(jrch)%out2 = rcurv
 
          !! add outflow to daily hydrograph and subdaily flow
          rto = outflo / tot_stor(jrch)%flo
          rto = Min (1., rto)
          ht2 = ht2 + rto * tot_stor(jrch)
          ob(icmd)%hyd_flo(1,irtstep) = ob(icmd)%hyd_flo(1,irtstep) + outflo
          !! subtract outflow from total storage
          tot_stor(jrch) = (1. - rto) * tot_stor(jrch)
        
          !! set rating curve for next time step
          ch_rcurv(jrch)%in1 = ch_rcurv(jrch)%in2
          ch_rcurv(jrch)%out1 = ch_rcurv(jrch)%out2
          
          !! partition channel and flood plain based on bankfull volume
          if (tot_stor(jrch)%flo > ch_rcurv(jrch)%elev(2)%vol_ch) then
            !! fill channel to bank full if below
            rto = (tot_stor(jrch)%flo - ch_rcurv(jrch)%elev(2)%vol_ch) / tot_stor(jrch)%flo
            fp_stor(jrch) = rto * tot_stor(jrch)
            ch_stor(jrch) = (1. - rto) * tot_stor(jrch)
          else
            ch_stor(jrch) = tot_stor(jrch)
            fp_stor(jrch) = hz
          end if
        
          tot_stor(jrch) = ch_stor(jrch) + fp_stor(jrch)
          
        end if  ! tot_stor(jrch)%flo < 1.e-6

      end do    ! end of sub-daily loop

      !! =========================================================
      !! Icejam:
      !! Explicit ice-jam water capture/release is handled by sd_channel_icejam
      !! before routing. No additional ch_stor storage-busting release is
      !! applied here.
      !! =========================================================

      !! compute water balance - evap and seep
      !! calculate transmission losses (seepage), only when gwflow is not activated
      if(bsn_cc%gwflow == 0) then
        if (ch_stor(jrch)%flo > 1.e-6) then
          !! mm/hr * km * m * 24. = m3
          trans_loss = sd_ch(jrch)%chk * sd_ch(jrch)%chl * rcurv%wet_perim * 24.
          !trans_loss = sd_ch(jrch)%chk * sd_ch(jrch)%chl * sd_ch(jrch)%chw * 24.
          trans_loss = Min(trans_loss, ch_stor(jrch)%flo)
          !! subtract transmission loses from outflow
          rto = trans_loss / ch_stor(jrch)%flo
          ch_stor(jrch) = (1. - rto) * ch_stor(jrch)
        end if
        ch_wat_d(jrch)%seep = trans_loss
      else
        !if gwflow active, seepage computed in gwflow routine
        ch_wat_d(jrch)%seep = 0.
      endif

      !! calculate evaporation losses
      if (ch_stor(jrch)%flo > 1.e-6) then
        !! calculate width of channel at water level - flood plain evap calculated in wetlands
        !if (dep_flo <= sd_ch(jrch)%chd) then
        !  topw = ch_rcurv(jrch)%out2%surf_area
        !else
          topw = 1000. * sd_ch(jrch)%chl * sd_ch(jrch)%chw
        !end if
        iwst = ob(icmd)%wst
        !! mm/day * m2 / 1000.
        evap = bsn_prm%evrch * wst(iwst)%weat%pet * topw / 1000.
        evap = Min(evap, ch_stor(jrch)%flo)
        rto = evap / ch_stor(jrch)%flo
        ch_stor(jrch)%flo = (1. - rto) * ch_stor(jrch)%flo
      end if
      ch_wat_d(jrch)%evap = evap
      
      tot_stor(jrch) = ch_stor(jrch) + fp_stor(jrch)

      !! check water balance at end of day
      sum_outflo = ht2%flo
      inout = sum_inflo - sum_outflo - trans_loss - evap
      !! total wetland volume at end of day
      wet_stor(jrch) = hz
      do ihru = 1, sd_ch(jrch)%fp%hru_tot
        iihru = sd_ch(jrch)%fp%hru(ihru)
        wet_stor(jrch) = wet_stor(jrch) + wet(iihru)
      end do
      del_stor = (ch_stor(jrch)%flo - ch_stor_init) + (fp_stor(jrch)%flo - fp_stor_init) +          &
                                                    (wet_stor(jrch)%flo - wet_stor_init)
      ch_fp_wb(jrch)%inflo = sum_inflo
      ch_fp_wb(jrch)%outflo = sum_outflo
      ch_fp_wb(jrch)%tl = trans_loss
      ch_fp_wb(jrch)%ev = evap
      ch_fp_wb(jrch)%ch_stor_init = ch_stor_init
      ch_fp_wb(jrch)%ch_stor = ch_stor(jrch)%flo
      ch_fp_wb(jrch)%fp_stor_init = fp_stor_init
      ch_fp_wb(jrch)%fp_stor = fp_stor(jrch)%flo
      ch_fp_wb(jrch)%tot_stor_init = tot_stor_init
      ch_fp_wb(jrch)%tot_stor = tot_stor(jrch)%flo
      ch_fp_wb(jrch)%wet_stor_init = wet_stor_init
      ch_fp_wb(jrch)%wet_stor = wet_stor(jrch)%flo

      !!conceptual ice-jam: transport downstream
      if(bsn_cc%icejam == 1) then
          call sd_channel_ice_advect(jrch)
      end if

!      if (jrch == 68) then
!          write(9003,*) time%yrc, time%day, jrch, &
!                  "phase", sd_ch(jrch)%ice_phase, &
!                  "jamming", sd_ch(jrch)%is_jamming, &
!                  "releasing", sd_ch(jrch)%is_releasing, &
!                  "Kmult", sd_ch(jrch)%ice_k_mult, &
!                  "Xice", sd_ch(jrch)%ice_x_current, &
!                  "ch_stor", ch_stor(jrch)%flo, &
!                  "fp_stor", fp_stor(jrch)%flo, &
!                  "tot_stor", tot_stor(jrch)%flo, &
!                  "bnk_stor", sd_ch(jrch)%bankfull_storage, &
!                  "wedge_stor", sd_ch(jrch)%ice_wedge_stor, &
!                  "wedge_cap", sd_ch(jrch)%ice_wedge_capacity, &
!                  "wedge_cap_day", sd_ch(jrch)%ice_wedge_capture, &
!                  "wedge_rel", sd_ch(jrch)%ice_wedge_release, &
!                  "wedge_leak", sd_ch(jrch)%ice_wedge_leak, &
!                  "wedge_ratio", sd_ch(jrch)%ice_wedge_stor / max(sd_ch(jrch)%ice_wedge_capacity, 1.e-6), &
!                  "stor_ratio", ch_stor(jrch)%flo / max(sd_ch(jrch)%bankfull_storage, 1.e-6), &
!                  "excess", sd_ch(jrch)%ice_excess_storage, &
!                  "shock_rel", sd_ch(jrch)%ice_shock_release, &
!                  "Khr", k_cur_hr, &
!                  "Xcur", x_cur, &
!                  "c1", c1_ice, &
!                  "c2", c2_ice, &
!                  "c3", c3_ice, &
!                  "force", sd_ch(jrch)%force_eff, &
!                  "resist", sd_ch(jrch)%resistance, &
!                  "F_R", sd_ch(jrch)%force_eff / max(sd_ch(jrch)%resistance, 1.e-6)
!      endif

      return
      end subroutine ch_rtmusk