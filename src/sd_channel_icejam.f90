subroutine sd_channel_icejam(j)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Icejam: strictly one-way seasonal phase plus repeatable jam-event module.
!!
!!    Scientific abstraction:
!!      * Solid channel ice is represented by sd_ch%ice_vol, with structural
!!        integrity sd_ch%ice_integrity.
!!      * Liquid-water backwater/impoundment caused by ice is represented by one
!!        conceptual wedge storage sd_ch%ice_wedge_stor.
!!      * Mobile/broken ice is represented by sd_ch%ice_mobile and
!!        sd_ch%ice_mobile_pass; it is ice mass, not liquid water.
!!      * Seasonal phase is a slow one-way background regime.  Dynamic jam/release
!!        episodes are separate boolean states that may repeat within FREEZEUP or BREAKUP.
!!      * preserves the event chain: jam formation -> wedge impoundment ->
!!        jam-break release.  Daily wedge leakage is never reported as release.
!!
!!    Mass-balance boundary:
!!      * This routine modifies only ht1%flo, ch_stor%flo, and sd_ch ice/wedge
!!        states before ch_rtmusk.
!!      * Muskingum routing remains in ch_rtmusk.  This routine only passes
!!        phase-dependent K/X directives through sd_ch%ice_hydro_active,
!!        sd_ch%ice_k_mult, and sd_ch%ice_x_current.

      use basin_module
      use time_module
      use hydrograph_module
      use sd_channel_module
      use climate_module
      use sd_channel_icejam_module

      implicit none

      integer, intent(in) :: j

      integer, parameter :: ICE_OPEN     = 0
      integer, parameter :: ICE_FREEZEUP = 1
      integer, parameter :: ICE_STABLE   = 2
      integer, parameter :: ICE_BREAKUP  = 3

      type(icejam_param_type), save :: prm
      type(icejam_reach_scale_type) :: reach
      logical, save :: prm_initialized = .false.

      integer :: old_phase
      integer :: new_phase
      real :: t_air, tmax, tw_ice, t_grow, t_decay
      real :: freeze_drive, thaw_drive, t_thaw
      real :: thaw_weak
      real :: struct_gain, struct_loss, struct_loss_cap
      real :: surface_weak_drive, surface_recovery, surface_integrity
      real :: area_ha, ros_area_frac, ros_water_mm
      real :: ht1_raw, ht1_adj, q_in, q_prev
      real :: q_bankfull, q_bnk_fallback
      real :: ch_len_m, bankfull_vol
      real :: ice_area, sim_ice_thick, ice_maturity
      real :: deep_winter_factor, ice_storage_factor, snowpack_factor
      real :: frozen_soil_factor, warm_flush_factor, discharge_factor
      real :: major_jam_factor, major_release_cap
      real :: major_base_factor, major_background_factor, major_trigger_factor, warm_memory_factor
      real :: snowpack_ante_factor, snowpack_peak_factor
      real :: warm_air_factor, meltwater_factor
      real :: fr_factor, qrise_factor, runoff_response_factor
      real :: post_release_capture_eff
      logical :: major_release_active
      logical :: deep_winter_ready, channel_ice_ready, snowpack_ready
      logical :: frozen_soil_ready, warm_flush_ready, discharge_ready
      logical :: major_release_gate, major_storage_ready
      real :: ice_target_thick, ice_target_vol, grow_pot
      real :: grow_act, melt_thick, melt_pot, melt_act
      real :: freeze_avail, freeze_from_stor, freeze_from_in, freeze_remain
      real :: local_shock_q, thermal_shock, damage_factor, force_mult
      real :: force_F, resistance_R, fr_ratio, alpha_ice
      real :: mobile_factor, mobile_gen, mobile_thermal, mobile_mech
      real :: mobile_supply, mobile_remain
      real :: mobile_melt_pot, mobile_melt_act, pass_melt_pot, pass_melt_act
      real :: q_rel, trans_flow_factor, trans_susc_factor
      real :: ice_transport_cap, ice_mobile_excess, mobile_assim
      real :: wedge_capacity, wedge_avail, wedge_ratio, wedge_ratio_for_major
      real :: wedge_capture_frac, capture_mult, retention_frac
      real :: underice_alpha, underice_capacity, underice_excess
      real :: background_capture, excess_capture, wedge_capture
      real :: wedge_release_frac, wedge_release, wedge_leak
      real :: tail_factor, tail_material_factor, leak_mult
      real :: stable_storage_relief, stable_flow_boost, stable_capture_factor
      real :: stable_leak_mult, q_rise_pos, q_fall_pos
      real :: capture_mobile_term
      real :: deepwinter_cover_factor, deepwinter_flow_factor
      real :: deepwinter_cover_index, deep_underice_alpha
      real :: deep_retention_frac, deep_block_capacity, deep_wedge_capture
      real :: winter_alpha_intact, winter_alpha_weak
      real :: winter_qcap_intact, winter_qcap_weak, winter_capacity_leak
      real :: winter_extra_qcap, winter_actual_excess_q, winter_pulse_factor
      real :: winter_additional_leak
      real :: ordinary_release_cap
      real :: release_force_eff
      real :: q_damp_eff
      logical :: cold_window, warm_cleanup_allowed
      logical :: freezeup_ready, stable_ready
      logical :: breakup_material_ready, breakup_weather_ready
      logical :: breakup_age_ready, breakup_force_ready, breakup_ready
      logical :: jam_material_ready, jam_transport_ready, jam_storage_ready
      logical :: jam_valid_today, release_hydro_ready, release_weak_ready
      logical :: jam_formed_today
      logical :: release_ready, no_active_ice_material, open_ready
      logical :: aged_jam, post_release_flush, deepwinter_cover_ready
      logical :: ordinary_release_allowed, jam_maturity_ready
      logical :: warm_flush_set_today, mechanical_breakup_ready
      logical :: winter_drain_set_today, winter_pulse_ready
      logical :: major_bg_set_today
      logical :: winter_thermal_ready, winter_hydro_ready
      logical :: force_open_ready, breakup_tail_allowed
      logical :: breakup_phase_allowed, jam_phase_allowed, stable_relax_ready
      logical :: ice_hydro_material

      if (.not. prm_initialized) then
          call icejam_default_params(prm)
          call icejam_validate_params(prm)
          prm_initialized = .true.
      endif

      ich = j
      if (ich <= 0) return

      iwst = ob(icmd)%wst
      t_air = wst(iwst)%weat%tave
      tmax = wst(iwst)%weat%tmax
      old_phase = sd_ch(ich)%ice_phase
      new_phase = old_phase
      ht1_raw = max(0.0, ht1%flo)
      ht1_adj = ht1_raw
      major_bg_set_today = .false.

      ! New ice season reset.  The lockout prevents unrealistic refreeze after
      ! spring breakup, but it must be cleared when the following cold season starts.
      cold_window = (time%day >= prm%new_ice_year_start_day) .or. &
                    (time%day <= prm%cold_start_freezeup_end_day)
      warm_cleanup_allowed = .not. cold_window
      breakup_tail_allowed = (time%day >= prm%breakup_onset_start_day) .and. &
                             (time%day <= prm%breakup_tail_end_day)
      breakup_phase_allowed = breakup_tail_allowed
      if (time%day >= prm%new_ice_year_start_day .and. t_air <= prm%ice_frz_tmp) then
          sd_ch(ich)%ice_season_done = .false.

          ! If the previous breakup tail failed to close numerically, start a new
          ! ice year from OPEN water/cover conditions instead of carrying an old
          ! breakup episode into freezeup.
          if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
              sd_ch(ich)%ice_phase = ICE_OPEN
              sd_ch(ich)%ice_phase_days = 0
              sd_ch(ich)%is_jamming = .false.
              sd_ch(ich)%is_releasing = .false.
              sd_ch(ich)%ice_release_active = 0
              sd_ch(ich)%jam_timer = 0
              sd_ch(ich)%wedge_release_timer = 0
              sd_ch(ich)%post_release_lock_timer = 0
              sd_ch(ich)%jam_inactive_days = 0
              sd_ch(ich)%warm_flush_timer = 0
              sd_ch(ich)%winter_drain_timer = 0
              sd_ch(ich)%major_release_pending_timer = 0
              sd_ch(ich)%major_bg_timer = 0
              sd_ch(ich)%major_release_done = 0
              if (sd_ch(ich)%ice_wedge_stor > 0.0) then
                  ht1_adj = ht1_adj + sd_ch(ich)%ice_wedge_stor
                  sd_ch(ich)%ice_wedge_leak = sd_ch(ich)%ice_wedge_leak + sd_ch(ich)%ice_wedge_stor
              endif
              if (sd_ch(ich)%ice_vol > 0.0) ht1_adj = ht1_adj + sd_ch(ich)%ice_vol
              if (sd_ch(ich)%ice_mobile > 0.0) ht1_adj = ht1_adj + sd_ch(ich)%ice_mobile
              if (sd_ch(ich)%ice_mobile_pass > 0.0) ht1_adj = ht1_adj + sd_ch(ich)%ice_mobile_pass
              sd_ch(ich)%ice_vol = 0.0
              sd_ch(ich)%ice = 0.0
              sd_ch(ich)%ice_mobile = 0.0
              sd_ch(ich)%ice_mobile_pass = 0.0
              sd_ch(ich)%ice_integrity = 0.0
              sd_ch(ich)%ice_integrity_peak = 0.0
              sd_ch(ich)%ice_surface_weak = 0.0
              sd_ch(ich)%ice_surface_int = 1.0
              sd_ch(ich)%ice_wedge_stor = 0.0
              old_phase = ICE_OPEN
              new_phase = ICE_OPEN
          endif
      endif

      ! Reset daily diagnostics.
      sd_ch(ich)%icejam_block = 0.0
      sd_ch(ich)%icejam_release = 0.0
      sd_ch(ich)%icejam_qraw = 0.0
      sd_ch(ich)%icejam_qadj = 0.0
      sd_ch(ich)%icejam_qratio = 0.0
      sd_ch(ich)%icejam_qrise = 0.0
      sd_ch(ich)%icejam_susc = 0.0
      sd_ch(ich)%ice_freeze_water = 0.0
      sd_ch(ich)%ice_melt_water = 0.0
      sd_ch(ich)%ice_wedge_capture = 0.0
      sd_ch(ich)%ice_wedge_release = 0.0
      sd_ch(ich)%ice_wedge_leak = 0.0
      sd_ch(ich)%ice_excess_storage = 0.0
      sd_ch(ich)%ice_shock_release = 0.0
      sd_ch(ich)%ice_mobile_in = 0.0
      sd_ch(ich)%release_active_today = 0
      sd_ch(ich)%jam_active_today = 0
      sd_ch(ich)%major_release_today = 0
      sd_ch(ich)%warm_flush_today = 0
      warm_flush_set_today = .false.
      winter_drain_set_today = .false.
      winter_pulse_ready = .false.
      jam_formed_today = .false.
      sd_ch(ich)%ice_hydro_active = 0
      sd_ch(ich)%ice_k_mult = 1.0
      sd_ch(ich)%ice_x_current = 0.20
      sd_ch(ich)%force_eff = 0.0
      sd_ch(ich)%resistance = 0.0
      if (sd_ch(ich)%post_release_lock_timer > 0) &
          sd_ch(ich)%post_release_lock_timer = sd_ch(ich)%post_release_lock_timer - 1

      ! After the breakup-tail window, any unfinished jam episode is closed so
      ! that old event flags cannot block forced OPEN conversion.
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. .not. breakup_tail_allowed) then
          sd_ch(ich)%is_jamming = .false.
          sd_ch(ich)%is_releasing = .false.
          sd_ch(ich)%ice_release_active = 0
          sd_ch(ich)%jam_timer = 0
          sd_ch(ich)%wedge_release_timer = 0
          sd_ch(ich)%jam_inactive_days = 0
      endif

      ! Channel-scale snow/ROS/frozen-soil diagnostics from upstream HRU aggregation.
      area_ha = max(1.0e-6, sd_ch(ich)%area_ha)
      if (area_ha > 1.0e-6) then
          sd_ch(ich)%snowpack = sd_ch(ich)%snowpack_m3 / (area_ha * 10.0)
          if (sd_ch(ich)%snowpack_peak < sd_ch(ich)%snowpack) then
              sd_ch(ich)%snowpack_peak = sd_ch(ich)%snowpack
          endif
          sd_ch(ich)%snowpack_ante = 0.90 * sd_ch(ich)%snowpack_ante + &
                                      0.10 * sd_ch(ich)%snowpack
          sd_ch(ich)%frz_surf_avg = sd_ch(ich)%frz_surf_avg / area_ha
          sd_ch(ich)%frz_area_frac = sd_ch(ich)%frz_area_frac / area_ha
          sd_ch(ich)%snow_melt = sd_ch(ich)%snow_melt_m3 / (area_ha * 10.0)
          sd_ch(ich)%ros_water = sd_ch(ich)%ros_water_m3 / (area_ha * 10.0)
          ros_area_frac = sd_ch(ich)%ros_area_ha / area_ha
      else
          sd_ch(ich)%snowpack = 0.0
          sd_ch(ich)%frz_surf_avg = 0.0
          sd_ch(ich)%frz_area_frac = 0.0
          sd_ch(ich)%snow_melt = 0.0
          sd_ch(ich)%ros_water = 0.0
          ros_area_frac = 0.0
      endif
      sd_ch(ich)%snowpack = max(0.0, sd_ch(ich)%snowpack)
      sd_ch(ich)%snowpack_peak = max(0.0, sd_ch(ich)%snowpack_peak)
      sd_ch(ich)%snowpack_ante = max(0.0, sd_ch(ich)%snowpack_ante)
      sd_ch(ich)%frz_surf_avg = icejam_clamp(sd_ch(ich)%frz_surf_avg, 0.0, 1.0)
      sd_ch(ich)%frz_area_frac = icejam_clamp(sd_ch(ich)%frz_area_frac, 0.0, 1.0)
      sd_ch(ich)%ros = ros_area_frac > 0.10
      ros_water_mm = max(0.0, sd_ch(ich)%ros_water)

      ! Bankfull flow and bankfull storage.
      if (sd_ch(ich)%bankfull_flo > 1.0e-6 .and. ch_rcurv(ich)%elev(2)%flo_rate > 1.0e-6) then
          q_bankfull = sd_ch(ich)%bankfull_flo * ch_rcurv(ich)%elev(2)%flo_rate
      else if (ch_rcurv(ich)%elev(2)%flo_rate > 1.0e-6) then
          q_bankfull = ch_rcurv(ich)%elev(2)%flo_rate
      else
          q_bnk_fallback = sd_ch(ich)%chw * sd_ch(ich)%chd * &
               max(0.01, (sd_ch(ich)%chd ** (2.0/3.0)) * sqrt(max(sd_ch(ich)%chs, 1.0e-6)) / &
               max(sd_ch(ich)%chn, 0.035))
          q_bankfull = max(0.05, q_bnk_fallback)
      endif
      sd_ch(ich)%q_bankfull = q_bankfull

      ch_len_m = max(1.0, 1000.0 * sd_ch(ich)%chl)
      if (ch_rcurv(ich)%elev(2)%vol_ch > 1.0e-6) then
          bankfull_vol = ch_rcurv(ich)%elev(2)%vol_ch
      else
          bankfull_vol = max(prm%bankfull_storage_min, sd_ch(ich)%chw * sd_ch(ich)%chd * ch_len_m)
      endif
      sd_ch(ich)%bankfull_storage = bankfull_vol

      call icejam_compute_reach_scale(prm, sd_ch(ich)%chw, sd_ch(ich)%chl, sd_ch(ich)%chd, &
          sd_ch(ich)%chs, sd_ch(ich)%sinu, ch_rcurv(ich)%elev(1)%flo_rate, &
          ch_rcurv(ich)%elev(2)%flo_rate, reach)
      sd_ch(ich)%icejam_susc = reach%jam_susc
      ice_area = reach%ice_area

      q_in = ht1_raw / 86400.0
      q_prev = sd_ch(ich)%q_prev
      sd_ch(ich)%icejam_qraw = q_in
      sd_ch(ich)%icejam_qratio = q_in / max(q_bankfull, 1.0e-6)
      sd_ch(ich)%icejam_qrise = (q_in - q_prev) / max(q_bankfull, 1.0e-6)
      q_rise_pos = max(0.0, q_in - q_prev) / max(q_bankfull, 1.0e-6)
      q_fall_pos = max(0.0, q_prev - q_in) / max(q_bankfull, 1.0e-6)
      sd_ch(ich)%q_prev = q_in

      ! 1. Thermal memory, fast surface weakening, and slow structural integrity.
      ! separates two ice-strength memories:
      !   ice_surface_weak : fast thermal/ROS weakening for daily hydraulics;
      !   ice_integrity    : slow structural integrity used by phase/release gates.
      ! This prevents a 1-2 day deep-winter warm spell from destroying the
      ! structural cover, while still allowing under-ice conveyance to respond
      ! rapidly to surface melt and rain-on-snow.
      tw_ice = sd_ch(ich)%tmp_prx
      if (tw_ice < -20.0 .or. tw_ice > 40.0) tw_ice = t_air
      t_grow = min(tw_ice, t_air)
      t_decay = 0.5 * (tw_ice + t_air)

      freeze_drive = max(0.0, prm%ice_frz_tmp - t_grow)
      thaw_drive = max(0.0, t_air - prm%ice_melt_tmp)

      surface_weak_drive = prm%surface_weak_loss_thaw * thaw_drive + &
                           prm%surface_weak_loss_ros * ros_water_mm
      surface_recovery = prm%surface_weak_recovery_freeze * freeze_drive
      sd_ch(ich)%ice_surface_weak = prm%surface_weak_memory * sd_ch(ich)%ice_surface_weak + &
                                    surface_weak_drive - surface_recovery
      sd_ch(ich)%ice_surface_weak = icejam_clamp(sd_ch(ich)%ice_surface_weak, 0.0, 1.0)
      sd_ch(ich)%ice_surface_int = 1.0 - sd_ch(ich)%ice_surface_weak

      struct_gain = min(prm%structural_max_gain, prm%integrity_gain_freeze * freeze_drive)
      struct_loss = prm%integrity_loss_thaw * thaw_drive + prm%integrity_loss_ros * ros_water_mm
      struct_loss_cap = prm%structural_max_loss
      if (sd_ch(ich)%ice_phase == ICE_STABLE .and. cold_window .and. &
          sd_ch(ich)%ice_freeze_dd >= prm%major_freeze_dd_min) then
          struct_loss_cap = prm%structural_max_loss_deepwinter
      endif
      struct_loss = min(struct_loss_cap, struct_loss)
      sd_ch(ich)%ice_integrity = sd_ch(ich)%ice_integrity + struct_gain - struct_loss
      sd_ch(ich)%ice_integrity = icejam_clamp(sd_ch(ich)%ice_integrity, 0.0, 1.0)
      if (sd_ch(ich)%ice_phase == ICE_STABLE .and. cold_window .and. &
          sd_ch(ich)%ice_freeze_dd >= prm%major_freeze_dd_min .and. &
          sd_ch(ich)%ice_vol / max(ice_area, 1.0e-6) >= prm%stable_ice_thick) then
          sd_ch(ich)%ice_integrity = max(sd_ch(ich)%ice_integrity, prm%deepwinter_integrity_floor)
      endif
      sd_ch(ich)%ice_integrity_peak = max(sd_ch(ich)%ice_integrity_peak, sd_ch(ich)%ice_integrity)

      sd_ch(ich)%ice_freeze_dd = prm%freeze_memory * sd_ch(ich)%ice_freeze_dd + freeze_drive
      if (ros_water_mm > 0.0) then
          t_thaw = max(0.0, tmax - prm%thaw_tmax_base_ros)
      else
          t_thaw = max(0.0, tmax - prm%thaw_tmax_base)
      endif
      if (t_air < prm%thaw_tave_base) t_thaw = 0.5 * t_thaw
      sd_ch(ich)%ice_thaw_dd = prm%thaw_memory * sd_ch(ich)%ice_thaw_dd + t_thaw
      if (freeze_drive > 0.0) sd_ch(ich)%ice_thaw_dd = 0.7 * sd_ch(ich)%ice_thaw_dd
      if (sd_ch(ich)%ice_freeze_dd < 1.0e-6) sd_ch(ich)%ice_freeze_dd = 0.0
      if (sd_ch(ich)%ice_thaw_dd < 1.0e-6) sd_ch(ich)%ice_thaw_dd = 0.0
      thaw_weak = sd_ch(ich)%ice_thaw_dd / max(sd_ch(ich)%ice_freeze_dd + sd_ch(ich)%ice_thaw_dd, 1.0e-6)
      thaw_weak = icejam_clamp(thaw_weak, 0.0, 1.0)

      ! 2. Solid ice growth and melt.
      grow_act = 0.0
      melt_act = 0.0
      if (t_grow < prm%ice_frz_tmp) then
          ice_target_thick = min(prm%ice_maturity_ref_thick, &
              prm%ice_growth_coeff * sqrt(max(0.0, sd_ch(ich)%ice_freeze_dd)))
          ice_target_vol = ice_target_thick * ice_area
          grow_pot = max(0.0, ice_target_vol - sd_ch(ich)%ice_vol)
          grow_pot = min(grow_pot, prm%max_daily_ice_growth_thick * ice_area)
          freeze_avail = prm%max_freeze_frac_stor * max(0.0, ch_stor(ich)%flo) + &
                         prm%ice_freeze_inflow_frac * ht1_adj
          grow_act = min(grow_pot, freeze_avail)
          if (grow_act > prm%ice_min_vol) then
              freeze_from_stor = min(prm%max_freeze_frac_stor * max(0.0, ch_stor(ich)%flo), grow_act)
              ch_stor(ich)%flo = max(0.0, ch_stor(ich)%flo - freeze_from_stor)
              freeze_remain = grow_act - freeze_from_stor
              freeze_from_in = min(prm%ice_freeze_inflow_frac * ht1_adj, freeze_remain)
              ht1_adj = max(0.0, ht1_adj - freeze_from_in)
              sd_ch(ich)%ice_vol = sd_ch(ich)%ice_vol + freeze_from_stor + freeze_from_in
              sd_ch(ich)%ice_freeze_water = freeze_from_stor + freeze_from_in
          endif
      endif

      if (t_decay > prm%ice_melt_tmp .and. sd_ch(ich)%ice_vol > prm%ice_min_vol) then
          melt_thick = prm%ice_decay_coeff * (t_decay - prm%ice_melt_tmp)
          melt_pot = min(sd_ch(ich)%ice_vol, melt_thick * ice_area)
          melt_act = min(melt_pot, prm%max_melt_frac_ice * sd_ch(ich)%ice_vol)
          if (melt_act > prm%ice_min_vol) then
              sd_ch(ich)%ice_vol = max(0.0, sd_ch(ich)%ice_vol - melt_act)
              ht1_adj = ht1_adj + melt_act
              sd_ch(ich)%ice_melt_water = melt_act
          endif
      endif

      if (sd_ch(ich)%ice_vol < prm%ice_min_vol .and. warm_cleanup_allowed) sd_ch(ich)%ice_vol = 0.0
      sd_ch(ich)%ice = sd_ch(ich)%ice_vol
      sim_ice_thick = sd_ch(ich)%ice_vol / max(ice_area, 1.0e-6)
      ice_maturity = icejam_clamp(sim_ice_thick / max(prm%ice_maturity_ref_thick, 1.0e-6), 0.0, 1.0)
      if (sim_ice_thick <= prm%warm_ice_thick .and. sd_ch(ich)%ice_vol <= prm%ice_min_vol) then
          sd_ch(ich)%ice_integrity = 0.0
          sd_ch(ich)%ice_surface_weak = 0.0
          sd_ch(ich)%ice_surface_int = 1.0
      endif

      if (sd_ch(ich)%ice_phase == ICE_OPEN .and. warm_cleanup_allowed .and. &
          sd_ch(ich)%ice_vol > prm%ice_min_vol) then
          sd_ch(ich)%ice_vol = prm%ice_tail_decay * sd_ch(ich)%ice_vol
          if (sd_ch(ich)%ice_vol < prm%ice_min_vol) sd_ch(ich)%ice_vol = 0.0
          sd_ch(ich)%ice = sd_ch(ich)%ice_vol
      endif

      ! 3. Warm-season decay of mobile ice and pass-through mobile ice.
      ! intentionally reuses the ice-cover melt parameters instead of adding
      ! an independent mobile-ice melt parameter family.  Broken/mobile ice is
      ! assumed to melt or lose hydraulic relevance at least as fast as local cover
      ! ice under warm open-water conditions.
      if (warm_cleanup_allowed .and. t_decay > prm%ice_melt_tmp) then
          if (sd_ch(ich)%ice_mobile > prm%ice_min_vol) then
              mobile_melt_pot = min(sd_ch(ich)%ice_mobile, &
                  prm%ice_decay_coeff * (t_decay - prm%ice_melt_tmp) * ice_area)
              mobile_melt_act = min(mobile_melt_pot, prm%max_melt_frac_ice * sd_ch(ich)%ice_mobile)
              if (mobile_melt_act > prm%ice_min_vol) then
                  sd_ch(ich)%ice_mobile = max(0.0, sd_ch(ich)%ice_mobile - mobile_melt_act)
                  ht1_adj = ht1_adj + mobile_melt_act
                  sd_ch(ich)%ice_melt_water = sd_ch(ich)%ice_melt_water + mobile_melt_act
              endif
          endif

          if (sd_ch(ich)%ice_mobile_pass > prm%ice_min_vol) then
              pass_melt_pot = min(sd_ch(ich)%ice_mobile_pass, &
                  prm%ice_decay_coeff * (t_decay - prm%ice_melt_tmp) * ice_area)
              pass_melt_act = min(pass_melt_pot, prm%max_melt_frac_ice * sd_ch(ich)%ice_mobile_pass)
              if (pass_melt_act > prm%ice_min_vol) then
                  sd_ch(ich)%ice_mobile_pass = max(0.0, sd_ch(ich)%ice_mobile_pass - pass_melt_act)
                  ht1_adj = ht1_adj + pass_melt_act
                  sd_ch(ich)%ice_melt_water = sd_ch(ich)%ice_melt_water + pass_melt_act
              endif
          endif

          if (sd_ch(ich)%ice_mobile < prm%ice_min_vol) sd_ch(ich)%ice_mobile = 0.0
          if (sd_ch(ich)%ice_mobile_pass < prm%ice_min_vol) sd_ch(ich)%ice_mobile_pass = 0.0
      endif

      ! 4. Incoming mobile ice assimilation.  Captured mobile ice remains mobile
      ! jam material unless an active jam is forming.
      if (sd_ch(ich)%ice_mobile > prm%ice_min_vol) then
          sd_ch(ich)%ice_mobile_in = sd_ch(ich)%ice_mobile
          if (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. sd_ch(ich)%is_jamming) then
              mobile_assim = min(sd_ch(ich)%ice_mobile, &
                  prm%jam_mobile_excess_capture_frac * sd_ch(ich)%ice_mobile)
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile - mobile_assim
              sd_ch(ich)%ice_vol = sd_ch(ich)%ice_vol + mobile_assim
              sd_ch(ich)%ice = sd_ch(ich)%ice_vol
          endif
      endif

      ! 5. Hydraulic forcing and resistance.
      local_shock_q = ros_water_mm * area_ha / 8640.0
      thermal_shock = max(0.0, tmax - prm%breakup_tmax_base) + &
                      prm%tave_weight * max(0.0, t_air - prm%breakup_tave_base)
      damage_factor = sd_ch(ich)%ice_surface_weak
      force_mult = 1.0 + prm%thermal_force_weight * thermal_shock + &
                   prm%damage_force_weight * damage_factor
      force_F = (q_in + prm%shock_lambda * local_shock_q) * max(1.0, force_mult)

      mobile_factor = min(1.0, sd_ch(ich)%ice_mobile / max(prm%reference_ice_vol, 1.0e-6))
      alpha_ice = prm%alpha_min + (prm%alpha_max - prm%alpha_min) * sd_ch(ich)%ice_integrity
      resistance_R = q_bankfull * alpha_ice * (1.0 + prm%mobile_resistance_weight * mobile_factor)
      sd_ch(ich)%force_eff = force_F
      sd_ch(ich)%resistance = resistance_R
      fr_ratio = force_F / max(resistance_R, 1.0e-6)

      wedge_capacity = max(prm%bankfull_storage_min, &
          prm%wedge_capacity_bankfull_mult * bankfull_vol * reach%jam_storage_modifier * &
          (1.0 + prm%mobile_wedge_capacity_weight * mobile_factor))
      sd_ch(ich)%ice_wedge_capacity = wedge_capacity

      ! preliminary major-release gate, computed before the phase machine
      ! so STABLE mechanical breakup can use the same deep-winter/snow/frozen/
      ! hydro background factors.  The full values are recomputed after mobile-
      ! ice generation for diagnostics and release decisions.
      deep_winter_factor = icejam_clamp(sd_ch(ich)%ice_freeze_dd / &
                           max(prm%major_freeze_dd_min, 1.0e-6), 0.0, 1.0)
      ice_storage_factor = min(icejam_clamp(ice_maturity / &
                           max(prm%major_ice_maturity_min, 1.0e-6), 0.0, 1.0), &
                           icejam_clamp(sd_ch(ich)%ice_integrity_peak / &
                           max(prm%major_integrity_peak_min, 1.0e-6), 0.0, 1.0))
      snowpack_ante_factor = icejam_clamp(sd_ch(ich)%snowpack_ante / &
                             max(prm%major_snowpack_ante_min, 1.0e-6), 0.0, 1.0)
      snowpack_peak_factor = icejam_clamp(sd_ch(ich)%snowpack_peak / &
                             max(prm%major_snowpack_peak_min, 1.0e-6), 0.0, 1.0)
      snowpack_factor = min(snowpack_ante_factor, snowpack_peak_factor)
      frozen_soil_factor = max(icejam_clamp(sd_ch(ich)%frz_surf_avg / &
                           max(prm%major_frz_surf_min, 1.0e-6), 0.0, 1.0), &
                           icejam_clamp(sd_ch(ich)%frz_area_frac / &
                           max(prm%major_frz_area_min, 1.0e-6), 0.0, 1.0))
      warm_air_factor = max(icejam_clamp((tmax - prm%major_warm_tmax_min) / 4.0, 0.0, 1.0), &
                        icejam_clamp((t_air - prm%major_warm_tave_min) / 4.0, 0.0, 1.0))
      meltwater_factor = max(icejam_clamp(sd_ch(ich)%snow_melt / max(prm%major_snomelt_min, 1.0e-6), 0.0, 1.0), &
                         icejam_clamp(ros_water_mm / max(prm%major_ros_min, 1.0e-6), 0.0, 1.0))
      warm_flush_factor = min(warm_air_factor, meltwater_factor)
      fr_factor = icejam_clamp(fr_ratio / max(prm%major_fr_min, 1.0e-6), 0.0, 1.0)
      qrise_factor = icejam_clamp(q_rise_pos / max(prm%major_qrise_min, 1.0e-6), 0.0, 1.0)
      runoff_response_factor = max(qrise_factor, meltwater_factor)
      discharge_factor = min(fr_factor, runoff_response_factor)
      warm_flush_ready = warm_flush_factor >= 1.0
      major_background_factor = min(min(deep_winter_factor, ice_storage_factor), &
                                min(snowpack_factor, frozen_soil_factor))
      major_background_factor = icejam_clamp(major_background_factor, 0.0, 1.0)
      warm_memory_factor = warm_flush_factor
      if (sd_ch(ich)%warm_flush_timer > 0) warm_memory_factor = 1.0
      major_trigger_factor = min(warm_memory_factor, discharge_factor)
      major_trigger_factor = icejam_clamp(major_trigger_factor, 0.0, 1.0)
      major_base_factor = min(major_background_factor, discharge_factor)
      major_base_factor = icejam_clamp(major_base_factor, 0.0, 1.0)
      if ((sd_ch(ich)%ice_phase == ICE_STABLE .or. sd_ch(ich)%ice_phase == ICE_BREAKUP) .and. &
          sd_ch(ich)%major_release_done == 0 .and. warm_flush_ready .and. &
          major_background_factor >= prm%warm_flush_memory_base_min) then
          sd_ch(ich)%warm_flush_today = 1
          sd_ch(ich)%warm_flush_timer = max(sd_ch(ich)%warm_flush_timer, &
                                            prm%warm_flush_release_days)
          warm_flush_set_today = .true.
      endif

      ! 6. Seasonal phase machine.  It is intentionally only one-way within an
      ! ice year; dynamic jam/release states are not phases.
      select case (sd_ch(ich)%ice_phase)
      case (ICE_OPEN)
          sd_ch(ich)%is_jamming = .false.
          sd_ch(ich)%is_releasing = .false.
          sd_ch(ich)%jam_timer = 0
          freezeup_ready = (.not. sd_ch(ich)%ice_season_done) .and. &
              ((cold_window) .or. (sd_ch(ich)%ice_phase_days >= prm%warm_min_days_before_freezeup)) .and. &
              sd_ch(ich)%ice_freeze_dd >= prm%freezeup_freeze_dd .and. &
              sim_ice_thick >= prm%freezeup_ice_thick .and. &
              thaw_weak <= prm%freezeup_strong_index
          if (freezeup_ready) then
              new_phase = ICE_FREEZEUP
              sd_ch(ich)%snowpack_peak = sd_ch(ich)%snowpack
              sd_ch(ich)%snowpack_ante = sd_ch(ich)%snowpack
              sd_ch(ich)%ice_integrity_peak = sd_ch(ich)%ice_integrity
              sd_ch(ich)%warm_flush_timer = 0
              sd_ch(ich)%winter_drain_timer = 0
              sd_ch(ich)%major_release_pending_timer = 0
              sd_ch(ich)%major_bg_timer = 0
              sd_ch(ich)%major_release_done = 0
          endif

      case (ICE_FREEZEUP)
          ! FREEZEUP is a seasonal background phase.  Jam/release events may
          ! occur inside it, but post-release lockout must not indefinitely
          ! block the background transition to stable cover.  Only an active
          ! release episode prevents STABLE conversion.
          stable_ready = .not. sd_ch(ich)%is_releasing .and. &
                         ((sd_ch(ich)%ice_phase_days >= prm%freezeup_min_days .and. &
                         sim_ice_thick >= prm%stable_ice_thick .and. &
                         thaw_weak <= prm%freezeup_strong_index) .or. &
                         (sd_ch(ich)%ice_phase_days >= prm%freezeup_max_days .and. &
                         sim_ice_thick >= prm%freezeup_ice_thick))
          if (stable_ready) new_phase = ICE_STABLE

      case (ICE_STABLE)
          ! STABLE is a seasonal background state.  High flow, snowmelt, or a
          ! short warm spell may weaken the cover and relax its hydraulic
          ! restriction, but it must not by itself trigger BREAKUP or releasing.
          ! STABLE -> BREAKUP is allowed only inside the spring breakup window.
          breakup_material_ready = sim_ice_thick >= prm%warm_ice_thick .or. &
                                   sd_ch(ich)%ice_vol > prm%ice_min_vol .or. &
                                   sd_ch(ich)%ice_wedge_stor > prm%warm_storage_exit_ratio * wedge_capacity
          breakup_weather_ready = (tmax >= prm%breakup_tmax_min) .or. &
                                  (t_air >= prm%breakup_tave_min) .or. &
                                  (ros_water_mm >= prm%breakup_ros_min)
          breakup_age_ready = sd_ch(ich)%ice_phase_days >= prm%stable_min_days_before_breakup .or. &
                              sd_ch(ich)%ice_phase_days >= prm%stable_max_days_before_breakup
          mechanical_breakup_ready = breakup_phase_allowed .and. &
              breakup_material_ready .and. breakup_age_ready .and. &
              (sd_ch(ich)%warm_flush_timer > 0 .or. warm_flush_ready .or. &
               warm_flush_factor >= prm%mechanical_breakup_warm_min) .and. &
              major_background_factor >= prm%mechanical_breakup_base_min .and. &
              sd_ch(ich)%ice_surface_weak >= prm%mechanical_breakup_surface_weak_min .and. &
              fr_ratio >= prm%breakup_force_ratio_min
          breakup_ready = breakup_phase_allowed .and. &
                          breakup_material_ready .and. breakup_weather_ready .and. &
                          breakup_age_ready .and. &
                          ((thaw_weak >= prm%breakup_onset_weakening_index .or. &
                           sd_ch(ich)%ice_integrity <= prm%release_integrity_max) .or. &
                           mechanical_breakup_ready)
          if (breakup_ready) then
              new_phase = ICE_BREAKUP
              sd_ch(ich)%is_jamming = .false.
              sd_ch(ich)%is_releasing = .false.
              sd_ch(ich)%ice_release_active = 0
              sd_ch(ich)%jam_timer = 0
              sd_ch(ich)%wedge_release_timer = 0
              sd_ch(ich)%breakup_intensity = icejam_clamp((force_F - resistance_R) / max(q_bankfull, 1.0e-6), 0.0, 1.0)
          endif

      case (ICE_BREAKUP)
          ! No transition back to STABLE.  Leave breakup only after the warm-season
          ! cleanup gate and after active ice material has mostly disappeared.
          no_active_ice_material = (sd_ch(ich)%ice_vol <= prm%jam_material_min_vol) .and. &
                                   (sd_ch(ich)%ice_mobile <= prm%jam_mobile_min_vol) .and. &
                                   (sd_ch(ich)%ice_mobile_pass <= prm%jam_mobile_min_vol)
          open_ready = time%day >= prm%breakup_warm_exit_start_day .and. &
                       sd_ch(ich)%ice_phase_days >= prm%breakup_min_days_before_open .and. &
                       thaw_weak >= prm%warm_season_weakening_index .and. &
                       sd_ch(ich)%ice_thaw_dd >= prm%flush_thaw_dd .and. &
                       sim_ice_thick <= prm%warm_ice_thick .and. &
                       sd_ch(ich)%ice_wedge_stor <= prm%warm_storage_exit_ratio * wedge_capacity .and. &
                       no_active_ice_material
          if (sd_ch(ich)%ice_phase_days >= prm%breakup_max_days_before_open .and. &
              sim_ice_thick <= prm%warm_ice_thick .and. no_active_ice_material) open_ready = .true.
          force_open_ready = time%day >= prm%breakup_force_open_day .and. &
                             sim_ice_thick <= prm%force_open_ice_thick .and. &
                             sd_ch(ich)%ice_integrity <= prm%force_open_integrity
          if ((open_ready .or. force_open_ready) .and. .not. sd_ch(ich)%is_jamming .and. &
              .not. sd_ch(ich)%is_releasing) new_phase = ICE_OPEN

      case default
          new_phase = ICE_OPEN
      end select

      if (new_phase /= sd_ch(ich)%ice_phase) then
          sd_ch(ich)%ice_phase = new_phase
          sd_ch(ich)%ice_phase_days = 1
      else
          sd_ch(ich)%ice_phase_days = sd_ch(ich)%ice_phase_days + 1
      endif
      if (old_phase == ICE_FREEZEUP .and. sd_ch(ich)%ice_phase == ICE_STABLE) then
          ! A freeze-up jam is assimilated into the stable-cover background.
          ! It should not keep the seasonal phase in FREEZEUP indefinitely.
          sd_ch(ich)%is_jamming = .false.
          sd_ch(ich)%jam_timer = 0
          sd_ch(ich)%jam_inactive_days = 0
      endif

      if (old_phase /= ICE_OPEN .and. sd_ch(ich)%ice_phase == ICE_OPEN) then
          sd_ch(ich)%ice_season_done = .true.
          sd_ch(ich)%is_jamming = .false.
          sd_ch(ich)%is_releasing = .false.
          sd_ch(ich)%jam_timer = 0
          sd_ch(ich)%wedge_release_timer = 0
          sd_ch(ich)%post_release_lock_timer = 0
          sd_ch(ich)%jam_inactive_days = 0
          sd_ch(ich)%warm_flush_timer = 0
          ! OPEN cleanup is mass-conservative: all remaining equivalent ice
          ! volume is returned to the current channel inflow before states are
          ! cleared.
          if (sd_ch(ich)%ice_mobile > 0.0) then
              ht1_adj = ht1_adj + sd_ch(ich)%ice_mobile
              sd_ch(ich)%ice_melt_water = sd_ch(ich)%ice_melt_water + sd_ch(ich)%ice_mobile
              sd_ch(ich)%ice_mobile = 0.0
          endif
          if (sd_ch(ich)%ice_mobile_pass > 0.0) then
              ht1_adj = ht1_adj + sd_ch(ich)%ice_mobile_pass
              sd_ch(ich)%ice_melt_water = sd_ch(ich)%ice_melt_water + sd_ch(ich)%ice_mobile_pass
              sd_ch(ich)%ice_mobile_pass = 0.0
          endif
          if (sd_ch(ich)%ice_vol > 0.0) then
              ht1_adj = ht1_adj + sd_ch(ich)%ice_vol
              sd_ch(ich)%ice_melt_water = sd_ch(ich)%ice_melt_water + sd_ch(ich)%ice_vol
              sd_ch(ich)%ice_vol = 0.0
              sd_ch(ich)%ice = 0.0
          endif
          sd_ch(ich)%ice_integrity = 0.0
      endif

      ! 7. Mobile ice generation and repeatable jam episode state.
      ! Jam/release events are allowed in FREEZEUP and in the spring BREAKUP
      ! window.  They do not change the seasonal phase.
      mobile_thermal = 0.0
      mobile_mech = 0.0
      jam_phase_allowed = (sd_ch(ich)%ice_phase == ICE_FREEZEUP) .or. &
                          (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. breakup_tail_allowed)
      if (jam_phase_allowed .and. .not. sd_ch(ich)%is_releasing .and. &
          sd_ch(ich)%ice_vol > prm%ice_min_vol) then
          if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
              if (thaw_weak >= prm%breakup_onset_weakening_index .or. sd_ch(ich)%snow_melt > 0.0) then
                  mobile_thermal = prm%thermal_mobile_frac * min(1.0, max(thaw_weak, sd_ch(ich)%snow_melt / 25.0)) * &
                                   sd_ch(ich)%ice_vol
              endif
          else if (sd_ch(ich)%ice_phase == ICE_FREEZEUP) then
              mobile_thermal = 0.50 * prm%thermal_mobile_frac * min(1.0, freeze_drive / 5.0) * &
                               sd_ch(ich)%ice_vol
          endif
          if (fr_ratio > prm%mechanical_breakup_fr) then
              mobile_mech = prm%mechanical_mobile_frac * &
                            min(1.0, (fr_ratio - prm%mechanical_breakup_fr) / &
                                     max(prm%mechanical_breakup_fr_scale, 1.0e-6)) * sd_ch(ich)%ice_vol
          endif
          mobile_gen = min(sd_ch(ich)%ice_vol, max(0.0, mobile_thermal) + max(0.0, mobile_mech))
          if (mobile_gen > prm%ice_min_vol) then
              sd_ch(ich)%ice_vol = sd_ch(ich)%ice_vol - mobile_gen
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile + mobile_gen
              sd_ch(ich)%ice = sd_ch(ich)%ice_vol
          endif
      endif

      q_rel = q_in / max(q_bankfull, 1.0e-6)
      trans_flow_factor = 0.2 + max(0.0, q_rel) ** prm%ice_transport_q_exp
      trans_susc_factor = max(0.10, 1.0 - prm%jam_susc_transport_weight * reach%jam_susc)
      ice_transport_cap = max(prm%ice_transport_cap_min, prm%ice_transport_cap_base * trans_flow_factor * trans_susc_factor)
      mobile_supply = sd_ch(ich)%ice_mobile + sd_ch(ich)%ice_mobile_pass
      ice_mobile_excess = mobile_supply - ice_transport_cap
      tail_material_factor = icejam_clamp(max(0.0, ice_mobile_excess) / &
                             max(prm%reference_ice_vol, 1.0e-6), 0.0, 1.0)

      ! major-release gate.  Jam formation remains generic.  Strong
      ! jam-break release requires deep winter, strong channel ice, antecedent
      ! AND peak snowpack, frozen soil, hydraulic response, and a remembered
      ! warm-flush event.
      deep_winter_factor = icejam_clamp(sd_ch(ich)%ice_freeze_dd / &
                           max(prm%major_freeze_dd_min, 1.0e-6), 0.0, 1.0)
      ice_storage_factor = min(icejam_clamp(ice_maturity / &
                           max(prm%major_ice_maturity_min, 1.0e-6), 0.0, 1.0), &
                           icejam_clamp(sd_ch(ich)%ice_integrity_peak / &
                           max(prm%major_integrity_peak_min, 1.0e-6), 0.0, 1.0))
      snowpack_ante_factor = icejam_clamp(sd_ch(ich)%snowpack_ante / &
                             max(prm%major_snowpack_ante_min, 1.0e-6), 0.0, 1.0)
      snowpack_peak_factor = icejam_clamp(sd_ch(ich)%snowpack_peak / &
                             max(prm%major_snowpack_peak_min, 1.0e-6), 0.0, 1.0)
      snowpack_factor = min(snowpack_ante_factor, snowpack_peak_factor)
      frozen_soil_factor = max(icejam_clamp(sd_ch(ich)%frz_surf_avg / &
                           max(prm%major_frz_surf_min, 1.0e-6), 0.0, 1.0), &
                           icejam_clamp(sd_ch(ich)%frz_area_frac / &
                           max(prm%major_frz_area_min, 1.0e-6), 0.0, 1.0))
      warm_air_factor = max(icejam_clamp((tmax - prm%major_warm_tmax_min) / 4.0, 0.0, 1.0), &
                        icejam_clamp((t_air - prm%major_warm_tave_min) / 4.0, 0.0, 1.0))
      meltwater_factor = max(icejam_clamp(sd_ch(ich)%snow_melt / max(prm%major_snomelt_min, 1.0e-6), 0.0, 1.0), &
                         icejam_clamp(ros_water_mm / max(prm%major_ros_min, 1.0e-6), 0.0, 1.0))
      warm_flush_factor = min(warm_air_factor, meltwater_factor)
      fr_factor = icejam_clamp(fr_ratio / max(prm%major_fr_min, 1.0e-6), 0.0, 1.0)
      qrise_factor = icejam_clamp(q_rise_pos / max(prm%major_qrise_min, 1.0e-6), 0.0, 1.0)
      runoff_response_factor = max(qrise_factor, meltwater_factor)
      discharge_factor = min(fr_factor, runoff_response_factor)

      deep_winter_ready = deep_winter_factor >= 1.0
      channel_ice_ready = ice_storage_factor >= 1.0
      snowpack_ready = snowpack_factor >= 1.0
      frozen_soil_ready = frozen_soil_factor >= 1.0
      warm_flush_ready = warm_flush_factor >= 1.0
      discharge_ready = discharge_factor >= 1.0
      major_background_factor = min(min(deep_winter_factor, ice_storage_factor), &
                                min(snowpack_factor, frozen_soil_factor))
      major_background_factor = icejam_clamp(major_background_factor, 0.0, 1.0)
      warm_memory_factor = warm_flush_factor
      if (sd_ch(ich)%warm_flush_timer > 0) warm_memory_factor = 1.0
      major_trigger_factor = min(warm_memory_factor, discharge_factor)
      major_trigger_factor = icejam_clamp(major_trigger_factor, 0.0, 1.0)
      major_base_factor = min(major_background_factor, discharge_factor)
      major_base_factor = icejam_clamp(major_base_factor, 0.0, 1.0)

      ! Warm-flush is an event memory, not only a same-day diagnostic.  A strong
      ! warm-air + meltwater/ROS trigger can support the following 1-3 release
      ! days, which is required for daily time-step representation of jam failure.
      if ((sd_ch(ich)%ice_phase == ICE_STABLE .or. sd_ch(ich)%ice_phase == ICE_BREAKUP) .and. &
          sd_ch(ich)%major_release_done == 0 .and. warm_flush_ready .and. &
          major_background_factor >= prm%warm_flush_memory_base_min) then
          sd_ch(ich)%warm_flush_today = 1
          sd_ch(ich)%warm_flush_timer = max(sd_ch(ich)%warm_flush_timer, &
                                            prm%warm_flush_release_days)
          warm_flush_set_today = .true.
      endif
      warm_memory_factor = warm_flush_factor
      if (sd_ch(ich)%warm_flush_timer > 0) warm_memory_factor = 1.0
      major_trigger_factor = min(warm_memory_factor, discharge_factor)
      major_trigger_factor = icejam_clamp(major_trigger_factor, 0.0, 1.0)
      major_jam_factor = min(major_background_factor, major_trigger_factor)
      major_jam_factor = icejam_clamp(major_jam_factor, 0.0, 1.0)
      sd_ch(ich)%major_jam_factor = major_jam_factor

      ! compute a dedicated storage ratio here.  Do not use
      ! wedge_ratio before its later diagnostic/leakage assignments.
      ! This scale-free major-storage readiness prevents small-wedge ordinary
      ! spring ice effects from being upgraded to major release candidates.
      wedge_ratio_for_major = sd_ch(ich)%ice_wedge_stor / max(wedge_capacity, 1.0e-6)
      major_storage_ready = wedge_ratio_for_major >= prm%major_wedge_ratio_min

      ! retain a short memory of slow major-event background
      ! conditions before the warm-flush/hydraulic trigger arrives.  This
      ! memory does not trigger major release by itself; release-day warm
      ! memory and F/R are still required.  It prevents events such as 2011
      ! from losing the major candidate simply because ice/snow background
      ! variables are partly consumed by the first warm/hydraulic pulse.
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. &
          sd_ch(ich)%major_release_done == 0 .and. &
          major_background_factor >= 0.999 .and. major_storage_ready .and. &
          (sd_ch(ich)%is_jamming .or. mobile_supply >= prm%jam_mobile_min_vol)) then
          sd_ch(ich)%major_bg_timer = max(sd_ch(ich)%major_bg_timer, &
                                          prm%major_release_pending_days)
          major_bg_set_today = .true.
      endif

      ! Keep the older pending diagnostic as a stricter subset: background and
      ! warm-flush memory are present, but hydraulic force is not yet large
      ! enough.  It remains useful for debugging and ordinary-release
      ! suppression, while major_bg_timer carries the slower background memory.
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. &
          sd_ch(ich)%warm_flush_timer > 0 .and. &
          sd_ch(ich)%major_release_done == 0 .and. &
          major_background_factor >= 0.999 .and. &
          fr_ratio < prm%major_release_start_fr_min) then
          sd_ch(ich)%major_release_pending_timer = max( &
              sd_ch(ich)%major_release_pending_timer, &
              prm%major_release_pending_days)
      endif

      major_release_gate = sd_ch(ich)%ice_phase == ICE_BREAKUP .and. &
                           sd_ch(ich)%warm_flush_timer > 0 .and. &
                           sd_ch(ich)%major_release_done == 0 .and. &
                           fr_ratio >= prm%major_release_start_fr_min .and. &
                           (major_background_factor >= 0.999 .or. &
                            sd_ch(ich)%major_bg_timer > 0)

      jam_material_ready = mobile_supply >= prm%jam_mobile_min_vol
      jam_transport_ready = ice_mobile_excess > 0.0
      jam_storage_ready = sd_ch(ich)%ice_wedge_stor >= prm%wedge_release_storage_frac * wedge_capacity
      jam_valid_today = (sd_ch(ich)%ice_vol >= prm%jam_material_min_vol) .or. &
                        (sd_ch(ich)%ice_mobile >= prm%jam_mobile_min_vol) .or. &
                        (sd_ch(ich)%ice_mobile_pass >= prm%jam_mobile_min_vol)

      if (jam_phase_allowed) then
          if (sd_ch(ich)%is_jamming) then
              sd_ch(ich)%jam_active_today = 1
              sd_ch(ich)%jam_timer = sd_ch(ich)%jam_timer + 1
              if (jam_valid_today) then
                  sd_ch(ich)%jam_inactive_days = 0
              else
                  sd_ch(ich)%jam_inactive_days = sd_ch(ich)%jam_inactive_days + 1
              endif
              if (sd_ch(ich)%jam_inactive_days >= prm%jam_inactive_max_days) then
                  sd_ch(ich)%is_jamming = .false.
                  sd_ch(ich)%jam_timer = 0
                  sd_ch(ich)%jam_inactive_days = 0
              endif
          else if (.not. sd_ch(ich)%is_releasing) then
              if (jam_phase_allowed .and. sd_ch(ich)%post_release_lock_timer <= 0 .and. &
                  jam_material_ready .and. jam_transport_ready .and. &
                  reach%jam_susc >= prm%jam_susc_min) then
                  sd_ch(ich)%is_jamming = .true.
                  sd_ch(ich)%jam_active_today = 1
                  jam_formed_today = .true.
                  sd_ch(ich)%jam_timer = 1
                  sd_ch(ich)%jam_inactive_days = 0
                  mobile_assim = min(mobile_supply, ice_mobile_excess * prm%jam_mobile_excess_capture_frac)
                  if (mobile_assim > prm%ice_min_vol) then
                      mobile_remain = mobile_assim
                      if (sd_ch(ich)%ice_mobile > 0.0) then
                          mobile_gen = min(sd_ch(ich)%ice_mobile, mobile_remain)
                          sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile - mobile_gen
                          mobile_remain = mobile_remain - mobile_gen
                      endif
                      if (mobile_remain > 0.0 .and. sd_ch(ich)%ice_mobile_pass > 0.0) then
                          mobile_gen = min(sd_ch(ich)%ice_mobile_pass, mobile_remain)
                          sd_ch(ich)%ice_mobile_pass = sd_ch(ich)%ice_mobile_pass - mobile_gen
                      endif
                      sd_ch(ich)%ice_vol = sd_ch(ich)%ice_vol + mobile_assim
                      sd_ch(ich)%ice = sd_ch(ich)%ice_vol
                  endif
              endif
          endif

          wedge_ratio = sd_ch(ich)%ice_wedge_stor / max(wedge_capacity, 1.0e-6)
          aged_jam = sd_ch(ich)%is_jamming .and. sd_ch(ich)%jam_timer >= prm%max_jam_days
          release_force_eff = prm%release_force_ratio
          if (aged_jam) release_force_eff = 0.75 * prm%release_force_ratio
          if (major_release_gate) release_force_eff = 0.70 * prm%release_force_ratio

          ! max_jam_days weakens the jam and enhances leakage, but it is no
          ! longer a sufficient condition for a dam-break release.  A release
          ! still requires hydraulic forcing and/or a large impoundment under
          ! physically plausible forcing.
          release_hydro_ready = fr_ratio >= release_force_eff .or. &
                                (jam_storage_ready .and. fr_ratio >= 0.50 * release_force_eff)
          if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
              release_hydro_ready = release_hydro_ready .or. &
                                    (jam_storage_ready .and. &
                                     thaw_weak >= prm%jam_release_weakening_index .and. &
                                     fr_ratio >= 0.25 * release_force_eff)
          endif
          release_weak_ready = thaw_weak >= prm%jam_release_weakening_index .or. &
                               sd_ch(ich)%ice_integrity <= prm%release_integrity_max .or. &
                               t_air >= prm%breakup_tave_min
          ordinary_release_allowed = .true.
          if ((sd_ch(ich)%major_release_pending_timer > 0 .or. &
              (sd_ch(ich)%major_bg_timer > 0 .and. sd_ch(ich)%warm_flush_timer > 0)) .and. &
              .not. major_release_gate) then
              ordinary_release_allowed = .false.
          endif
          if (sd_ch(ich)%major_release_done == 1 .and. &
              sd_ch(ich)%ice_phase == ICE_BREAKUP .and. .not. major_release_gate) then
              ordinary_release_allowed = fr_ratio >= 2.0 * release_force_eff .and. &
                                         mobile_supply >= 2.0 * prm%jam_mobile_min_vol .and. &
                                         ice_mobile_excess >= 2.0 * prm%ice_transport_cap_min
          endif

          ! ordinary release requires a more mature jam.
          ! Major release can still occur after the configured minimum jam age;
          ! ordinary release must wait until the jam_timer reaches at least 3.
          ! Because jam_timer is initialized to 1 on the formation day, this
          ! prevents a newly entered BREAKUP phase from producing a next-day
          ! minor dam-break that behaves like a weak major release.
          if (major_release_gate) then
              jam_maturity_ready = sd_ch(ich)%jam_timer >= prm%min_jam_days
          else
              jam_maturity_ready = sd_ch(ich)%jam_timer >= max(prm%min_jam_days, 3)
          endif

          release_ready = jam_phase_allowed .and. ordinary_release_allowed .and. &
                          (sd_ch(ich)%ice_phase == ICE_FREEZEUP .or. &
                           time%day >= prm%breakup_release_start_day) .and. &
                          sd_ch(ich)%is_jamming .and. &
                          .not. jam_formed_today .and. &
                          jam_maturity_ready .and. &
                          release_hydro_ready .and. release_weak_ready
          if (release_ready) then
              sd_ch(ich)%is_jamming = .false.
              sd_ch(ich)%is_releasing = .true.
              sd_ch(ich)%wedge_release_timer = 0
              sd_ch(ich)%jam_timer = 0
              if (major_release_gate) then
                  sd_ch(ich)%ice_release_active = 2
                  sd_ch(ich)%major_release_pending_timer = 0
                  sd_ch(ich)%major_bg_timer = 0
              else
                  sd_ch(ich)%ice_release_active = 1
              endif
              sd_ch(ich)%breakup_intensity = icejam_clamp((force_F - resistance_R) / max(q_bankfull, 1.0e-6), 0.0, 1.0)
          endif
      endif

      ! 8. Wedge storage capture/release.  This is the only conceptual liquid storage.
      ! wedge_release is permitted only during a
      ! jam-break release episode.  All non-event drainage is wedge_leak.
      if (jam_phase_allowed .and. sd_ch(ich)%is_releasing) then
          ! This block must only be reached after release_ready set is_releasing
          ! from a pre-existing jam.  It is intentionally separate from daily leakage.
          major_release_active = sd_ch(ich)%ice_release_active == 2
          major_release_cap = prm%wedge_release_max
          if (major_release_active) then
              major_release_cap = prm%major_release_max
          endif
          wedge_release_frac = prm%wedge_release_min + &
              (major_release_cap - prm%wedge_release_min) * sd_ch(ich)%breakup_intensity
          wedge_release_frac = icejam_clamp(wedge_release_frac, prm%wedge_release_min, major_release_cap)
          if (sd_ch(ich)%ice_phase == ICE_FREEZEUP) then
              wedge_release_frac = prm%freezeup_release_frac * wedge_release_frac
          endif
          if (major_release_active) then
              ! a major jam-break is a short multi-day burst, not a
              ! same-day complete emptying of the wedge.  The configured
              ! major_release_max controls the per-day fraction; storage
              ! boost defaults to zero so the 3-day episode drains the wedge
              ! progressively instead of flushing it all at once.
              wedge_release = min(sd_ch(ich)%ice_wedge_stor, &
                              wedge_release_frac * sd_ch(ich)%ice_wedge_stor * &
                              (1.0 + prm%major_release_storage_boost))
          else
              ! Ordinary jam-break release must not behave like a major burst.
              ! It can drain only a moderate fraction of the existing wedge and
              ! is also capped by reach-scaled wedge capacity.  This prevents
              ! large winter-accumulated wedge storage from being released by a
              ! minor/residual jam episode.
              ordinary_release_cap = min(prm%ordinary_release_max_frac * &
                                         sd_ch(ich)%ice_wedge_stor, &
                                         prm%ordinary_release_capacity_frac * &
                                         wedge_capacity)
              wedge_release = min(sd_ch(ich)%ice_wedge_stor, &
                                  wedge_release_frac * sd_ch(ich)%ice_wedge_stor, &
                                  ordinary_release_cap)
          endif
          if (wedge_release > 1.0e-6) then
              ht1_adj = ht1_adj + wedge_release
              sd_ch(ich)%ice_wedge_stor = sd_ch(ich)%ice_wedge_stor - wedge_release
              sd_ch(ich)%release_active_today = 1
              if (major_release_active) sd_ch(ich)%major_release_today = 1
          endif
          sd_ch(ich)%ice_wedge_release = wedge_release
          sd_ch(ich)%ice_shock_release = wedge_release
          sd_ch(ich)%icejam_release = wedge_release

          if (wedge_release > 1.0e-6 .and. sd_ch(ich)%ice_vol > prm%ice_min_vol) then
              mobile_gen = min(sd_ch(ich)%ice_vol, prm%release_ice_to_mobile_frac * wedge_release_frac * sd_ch(ich)%ice_vol)
              if (mobile_gen > prm%ice_min_vol) then
                  sd_ch(ich)%ice_vol = sd_ch(ich)%ice_vol - mobile_gen
                  sd_ch(ich)%ice_mobile_pass = sd_ch(ich)%ice_mobile_pass + mobile_gen
                  sd_ch(ich)%ice = sd_ch(ich)%ice_vol
              endif
          endif

          sd_ch(ich)%wedge_release_timer = sd_ch(ich)%wedge_release_timer + 1
          ! ordinary release is a one-day local event.  Only
          ! major release uses the configured multi-day burst duration.
          if (((.not. major_release_active) .and. sd_ch(ich)%wedge_release_timer >= 1) .or. &
              (major_release_active .and. sd_ch(ich)%wedge_release_timer >= prm%release_duration_days) .or. &
              sd_ch(ich)%ice_wedge_stor <= prm%warm_storage_exit_ratio * wedge_capacity) then
              sd_ch(ich)%is_releasing = .false.
              sd_ch(ich)%ice_integrity = max(0.0, 0.35 * sd_ch(ich)%ice_integrity)
              sd_ch(ich)%ice_surface_weak = max(sd_ch(ich)%ice_surface_weak, 0.75)
              sd_ch(ich)%ice_surface_int = 1.0 - sd_ch(ich)%ice_surface_weak
              sd_ch(ich)%ice_vol = prm%post_release_ice_retention * sd_ch(ich)%ice_vol
              if (sd_ch(ich)%ice_vol < prm%ice_min_vol .and. warm_cleanup_allowed) sd_ch(ich)%ice_vol = 0.0
              sd_ch(ich)%ice = sd_ch(ich)%ice_vol
              sd_ch(ich)%wedge_release_timer = 0
              if (sd_ch(ich)%ice_release_active == 2) then
                  sd_ch(ich)%post_release_lock_timer = prm%major_post_release_lock_days
                  sd_ch(ich)%major_release_done = 1
                  sd_ch(ich)%major_release_pending_timer = 0
                  sd_ch(ich)%major_bg_timer = 0
                  sd_ch(ich)%warm_flush_timer = 0
              else
                  sd_ch(ich)%post_release_lock_timer = prm%post_release_lock_days
              endif
              sd_ch(ich)%ice_release_active = 0
              sd_ch(ich)%jam_inactive_days = 0
          endif
      endif

      ice_hydro_material = (sd_ch(ich)%ice_vol >= prm%jam_material_min_vol) .or. &
                           (sd_ch(ich)%ice_mobile >= prm%jam_mobile_min_vol) .or. &
                           (sd_ch(ich)%ice_mobile_pass >= prm%jam_mobile_min_vol)

      if (.not. sd_ch(ich)%is_releasing .and. sd_ch(ich)%release_active_today == 0) then
          post_release_flush = sd_ch(ich)%post_release_lock_timer > 0
          post_release_capture_eff = prm%post_release_capture_frac
          if (post_release_flush .and. sd_ch(ich)%major_release_done == 1) then
              post_release_capture_eff = prm%major_post_release_capture_frac
          endif

          select case (sd_ch(ich)%ice_phase)
          case (ICE_FREEZEUP)
              if (post_release_flush) then
                  wedge_capture_frac = post_release_capture_eff * prm%wedge_capture_jam_frac
              else if (sd_ch(ich)%is_jamming) then
                  if (aged_jam) then
                      wedge_capture_frac = 0.50 * prm%wedge_capture_jam_frac
                  else
                      wedge_capture_frac = prm%wedge_capture_jam_frac
                  endif
              else
                  wedge_capture_frac = 0.0
              endif
          case (ICE_STABLE)
              ! Stable-cover restriction is treated as event-scale storage routing,
              ! not as a spring breakup event.  On the rising limb, an intact cover
              ! can impound more water; as wedge storage builds, hydraulic head and
              ! through-flow relieve further capture.
              wedge_ratio = sd_ch(ich)%ice_wedge_stor / max(wedge_capacity, 1.0e-6)
              stable_storage_relief = 1.0 / (1.0 + 3.0 * max(0.0, wedge_ratio))
              stable_flow_boost = 1.0 + min(0.75, 2.0 * q_rise_pos)
              ! on the falling limb, high wedge storage should drain through
              ! the ice cover rather than continue as net impoundment.
              stable_capture_factor = stable_flow_boost * stable_storage_relief / &
                                      (1.0 + 10.0 * q_fall_pos)
              stable_capture_factor = icejam_clamp(stable_capture_factor, 0.05, 1.50)
              wedge_capture_frac = prm%wedge_capture_cover_frac * stable_capture_factor
              ! stable-cover protection must not be vetoed by a
              ! single-day drop in simulated thickness.  Structural integrity
              ! plus either sufficient thickness or seasonal maturity is enough
              ! to keep the intact-cover storage/leakage regime active.
              deepwinter_cover_ready = sd_ch(ich)%ice_phase == ICE_STABLE .and. &
                  sd_ch(ich)%ice_integrity >= prm%deepwinter_integrity_floor .and. &
                  (sim_ice_thick >= prm%stable_ice_thick .or. &
                   ice_maturity >= prm%major_ice_maturity_min)
              if (deepwinter_cover_ready) then
                  ! winter drainage is a non-event under-ice through-flow
                  ! response.  It requires both thermal/surface weakening and a
                  ! hydraulic driver.  A cold high-flow or snowfall day alone
                  ! should not start controlled winter drainage.
                  winter_thermal_ready = sd_ch(ich)%ice_surface_weak >= 0.35 .or. &
                                         t_air >= 0.0 .or. tmax >= 3.0 .or. &
                                         ros_water_mm >= prm%major_ros_min
                  winter_hydro_ready = q_rise_pos >= prm%major_qrise_min .or. &
                                       q_rel >= 0.20 .or. &
                                       wedge_ratio >= 0.30
                  winter_pulse_ready = winter_thermal_ready .and. winter_hydro_ready
                  if (winter_pulse_ready) then
                      sd_ch(ich)%winter_drain_timer = max(sd_ch(ich)%winter_drain_timer, &
                                                           prm%winter_pulse_drain_days)
                      winter_drain_set_today = .true.
                  endif

                  deepwinter_flow_factor = max(0.0, q_rel) / (max(0.0, q_rel) + prm%deepwinter_cover_q_ref)
                  deepwinter_cover_factor = ice_maturity * sd_ch(ich)%ice_integrity * &
                      (1.0 + prm%deepwinter_cover_flow_boost * deepwinter_flow_factor) / &
                      (1.0 + prm%deepwinter_cover_q_damp * max(0.0, q_rel))
                  wedge_capture_frac = max(wedge_capture_frac, &
                      prm%deepwinter_cover_capture_frac * deepwinter_cover_factor)
              endif
          case (ICE_BREAKUP)
              if (.not. breakup_tail_allowed) then
                  wedge_capture_frac = 0.0
              else if (post_release_flush) then
                  ! A jam-break release temporarily opens conveyance.  During
                  ! lockout, avoid immediately rebuilding a large tail impoundment.
                  wedge_capture_frac = post_release_capture_eff * prm%wedge_capture_tail_frac * &
                                       tail_material_factor
              else if (sd_ch(ich)%is_jamming) then
                  if (aged_jam) then
                      wedge_capture_frac = 0.50 * prm%wedge_capture_jam_frac
                  else
                      wedge_capture_frac = prm%wedge_capture_jam_frac
                  endif
              else if (ice_hydro_material .and. ice_mobile_excess > 0.0) then
                  ! Residual ice alone is not enough for tail impoundment.  Tail
                  ! capture requires a positive local mobile-ice transport deficit.
                  wedge_capture_frac = prm%wedge_capture_tail_frac * tail_material_factor
              else
                  wedge_capture_frac = 0.0
              endif
          case default
              wedge_capture_frac = 0.0
          end select

          if (wedge_capture_frac > 0.0 .and. ht1_adj > 1.0e-6) then
              ! Only an active ice jam receives the extra mobile-ice capture boost.
              ! Stable-cover and non-jam breakup-tail capture remain cover/tail
              ! hydraulics and are not amplified by residual mobile ice alone.
              if (sd_ch(ich)%is_jamming) then
                  if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
                      wedge_capture_frac = min(0.95, wedge_capture_frac + &
                                           prm%major_capture_boost * major_jam_factor)
                  endif
                  capture_mobile_term = prm%mobile_wedge_capture_weight * mobile_factor * &
                                        (1.0 - wedge_capture_frac)
              else
                  capture_mobile_term = 0.0
              endif
              capture_mult = min(0.95, reach%jam_block_modifier * &
                  (wedge_capture_frac + capture_mobile_term))
              q_damp_eff = prm%wedge_capture_q_damp
              if (sd_ch(ich)%ice_phase == ICE_STABLE .and. deepwinter_cover_ready) then
                  q_damp_eff = prm%wedge_capture_q_damp * prm%deepwinter_capture_q_damp_frac
              endif
              retention_frac = capture_mult * ice_maturity / &
                               (1.0 + q_damp_eff * max(0.0, q_rel))
              retention_frac = icejam_clamp(retention_frac, 0.0, 0.95)
              underice_alpha = prm%underice_alpha_min + &
                  (prm%underice_alpha_max - prm%underice_alpha_min) * sd_ch(ich)%ice_surface_weak
              underice_alpha = icejam_clamp(underice_alpha, prm%underice_alpha_min, prm%underice_alpha_max)
              underice_capacity = q_bankfull * underice_alpha * 86400.0
              underice_excess = max(0.0, ht1_adj - underice_capacity)
              background_capture = prm%wedge_base_capture_frac * ht1_adj
              excess_capture = retention_frac * underice_excess
              wedge_avail = max(0.0, wedge_capacity - sd_ch(ich)%ice_wedge_stor)
              wedge_capture = min(wedge_avail, background_capture + excess_capture, ht1_adj)

              ! under intact deep-winter cover, use the cover-retention concept.  Ice maturity and structural
              ! integrity reduce under-ice conveyance; only flow above that
              ! capacity is impounded, with an explicit deep-winter retention
              ! floor.  This avoids adding more tuning knobs while returning
              ! to the physical mechanism: thicker/stronger cover passes less
              ! water and stores more backwater.
              if (sd_ch(ich)%ice_phase == ICE_STABLE .and. deepwinter_cover_ready) then
                  deepwinter_cover_index = icejam_clamp(ice_maturity * &
                                            sd_ch(ich)%ice_integrity, 0.0, 1.0)
                  deep_underice_alpha = max(prm%underice_alpha_min, &
                      prm%underice_alpha_max * (1.0 - deepwinter_cover_index)**2)
                  underice_capacity = q_bankfull * deep_underice_alpha * 86400.0
                  underice_excess = max(0.0, ht1_adj - underice_capacity)
                  deep_retention_frac = prm%deepwinter_cover_capture_frac * &
                      deepwinter_cover_index * reach%jam_block_modifier / &
                      (1.0 + prm%wedge_capture_q_damp * &
                       prm%deepwinter_capture_q_damp_frac * max(0.0, q_rel))
                  deep_retention_frac = icejam_clamp(deep_retention_frac, 0.25, 0.90)
                  deep_block_capacity = max(0.0, 3.0 * sd_ch(ich)%ice_vol * &
                                            reach%jam_block_modifier)
                  deep_wedge_capture = min(wedge_avail, underice_excess, &
                                           deep_retention_frac * ht1_adj, &
                                           deep_block_capacity)
                  wedge_capture = max(wedge_capture, deep_wedge_capture)
              endif

              if (wedge_capture > 1.0e-6) then
                  ht1_adj = ht1_adj - wedge_capture
                  sd_ch(ich)%ice_wedge_stor = sd_ch(ich)%ice_wedge_stor + wedge_capture
              else
                  wedge_capture = 0.0
              endif
              sd_ch(ich)%ice_wedge_capture = wedge_capture
              sd_ch(ich)%icejam_block = wedge_capture
          endif
      endif

      if (.not. sd_ch(ich)%is_releasing .and. sd_ch(ich)%release_active_today == 0 .and. &
          sd_ch(ich)%ice_wedge_stor > 1.0e-6) then
          if (sd_ch(ich)%ice_phase == ICE_OPEN) then
              wedge_leak = min(sd_ch(ich)%ice_wedge_stor, prm%open_wedge_leak_frac * sd_ch(ich)%ice_wedge_stor)
          else
              ! keep non-event leakage below event-release scale.
              ! BREAKUP background drainage is weaker than generic daily leakage;
              ! post-release leakage distinguishes major vs ordinary releases.
              if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
                  leak_mult = prm%breakup_background_leak_mult
              else
                  leak_mult = 1.0
              endif
              if (sd_ch(ich)%post_release_lock_timer > 0) then
                  if (sd_ch(ich)%major_release_done == 1 .and. &
                      sd_ch(ich)%ice_phase == ICE_BREAKUP) then
                      leak_mult = max(leak_mult, prm%post_release_leak_mult)
                  else
                      leak_mult = max(leak_mult, prm%ordinary_post_release_leak_mult)
                  endif
              endif
              ! aged-jam leakage is a hydraulic seepage/through-flow
              ! enhancement, not a release surrogate.  It is only enabled
              ! when the aged jam is also hydraulically forced (F_R >= 1).
              if (sd_ch(ich)%is_jamming .and. sd_ch(ich)%jam_timer >= prm%max_jam_days .and. &
                  fr_ratio >= 1.0) then
                  leak_mult = max(leak_mult, prm%aged_jam_leak_mult)
              endif
              if (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. &
                  .not. post_release_flush .and. .not. sd_ch(ich)%is_jamming) then
                  leak_mult = leak_mult * max(0.25, tail_material_factor)
              endif
              if (sd_ch(ich)%ice_phase == ICE_STABLE) then
                  ! Stable ice-cover storage should behave as temporary under-ice
                  ! detention.  Under intact deep-winter cover, daily leakage is
                  ! deliberately suppressed so wedge storage can accumulate before
                  ! the spring mechanical breakup/warm-flush release.  Outside the
                  ! protected deep-winter state, high storage/head and falling limb
                  ! still enhance through-flow.
                  wedge_ratio = sd_ch(ich)%ice_wedge_stor / max(wedge_capacity, 1.0e-6)
                  if (deepwinter_cover_ready) then
                      stable_leak_mult = prm%deepwinter_leak_mult * &
                                         (1.0 + 0.50 * max(0.0, wedge_ratio) + 2.0 * q_fall_pos)
                      ! winter_drain is handled as capacity-limited
                      ! additional leakage below. Do not raise the whole stable
                      ! storage-based leakage multiplier to a pulse value.
                      stable_leak_mult = icejam_clamp(stable_leak_mult, prm%deepwinter_leak_mult, 1.0)
                      leak_mult = min(leak_mult, stable_leak_mult)
                  else
                      stable_leak_mult = 1.0 + 4.0 * max(0.0, wedge_ratio) + &
                                         0.10 * max(0.0, fr_ratio - 1.0) + 10.0 * q_fall_pos
                      stable_leak_mult = icejam_clamp(stable_leak_mult, &
                                         prm%deepwinter_leak_mult, &
                                         prm%stable_unprotected_leak_max_mult)
                      leak_mult = max(leak_mult, stable_leak_mult)
                  endif
              endif
              wedge_leak = min(sd_ch(ich)%ice_wedge_stor, &
                               leak_mult * prm%tail_wedge_leak_frac * sd_ch(ich)%ice_wedge_stor)
              if (sd_ch(ich)%ice_phase == ICE_STABLE .and. deepwinter_cover_ready .and. &
                  sd_ch(ich)%winter_drain_timer > 0) then
                  ! controlled winter drainage is limited by the extra
                  ! under-ice conveyance opened by surface weakening.  This is
                  ! scale-adaptive through q_bankfull and avoids an arbitrary
                  ! absolute leak cap.
                  deepwinter_cover_index = icejam_clamp(ice_maturity * &
                                            sd_ch(ich)%ice_integrity, 0.0, 1.0)
                  winter_alpha_intact = max(prm%underice_alpha_min, &
                      prm%underice_alpha_max * (1.0 - deepwinter_cover_index)**2)
                  ! once a winter drainage pulse is active, its
                  ! conveyance memory should not collapse solely because the
                  ! same-day surface-weakening diagnostic has returned to zero.
                  ! Use the timer itself as event memory with a modest minimum
                  ! effective weakening.
                  winter_alpha_weak = winter_alpha_intact + &
                      (prm%underice_alpha_max - winter_alpha_intact) * &
                      max(sd_ch(ich)%ice_surface_weak, 0.35)
                  winter_alpha_weak = icejam_clamp(winter_alpha_weak, &
                                      winter_alpha_intact, prm%underice_alpha_max)
                  winter_qcap_intact = q_bankfull * winter_alpha_intact
                  winter_qcap_weak = q_bankfull * winter_alpha_weak
                  ! controlled winter drainage is non-event drainage. It is
                  ! limited by the theoretical extra under-ice conveyance opened
                  ! by surface weakening and by only a fraction of the actual
                  ! inflow excess over intact-cover conveyance. Large drainage
                  ! should not be completed within STABLE as if it were a jam-break
                  ! release. The pulse factor ties the opportunity to the current
                  ! or recent thermal/melt signal while preserving timer memory.
                  winter_extra_qcap = max(0.0, winter_qcap_weak - winter_qcap_intact)
                  winter_actual_excess_q = max(0.0, q_in - winter_qcap_intact)
                  winter_pulse_factor = max(max(sd_ch(ich)%ice_surface_weak, 0.35), &
                                            max(warm_air_factor, meltwater_factor))
                  winter_pulse_factor = icejam_clamp(winter_pulse_factor, 0.0, 1.0)
                  winter_capacity_leak = min(winter_extra_qcap, &
                                         prm%winter_drain_excess_frac * winter_actual_excess_q) * &
                                         winter_pulse_factor * 86400.0
                  ! winter_drain is additional controlled under-ice
                  ! drainage, not a replacement for the whole stable leakage
                  ! formula. It is also capped as a small fraction of current
                  ! wedge storage so non-event drainage cannot function as a
                  ! jam-break release.
                  winter_capacity_leak = min(winter_capacity_leak, &
                                             prm%winter_drain_storage_frac * &
                                             sd_ch(ich)%ice_wedge_stor)
                  winter_additional_leak = min(max(0.0, sd_ch(ich)%ice_wedge_stor - wedge_leak), &
                                               winter_capacity_leak)
                  wedge_leak = wedge_leak + winter_additional_leak
              endif
          endif
          ht1_adj = ht1_adj + wedge_leak
          sd_ch(ich)%ice_wedge_stor = sd_ch(ich)%ice_wedge_stor - wedge_leak
          sd_ch(ich)%ice_wedge_leak = wedge_leak
      endif

      if (sd_ch(ich)%warm_flush_timer > 0 .and. .not. warm_flush_set_today) then
          sd_ch(ich)%warm_flush_timer = sd_ch(ich)%warm_flush_timer - 1
      endif
      if (sd_ch(ich)%winter_drain_timer > 0 .and. .not. winter_drain_set_today) then
          sd_ch(ich)%winter_drain_timer = sd_ch(ich)%winter_drain_timer - 1
      endif
      if (sd_ch(ich)%major_release_pending_timer > 0 .and. .not. major_release_gate) then
          sd_ch(ich)%major_release_pending_timer = sd_ch(ich)%major_release_pending_timer - 1
      endif
      if (sd_ch(ich)%major_bg_timer > 0 .and. .not. major_bg_set_today .and. &
          .not. major_release_gate) then
          sd_ch(ich)%major_bg_timer = sd_ch(ich)%major_bg_timer - 1
      endif

      sd_ch(ich)%ice_excess_storage = max(0.0, sd_ch(ich)%ice_wedge_stor - &
                                          prm%wedge_release_storage_frac * wedge_capacity)
      sd_ch(ich)%icejam_qadj = ht1_adj / 86400.0
      wedge_ratio = sd_ch(ich)%ice_wedge_stor / max(wedge_capacity, 1.0e-6)

      ! 9. Hydraulic routing directives for ch_rtmusk.
      tail_factor = min(1.0, (sd_ch(ich)%ice_vol + sd_ch(ich)%ice_mobile) / max(prm%reference_ice_vol, 1.0e-6))
      select case (sd_ch(ich)%ice_phase)
      case (ICE_FREEZEUP)
          if (sd_ch(ich)%is_releasing .or. sd_ch(ich)%release_active_today == 1) then
              sd_ch(ich)%ice_hydro_active = 1
              sd_ch(ich)%ice_k_mult = prm%k_release_mult
              sd_ch(ich)%ice_x_current = prm%x_release
          else if (sd_ch(ich)%post_release_lock_timer > 0) then
              sd_ch(ich)%ice_hydro_active = 0
              sd_ch(ich)%ice_k_mult = 1.0
              sd_ch(ich)%ice_x_current = 0.20
          else if (sd_ch(ich)%is_jamming) then
              sd_ch(ich)%ice_hydro_active = 1
              ! Freeze-up jams usually represent partial blockage and frazil/border-ice
              ! congestion.  Use a weaker routing modifier than breakup jams.
              sd_ch(ich)%ice_k_mult = 1.0 + 0.50 * (prm%k_jam_mult - 1.0)
              sd_ch(ich)%ice_x_current = 0.5 * (prm%x_jam + prm%x_cover)
          else if (ice_hydro_material) then
              sd_ch(ich)%ice_hydro_active = 1
              sd_ch(ich)%ice_k_mult = prm%k_cover_mult
              sd_ch(ich)%ice_x_current = prm%x_cover
          endif
      case (ICE_STABLE)
          ! Stable-cover restriction is already represented by wedge capture/leakage.
          ! disables dynamic Muskingum K/X changes in STABLE to avoid double
          ! counting ice-cover storage effects.
          sd_ch(ich)%ice_hydro_active = 0
          sd_ch(ich)%ice_k_mult = 1.0
          sd_ch(ich)%ice_x_current = 0.20
      case (ICE_BREAKUP)
          if ((sd_ch(ich)%is_releasing .or. sd_ch(ich)%release_active_today == 1) .and. &
              breakup_tail_allowed) then
              sd_ch(ich)%ice_hydro_active = 1
              sd_ch(ich)%ice_k_mult = prm%k_release_mult
              sd_ch(ich)%ice_x_current = prm%x_release
          else if (sd_ch(ich)%post_release_lock_timer > 0 .and. breakup_tail_allowed) then
              sd_ch(ich)%ice_hydro_active = 0
              sd_ch(ich)%ice_k_mult = 1.0
              sd_ch(ich)%ice_x_current = 0.20
          else if (sd_ch(ich)%is_jamming .and. breakup_tail_allowed) then
              sd_ch(ich)%ice_hydro_active = 1
              sd_ch(ich)%ice_k_mult = 1.0 + prm%breakup_jam_k_frac * (prm%k_jam_mult - 1.0)
              sd_ch(ich)%ice_x_current = prm%x_jam
          else if (ice_hydro_material .and. breakup_tail_allowed .and. ice_mobile_excess > 0.0) then
              sd_ch(ich)%ice_hydro_active = 1
              sd_ch(ich)%ice_k_mult = 1.0 + prm%k_tail_max * tail_factor * tail_material_factor
              sd_ch(ich)%ice_x_current = 0.20
          endif
      case default
          ! OPEN phase has no ice-related routing modifier.
      end select
      sd_ch(ich)%ice_k_mult = icejam_clamp(sd_ch(ich)%ice_k_mult, prm%k_min_mult, prm%k_max_mult)
      sd_ch(ich)%ice_x_current = icejam_clamp(sd_ch(ich)%ice_x_current, 0.0, 0.49)

      ! Optional switch: keep all icejam volume adjustments but prevent dynamic
      ! Muskingum K/X modification in ch_rtmusk.
      if (prm%icejam_msk_dynamic == 0) then
          sd_ch(ich)%ice_hydro_active = 0
          sd_ch(ich)%ice_k_mult = 1.0
          sd_ch(ich)%ice_x_current = 0.20
      endif

      ! Optional diagnostic trace for outlet channel 68.
      if (ich == 68) then
          write(9003,*) time%yrc, time%day, ich, &
              "phase", sd_ch(ich)%ice_phase, &
              "jamming", sd_ch(ich)%is_jamming, &
              "releasing", sd_ch(ich)%is_releasing, &
              "jam_today", sd_ch(ich)%jam_active_today, &
              "rel_today", sd_ch(ich)%release_active_today, &
              "phase_days", sd_ch(ich)%ice_phase_days, &
              "tail_allowed", breakup_tail_allowed, &
              "jam_timer", sd_ch(ich)%jam_timer, &
              "lock", sd_ch(ich)%post_release_lock_timer, &
              "tmax", tmax, "tave", t_air, &
              "freeze_dd", sd_ch(ich)%ice_freeze_dd, &
              "thaw_dd", sd_ch(ich)%ice_thaw_dd, &
              "thaw_weak", thaw_weak, &
              "Iice", sd_ch(ich)%ice_integrity, &
              "Iice_pk", sd_ch(ich)%ice_integrity_peak, &
              "Iweak", sd_ch(ich)%ice_surface_weak, &
              "Isurf", sd_ch(ich)%ice_surface_int, &
              "snowpack", sd_ch(ich)%snowpack, &
              "snow_ante", sd_ch(ich)%snowpack_ante, &
              "snow_peak", sd_ch(ich)%snowpack_peak, &
              "frz_avg", sd_ch(ich)%frz_surf_avg, &
              "frz_frac", sd_ch(ich)%frz_area_frac, &
              "warm_flush", sd_ch(ich)%warm_flush_today, &
              "flush_timer", sd_ch(ich)%warm_flush_timer, &
              "winter_drain", sd_ch(ich)%winter_drain_timer, &
              "major_fac", sd_ch(ich)%major_jam_factor, &
              "major_rel", sd_ch(ich)%major_release_today, &
              "major_done", sd_ch(ich)%major_release_done, &
              "major_pending", sd_ch(ich)%major_release_pending_timer, &
              "major_bg_timer", sd_ch(ich)%major_bg_timer, &
              "major_stor", merge(1, 0, major_storage_ready), &
              "wedge_ratio_major", wedge_ratio_for_major, &
              "snow_fac", snowpack_factor, &
              "warm_fac", warm_flush_factor, &
              "warm_mem", warm_memory_factor, &
              "base_fac", major_base_factor, &
              "bg_fac", major_background_factor, &
              "trig_fac", major_trigger_factor, &
              "deep_fac", deep_winter_factor, &
              "ice_fac", ice_storage_factor, &
              "snow_ante_fac", snowpack_ante_factor, &
              "snow_peak_fac", snowpack_peak_factor, &
              "frz_fac", frozen_soil_factor, &
              "warm_air_fac", warm_air_factor, &
              "melt_fac", meltwater_factor, &
              "fr_fac", fr_factor, &
              "qrise_fac", qrise_factor, &
              "runoff_fac", runoff_response_factor, &
              "hydro_fac", discharge_factor, &
              "ice_vol", sd_ch(ich)%ice_vol, &
              "ice_mobile", sd_ch(ich)%ice_mobile, &
              "ice_pass", sd_ch(ich)%ice_mobile_pass, &
              "sim_thick", sim_ice_thick, &
              "q_in", q_in, "q_bnk", q_bankfull, &
              "force", force_F, "resist", resistance_R, "F_R", fr_ratio, &
              "susc", reach%jam_susc, &
              "trans_cap", ice_transport_cap, &
              "mob_excess", ice_mobile_excess, &
              "wedge", sd_ch(ich)%ice_wedge_stor, &
              "wedge_cap", wedge_capacity, &
              "wedge_ratio", wedge_ratio, &
              "wedge_cap_day", sd_ch(ich)%ice_wedge_capture, &
              "wedge_rel", sd_ch(ich)%ice_wedge_release, &
              "wedge_leak", sd_ch(ich)%ice_wedge_leak, &
              "ht1_raw", ht1_raw, "ht1_adj", ht1_adj, &
              "Kmult", sd_ch(ich)%ice_k_mult, &
              "Xice", sd_ch(ich)%ice_x_current
      endif

      ht1%flo = max(0.0, ht1_adj)

      return
end subroutine sd_channel_icejam

subroutine sd_channel_ice_advect(j)

!!    Advect pass-through mobile ice to downstream chandeg objects.
!!    Mobile ice is not liquid water and this routine does not modify ht1/ht2.

      use hydrograph_module
      use sd_channel_module
      use channel_module
      use sd_channel_icejam_module

      implicit none

      integer, intent(in) :: j

      integer, parameter :: ICE_OPEN     = 0
      integer, parameter :: ICE_FREEZEUP = 1
      integer, parameter :: ICE_STABLE   = 2
      integer, parameter :: ICE_BREAKUP  = 3

      type(icejam_param_type), save :: prm
      type(icejam_reach_scale_type) :: reach_dn
      logical, save :: prm_initialized = .false.
      real, parameter :: ice_eps = 1.0e-4

      integer :: iout, iob_dn, ich_dn
      real :: ice_out, ice_to_dn, ice_sent, ice_unsent
      real :: frac_dn, capture_frac, capture_capacity
      real :: ice_capture, ice_pass
      real :: sim_ice_thick_dn, ice_ratio_dn, ice_depth_ratio_dn
      logical :: has_downstream_channel, has_downstream_object
      logical :: source_can_export, downstream_can_store

      if (j <= 0) return
      if (sd_ch(j)%ice_mobile_pass <= ice_eps) then
          sd_ch(j)%ice_mobile_pass = 0.0
          return
      endif

      if (.not. prm_initialized) then
          call icejam_default_params(prm)
          call icejam_validate_params(prm)
          prm_initialized = .true.
      endif

      ! Do not allow residual mobile ice to be exported outside FREEZEUP/BREAKUP.
      ! This prevents a single upstream reach with stale mobile ice from relocking
      ! the downstream river network after the breakup season has ended.
      source_can_export = sd_ch(j)%ice_phase == ICE_FREEZEUP .or. &
                          sd_ch(j)%ice_phase == ICE_BREAKUP
      if (.not. source_can_export) then
          sd_ch(j)%ice_mobile_pass = 0.0
          return
      endif

      ice_out = sd_ch(j)%ice_mobile_pass
      ice_sent = 0.0
      has_downstream_channel = .false.
      has_downstream_object = .false.

      do iout = 1, ob(icmd)%src_tot
          iob_dn = ob(icmd)%obj_out(iout)
          if (iob_dn <= 0) cycle
          has_downstream_object = .true.
          frac_dn = max(0.0, min(1.0, ob(icmd)%frac_out(iout)))
          if (frac_dn <= ice_eps) cycle

          if (trim(ob(iob_dn)%typ) == "chandeg") then
              has_downstream_channel = .true.
              ich_dn = ob(iob_dn)%num
              if (ich_dn <= 0) cycle

              ice_to_dn = min(frac_dn * ice_out, max(0.0, ice_out - ice_sent))
              if (ice_to_dn <= ice_eps) cycle

              downstream_can_store = sd_ch(ich_dn)%ice_phase == ICE_FREEZEUP .or. &
                                     sd_ch(ich_dn)%ice_phase == ICE_BREAKUP
              if (.not. downstream_can_store) then
                  ! The downstream reach is not in FREEZEUP/BREAKUP.  The incoming mobile ice
                  ! is treated as hydraulically ineffective and is removed from
                  ! the mobile-ice routing chain rather than stored as new ice
                  ! material that could prevent OPEN from persisting.
                  ice_sent = ice_sent + ice_to_dn
                  cycle
              endif

              call icejam_compute_reach_scale(prm, sd_ch(ich_dn)%chw, sd_ch(ich_dn)%chl, &
                  sd_ch(ich_dn)%chd, sd_ch(ich_dn)%chs, sd_ch(ich_dn)%sinu, &
                  ch_rcurv(ich_dn)%elev(1)%flo_rate, ch_rcurv(ich_dn)%elev(2)%flo_rate, reach_dn)

              sim_ice_thick_dn = sd_ch(ich_dn)%ice / max(reach_dn%ice_area, 1.0e-6)
              ice_ratio_dn = icejam_clamp(sim_ice_thick_dn / max(prm%ice_maturity_ref_thick, 1.0e-6), 0.0, 1.0)
              ice_depth_ratio_dn = icejam_clamp(sim_ice_thick_dn / max(sd_ch(ich_dn)%chd, 1.0e-6), 0.0, 1.0)

              capture_frac = prm%mobile_capture_base + &
                  prm%mobile_capture_susc_weight * reach_dn%ice_capture_modifier + &
                  prm%mobile_capture_ice_weight * ice_ratio_dn + &
                  prm%mobile_capture_depth_weight * ice_depth_ratio_dn

              select case (sd_ch(ich_dn)%ice_phase)
              case (ICE_FREEZEUP)
                  capture_frac = max(capture_frac, prm%freezeup_capture_min)
              case (ICE_STABLE)
                  capture_frac = max(capture_frac, prm%stable_capture_min)
              case (ICE_BREAKUP)
                  capture_frac = max(capture_frac, prm%breakup_capture_min)
              case default
                  if (ice_ratio_dn <= 0.05) capture_frac = min(capture_frac, prm%warm_capture_max)
              end select
              capture_frac = icejam_clamp(capture_frac, prm%mobile_capture_min, prm%mobile_capture_max)

              capture_capacity = max(0.0, prm%mobile_capture_capacity_mult * reach_dn%ice_cap_vol)
              ice_capture = min(capture_frac * ice_to_dn, capture_capacity, ice_to_dn)
              ice_pass = max(0.0, ice_to_dn - ice_capture)

              sd_ch(ich_dn)%ice_mobile = sd_ch(ich_dn)%ice_mobile + ice_capture
              sd_ch(ich_dn)%ice_mobile_pass = sd_ch(ich_dn)%ice_mobile_pass + ice_pass
              ice_sent = ice_sent + ice_to_dn
          endif

          if (ice_sent >= ice_out - ice_eps) exit
      end do

      ice_unsent = max(0.0, ice_out - ice_sent)
      if (has_downstream_channel) then
          sd_ch(j)%ice_mobile_pass = ice_unsent
      else
          if (has_downstream_object) then
              sd_ch(j)%ice_mobile_pass = ice_unsent
          else
              sd_ch(j)%ice_mobile_pass = 0.0
          endif
      endif

      return
end subroutine sd_channel_ice_advect
