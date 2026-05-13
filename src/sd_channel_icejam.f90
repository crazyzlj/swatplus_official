subroutine sd_channel_icejam(j)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Conceptual daily river-ice and ice-jam storage-release module.
!!
!!    Annual process separation:
!!      ICE_WARM       : no active river-ice processes.
!!      ICE_FREEZEUP   : ice-cover growth and possible weak freeze-up obstruction.
!!      ICE_DEEPWINTER : stable ice-cover retention; no ordinary release.
!!      ICE_BREAKUP    : breakup release, mobile ice generation, and new jams.
!!
!!    Trigger logic:
!!      thermal_trigger    : accumulated thaw forcing weakens ice.
!!      mechanical_trigger : local + mobile ice load exceeds support capacity.
!!
!!    Water/ice balance:
!!      ice growth       : ch_stor/ht1%flo -> sd_ch%ice
!!      ice melt         : sd_ch%ice -> ht1%flo
!!      winter retention : ht1%flo -> sd_ch%ice_cover_stor
!!      jam block        : ht1%flo -> sd_ch%ice_jam_stor
!!      storage release  : source-specific release from ice_cover_stor and/or ice_jam_stor -> ht1%flo
!!      ice mobilization : breakup and drift convert sd_ch%ice -> sd_ch%ice_mobile
!!
!!    Required channel state variables:
!!      real    :: ice
!!      real    :: ice_cover_stor  !stable ice-cover retention storage
!!      real    :: ice_jam_stor    !freeze-up/breakup jam-blocking storage
!!      real    :: ice_mobile
!!      real    :: q_prev
!!      real    :: ice_thaw_dd
!!      real    :: ice_freeze_dd
!!      integer :: ice_phase
!!      integer :: ice_phase_days
!!      integer :: ice_release_active  !0 none, nonzero active release episode
!!      integer :: ice_release_days    !duration of current release episode in days
!!      integer :: ice_block_days      !duration of current blocking/buildup episode in days
!!
!!    Recommended placement:
!!      after gwflow-channel exchange
!!      before chsd_d/ch_in_d printing, sediment routing, and ch_rtmusk

      use basin_module
      use time_module
      use hydrograph_module
      use sd_channel_module
      use climate_module
      use hru_module, only : hru
      use sd_channel_icejam_module

      implicit none

      integer, intent(in) :: j

      integer, parameter :: JAM_NONE  = 0
      integer, parameter :: JAM_RELEASE = 1

      integer, parameter :: ICE_WARM       = 0
      integer, parameter :: ICE_FREEZEUP   = 1
      integer, parameter :: ICE_DEEPWINTER = 2
      integer, parameter :: ICE_BREAKUP    = 3

      integer, parameter :: BRK_DAY_NONE      = 0
      integer, parameter :: BRK_DAY_DRIFT     = 1
      integer, parameter :: BRK_DAY_BLOCK     = 2
      integer, parameter :: BRK_DAY_RELEASE   = 3

      integer :: ord = 1               !none    |stream/channel order
      integer :: old_phase = 0
      integer :: new_phase = 0

      real :: ch_vol_cap = 0.          !m3      |approximate bankfull channel volume
      real :: ice_cover_max = 0.       !m3      |alias of ice_cap_vol for backward-compatible calculations
      real :: ice_ratio = 0.           !none    |ice_maturity = mean ice thickness / ice_maturity_ref_thick
      real :: mobile_jam_ratio = 0.    !none    |ice_mobile / ice_cover_max
      real :: mobile_pass_jam_ratio = 0. !none |ice_mobile_pass / ice_cover_max
      real :: ice_stor_eps = 0.        !m3      |negligible ice-cover storage threshold
      
      real :: q_in_rate_raw = 0.       !m3/s    |raw inflow rate before ice-jam adjustment
      real :: q_jam_ref_rate = 0.      !m3/s    |ice-jam trigger reference flow rate
      real :: q_ratio = 0.             !none    |q_in_rate_raw / q_jam_ref_rate

      real :: tw_ice = 0.              !deg C   |water temperature proxy for ice processes
      real :: t_air = 0.               !deg C   |daily mean air temperature
      real :: tmax = 0.                !deg C   |daily maximum air temperature
      real :: t_ice_growth = 0.        !deg C   |temperature driver for ice growth
      real :: t_ice_decay = 0.         !deg C   |temperature driver for ice decay / breakup
      real :: t_freeze = 0.            !deg C   |daily freeze forcing for freeze-memory index
      real :: t_thaw = 0.              !deg C   |daily thaw forcing for breakup-memory index
      real :: thaw_tmax_base_eff = 0.  !deg C   |effective Tmax base for thaw forcing
      real :: snow_melt_mm = 0.        !mm      |channel-scale HRU-aggregated snowmelt used for ROS filtering
      
      real :: ice_growth = 0.          !m3/day  |actual ice-cover growth
      real :: ice_growth_pot = 0.      !m3/day  |potential ice-cover growth
      real :: ice_growth_cap = 0.      !m3/day  |daily cap on ice-cover growth
      real :: ice_target_thick = 0.    !m       |Stefan-type target ice thickness
      real :: ice_growth_pot_thick = 0.!m       |potential growth before cap
      real :: ice_growth_cap_thick = 0.!m/day   |daily growth cap
      real :: ice_growth_thick = 0.    !m/day   |actual growth thickness
      real :: ice_target = 0.          !m3      |target ice storage from freezing-degree index
      real :: ice_decay = 0.           !m3/day  |ice-cover melt/decay
      real :: mobile_ice_decay = 0.    !m3/day  |mobile ice melt/decay
      real :: mobile_pass_decay = 0.   !m3/day  |pass-through mobile ice melt/flush
      real :: ice_melt_thick = 0.       !m/day   |thickness-based cover-ice melt
      real :: mobile_melt_thick = 0.    !m/day   |thickness-based mobile-ice melt equivalent

      
      real :: ice_avail = 0.           !m3/day  |liquid water available for ice-cover formation
      real :: freeze_from_chstor = 0.  !m3/day  |ice growth supplied from channel storage
      real :: freeze_from_ht1 = 0.     !m3/day  |ice growth supplied from incoming flow
      real :: freeze_remain = 0.       !m3/day  |remaining ice growth demand after using channel storage

      real :: jam_susc = 0.            !none    |channel ice-jam susceptibility
      real :: q_underice_cap = 0.      !m3/s    |effective under-ice conveyance capacity
      real :: underice_excess = 0.     !m3/day  |water exceeding under-ice conveyance
      real :: block_capacity = 0.      !m3/day  |daily ice-limited blocking capacity
      real :: block_cap_coeff = 0.     !none    |water-storage capacity per ice volume
      real :: block_frac_max = 0.      !none    |daily upper fraction of incoming flow
      real :: jam_stor_max_frac_eff = 0. !none |effective maximum jam storage fraction for BREAKUP jam reservoir
      real :: jam_stor_max = 0.        !m3      |BREAKUP ice-jam virtual reservoir capacity
      real :: jam_remaining_capacity_step = 0. !m3 |remaining jam capacity available for today's cover_to_jam/block_jam
      real :: jam_stor_eps = 0.        !m3      |negligible jam storage threshold
      real :: active_release_min = 0.  !m3      |minimum release to suppress same-day new jam
      real :: cover_stor_max = 0.      !m3      |stable ice-cover retention storage capacity
      real :: cover_stor_ratio = 0.    !none    |ice_cover_stor / cover_stor_max
      real :: jam_stor_ratio = 0.      !none    |unified ice_jam_stor / jam_stor_max
      real :: total_ice_stor = 0.      !m3      |ice_cover_stor + ice_jam_stor
      real :: cover_stor_cap = 0.    !m3      |ice-supported capacity for cover storage
      real :: jam_stor_cap = 0.      !m3      |ice-load-supported capacity for jam storage
      real :: cover_overflow_release = 0. !m3/day |storage released because ice support is insufficient
      real :: jam_overflow_release = 0.   !m3/day |jam storage released because ice support is insufficient

      real :: cover_q_cap = 0.
      real :: cover_underice_excess = 0.
      real :: cover_block_capacity = 0.
      real :: cover_stor_max_dbg = 0.
      real :: cover_remaining_capacity = 0.

      real :: jam_block_capacity = 0.
      real :: ice_load_block_capacity = 0.
      real :: jam_constriction_capacity = 0.
      real :: jam_maturity_factor = 0.
      real :: jam_presence_factor = 0.  !none |jam storage presence scaled by hydraulic storage, independent of jam_stor_max
      real :: jam_material_factor = 0.  !none |available ice material support for mature jam blocking
      real :: flow_supply_factor = 0.  !none |qraw / q_jam_ref_rate factor for jam blocking capacity
      real :: ice_block_material = 0.   !m3   |mobile ice material available for BREAKUP blocking
      real :: onset_block_mult = 0.     !none |derived weak residual cover obstruction factor on BREAKUP onset
      real :: jam_stor_max_dbg = 0.
      real :: jam_remaining_capacity = 0.

      real :: blocked = 0.             !m3/day  |water blocked into ice-jam storage
      real :: blocked_cover = 0.       !m3/day  |water retained by stable ice cover
      real :: blocked_jam = 0.         !m3/day  |water blocked by freeze-up/breakup ice jam
      real :: blocked_total = 0.       !m3/day  |total daily ice-related blocking
      real :: released = 0.            !m3/day  |water released from ice-jam storage
      real :: released_event = 0.      !m3/day  |event-based release for diagnostics
      real :: released_jam_event = 0.  !m3/day  |event release from ice-jam storage
      real :: released_leak = 0.       !m3/day  |thaw-only leakage from ice-related storage
      real :: released_cover_leak = 0. !m3/day  |leakage from stable cover storage
      real :: released_jam_leak = 0.   !m3/day  |leakage from ice-jam storage
      real :: stor_before = 0.         !m3      |total ice-related storage before release
      real :: release_ratio = 0.       !none    |fraction of storage released by event
      real :: ice_mobilized = 0.       !m3      |total mobilized ice that can be transported downstream
      real :: ice_mobilized_drift = 0. !m3      |background drift mobilization before event release
      real :: ice_mobilized_dynamic = 0. !m3    |dynamic/mechanical mobilization before event release
      real :: ice_mobilized_pre = 0.   !m3      |pre-event mobile ice generation
      real :: ice_mobilized_event = 0. !m3      |event-breakup mobile ice generation
      real :: ice_mobilized_cover_break = 0. !m3 |cover ice mobilized during cover_to_jam / jam buildup
      real :: ice_mobilized_total = 0. !m3      |pre-event plus event mobilized ice
      real :: local_mobile_capture_frac = 1.0 !none |fraction of newly mobilized local ice retained as local jam material
      real :: local_mobile_pass_frac = 0.0    !none |fraction of newly mobilized local ice routed as pass-through next day
      real :: ice_mobile_generated_pass = 0.0 !m3   |new local mobile ice assigned to pass-through pool
      real :: mobile_flushed_by_release = 0.0 !m3 |existing mobile ice flushed to pass-through pool by jam release
      real :: mobile_flush_frac = 0.0        !none |fraction of existing mobile ice flushed by jam release
      real :: mobile_drift_pass = 0.0        !m3   |existing mobile ice drifting downstream during DRIFT days
      real :: mobile_drift_pass_frac = 0.0   !none |daily pass-through fraction for existing mobile ice during DRIFT
      real :: warm_cleanup_return = 0.0      !m3   |mass-balance return from BREAKUP->WARM cleanup
      real :: drift_frac_eff = 0.      !none    |effective daily drift mobilization fraction
      real :: dynamic_frac_eff = 0.    !none    |effective daily dynamic mobilization fraction
      real :: ice_load = 0.            !m3      |local ice plus mobile ice load
      real :: ice_support_capacity = 0. !m3     |ice-load threshold for mechanical trigger
      real :: ice_strength_factor = 1.  !none   |phase-dependent ice support factor
      
      real :: raw_flo = 0.             !m3/day  |raw inflow volume before ice-jam adjustment
      real :: adj_ratio = 1.           !none    |ratio for subdaily tsin adjustment
      real :: tsin_sum = 0.            !m3/day  |sum of subdaily inflow hydrograph
      real :: retention_frac_eff = 0.  !none    |effective stable ice-cover retention fraction

      logical :: ros_day = .false.     !logical |channel-scale HRU-diagnosed rain-on-snow/ice flag after snowmelt filter
      logical :: breakup_thermal_ready = .false.
      logical :: breakup_material_ready = .false.
      logical :: deepwinter_age_ready = .false.
      logical :: freezeup_max_ready = .false.
      logical :: deepwinter_max_ready = .false.
      logical :: breakup_max_ready = .false.
      logical :: warm_to_freezeup_ready = .false.
      logical :: breakup_long_enough = .false.
      logical :: thermal_warm_ready = .false.
      logical :: ice_small_enough = .false.
      logical :: storage_not_active = .false.
      logical :: drift_weak_zone = .false. !true when thaw_weak lies between stable and release thresholds
      logical :: thermal_trigger = .false.
      logical :: mechanical_trigger = .false.
      logical :: jam_formation_ready = .false.
      logical :: block_flow_ready = .false. !true when qraw supply is sufficient for mature BREAKUP jam blocking
      logical :: allow_new_jam_today = .false.
      logical :: do_jam_formation_today = .false. !true only when mature ice-jam blocking should be computed
      logical :: do_cover_to_jam_today = .false.   !true when cover storage can be reclassified as jam storage
      logical :: do_onset_cover_block_today = .false. !true for weak residual cover-controlled blocking on BREAKUP onset
      logical :: phase_changed_today = .false.
      logical :: breakup_onset_today = .false.
      logical :: seasonal_breakup_reset = .false.
      logical :: breakup_release_gate = .false. !true when thaw weakening is sufficient for BREAKUP jam-storage release
      integer :: breakup_day_type = BRK_DAY_NONE !0 none, 1 drift/open-flow, 2 block/buildup, 3 release/recession

      real :: cover_to_jam = 0.      !m3/day  |cover-controlled storage reclassified as jam-controlled storage
      real :: cover_to_jam_frac_eff = 0. !none |effective cover-to-jam transfer fraction
      real :: cover_to_jam_capacity = 0. !m3   |remaining jam-storage capacity after new blocking
      real :: mobile_order_mult = 1.0 !none   |channel-order multiplier for cover-to-mobile ice conversion

      real :: cover_stor_before = 0. !m3      |cover storage before event release
      real :: jam_stor_before = 0.   !m3      |jam storage before event release
      real :: jam_stor_start = 0.    !m3      |jam storage at start of the daily release/block sequence
      real :: jam_remain_capacity_start = 0. !m3 |jam capacity available at daily start; not refilled by same-day release
      real :: jam_capacity_used_today = 0. !m3 |cover_to_jam + block_jam constrained by day-start capacity
      real :: jam_release_ratio = 0.   !none  |event release fraction from jam storage
      real :: mobilization_ratio = 0.  !none  |fraction of local ice mobilized
      real :: recession_frac = 0.      !none  |active-release recession fraction for jam storage
      real :: release_frac_eff = 0.  !none  |effective continuous jam-storage release fraction
      real :: release_weak_eff = 0.  !none |nonlinear thaw-weakening release response
      real :: release_ramp_factor = 1.0 !none |episode-day release ramp multiplier
      logical :: release_recession_day = .false. !true when ice_release_active continues after trigger day
      type(icejam_param_type), save :: prm
      type(icejam_reach_scale_type) :: reach
      logical, save :: prm_initialized = .false.
      real :: thaw_weakening_index = 0. !none |thaw_dd / (freeze_dd + thaw_dd)
      real :: ice_area = 1.          !m2      |channel water-surface area for ice-thickness scaling
      real :: ice_cap_vol = 1.       !m3      |characteristic ice volume = ice_maturity_ref_thick * ice_area
      real :: sim_ice_thick = 0.     !m       |simulated reach-average ice thickness
      real :: ice_maturity = 0.      !none    |sim_ice_thick / ice_maturity_ref_thick
      real :: phase_ret_mult = 1.0    !none    |phase-dependent multiplier for stable-cover retention
      real :: ice_depth_ratio = 0.    !none    |sim_ice_thick / channel depth, hydraulic obstruction diagnostic
      real :: ice_state_eps = 1.0      !m3      |numerical threshold for local cover-ice state
      real :: mobile_state_eps = 1.0   !m3      |numerical threshold for mobile-ice state
      real :: stor_state_eps = 1.0     !m3      |numerical threshold for ice-related water storage
      real :: warm_storage_exit_threshold = 1.0 !m3 |relative storage threshold for BREAKUP -> WARM
      logical :: ice_absent = .false.
      logical :: mobile_absent = .false.


      if (.not. prm_initialized) then
          call icejam_default_params(prm)
          call icejam_validate_params(prm)
          prm_initialized = .true.
      endif

      ich = j
      iwst = ob(icmd)%wst
      t_air = wst(iwst)%weat%tave
      tmax = wst(iwst)%weat%tmax
      ros_day = .false.


      !! ice_jam_flag is a daily diagnostic event flag and is reset every day.
      !! The cross-day release-episode state is stored in ice_release_active.
      sd_ch(ich)%ice_jam_flag = JAM_NONE
      sd_ch(ich)%icejam_block = 0.
      sd_ch(ich)%icejam_release = 0.
      sd_ch(ich)%icejam_qraw = 0.
      sd_ch(ich)%icejam_qadj = 0.
      sd_ch(ich)%icejam_qratio = 0.
      sd_ch(ich)%icejam_qrise = 0.
      sd_ch(ich)%icejam_susc = 0.

      blocked = 0.
      blocked_cover = 0.
      blocked_jam = 0.
      blocked_total = 0.
      released = 0.
      released_event = 0.
      released_jam_event = 0.
      released_leak = 0.
      released_cover_leak = 0.
      released_jam_leak = 0.
      release_ratio = 0.
      cover_stor_ratio = 0.
      jam_stor_ratio = 0.
      total_ice_stor = 0.
      cover_stor_cap = 0.
      jam_stor_cap = 0.
      cover_overflow_release = 0.
      jam_overflow_release = 0.
      mobile_jam_ratio = 0.
      ice_mobilized = 0.
      ice_mobilized_drift = 0.
      ice_mobilized_dynamic = 0.
      ice_mobilized_pre = 0.
      ice_mobilized_event = 0.
      ice_mobilized_cover_break = 0.
      ice_mobilized_total = 0.
      local_mobile_capture_frac = 1.0
      local_mobile_pass_frac = 0.
      ice_mobile_generated_pass = 0.
      mobile_flushed_by_release = 0.
      mobile_flush_frac = 0.
      mobile_drift_pass = 0.
      mobile_drift_pass_frac = 0.
      warm_cleanup_return = 0.
      drift_frac_eff = 0.
      dynamic_frac_eff = 0.
      block_capacity = 0.
      block_cap_coeff = 0.
      block_frac_max = 0.
      jam_remaining_capacity_step = 0.

      cover_q_cap = 0.
      cover_underice_excess = 0.
      cover_block_capacity = 0.
      cover_stor_max_dbg = 0.
      cover_remaining_capacity = 0.
      jam_block_capacity = 0.
      jam_stor_max_dbg = 0.
      jam_remaining_capacity = 0.

      allow_new_jam_today = .false.
      do_jam_formation_today = .false.
      do_cover_to_jam_today = .false.
      phase_changed_today = .false.
      breakup_onset_today = .false.
      seasonal_breakup_reset = .false.
      ice_absent = .false.
      mobile_absent = .false.
      breakup_day_type = BRK_DAY_NONE
      release_recession_day = .false.
      cover_to_jam = 0.
      cover_to_jam_frac_eff = 0.
      cover_to_jam_capacity = 0.
      mobile_order_mult = 1.0
      cover_stor_before = 0.
      jam_stor_before = 0.
      jam_stor_start = 0.
      jam_remain_capacity_start = 0.
      jam_capacity_used_today = 0.
      jam_release_ratio = 0.
      mobilization_ratio = 0.
      q_underice_cap = 0.
      underice_excess = 0.
      thaw_weakening_index = 0.
      ice_growth_cap = 0.
      ice_target_thick = 0.
      ice_growth_pot_thick = 0.
      ice_growth_cap_thick = 0.
      ice_growth_thick = 0.
      ice_melt_thick = 0.
      mobile_melt_thick = 0.
      mobile_pass_decay = 0.
      retention_frac_eff = 0.
      jam_stor_max_frac_eff = 0.
      jam_stor_max = 0.
      ice_load_block_capacity = 0.
      jam_constriction_capacity = 0.
      jam_maturity_factor = 0.
      jam_presence_factor = 0.
      jam_material_factor = 0.
      onset_block_mult = 0.
      jam_formation_ready = .false.
      block_flow_ready = .false.
      thermal_trigger = .false.
      mechanical_trigger = .false.

      !! ------------------------------------------------------------
      !! Reach-scale variables.
      !! The old stream-order select-case is replaced by a continuous
      !! pedotransfer-like susceptibility function.  Sinuosity is passed as
      !! 1.0 here because not all SWAT+ channel data structures expose it;
      !! This version uses the channel-specific sinuosity stored in sd_ch%sinu.
      !! ------------------------------------------------------------
      ord = sd_ch(ich)%order
      call icejam_compute_reach_scale(prm, sd_ch(ich)%chw, sd_ch(ich)%chl, &
              sd_ch(ich)%chd, sd_ch(ich)%chs, sd_ch(ich)%sinu, &
              ch_rcurv(ich)%elev(1)%flo_rate, ch_rcurv(ich)%elev(2)%flo_rate, reach)

      ch_vol_cap = reach%hyd_storage_scale
      ice_area = reach%ice_area
      ice_cap_vol = reach%ice_cap_vol
      ice_cover_max = ice_cap_vol
      q_jam_ref_rate = reach%q_jam_ref_rate
      jam_susc = reach%jam_susc
      mobile_order_mult = reach%mobile_order_mult
      sd_ch(ich)%icejam_susc = jam_susc

      !raw inflow before ice-jam modification.
      q_in_rate_raw = max(0., ht1%flo) / 86400.
      sd_ch(ich)%icejam_qraw = q_in_rate_raw

      !! v43: qraw is not an event trigger, but low flow supply should make
      !! mature BREAKUP jam blocking unlikely/weak.  Compute the flow-supply
      !! factor early so it can also be used in the BREAKUP day-type decision.
      flow_supply_factor = q_in_rate_raw / max(q_jam_ref_rate, 1.e-6)
      flow_supply_factor = max(0., min(1., flow_supply_factor))
      block_flow_ready = flow_supply_factor >= max(0.50, prm%jam_mobile_trigger_ratio)

      !Reach hydraulic capacity and ice-thickness capacity were computed above.
      !ch_vol_cap remains as a hydraulic storage scale for water impoundment;
      !ice_cover_max is now an alias for ice_cap_vol, not ch_vol_cap * ice_max_frac.

      !! Separate storage capacities and maturity ratios.
      !! ice_cover_stor is produced only by stable under-ice retention.
      !! ice_jam_stor is produced only by explicit freeze-up/breakup jam blocking.
      cover_stor_max = max(1., prm%ice_cover_ret_stor_frac * ch_vol_cap * reach%jam_storage_modifier)
      jam_stor_max = max(1., prm%jam_form_stor_max_frac * ch_vol_cap * reach%jam_storage_modifier)
      cover_stor_ratio = sd_ch(ich)%ice_cover_stor / cover_stor_max
      cover_stor_ratio = max(0., min(1., cover_stor_ratio))
      jam_stor_ratio = sd_ch(ich)%ice_jam_stor / max(jam_stor_max, 1.e-6)
      jam_stor_ratio = max(0., min(1., jam_stor_ratio))
      total_ice_stor = sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor

      jam_stor_eps = max(1.0, 1.e-6 * ch_vol_cap)
      active_release_min = max(jam_stor_eps, 0.01 * max(1.0, ht1%flo))

      ice_cover_max = max(ice_cap_vol, 1.e-6)
      ice_stor_eps = max(1.0, 1.e-6 * ice_cover_max)
      ice_state_eps = max(1.0, 1.e-6 * ice_cover_max)
      !! Mobile ice is a transport/jam-material pool; tiny numerical remnants
      !! should not keep mobile_jam_ratio nonzero during ice-free periods.
      mobile_state_eps = max(1.0, 1.e-4 * ice_cover_max)
      stor_state_eps = max(jam_stor_eps, 1.e-6 * ch_vol_cap)
      !! Relative storage threshold for ending BREAKUP. This is intentionally
      !! not a new calibration parameter: it is 1% of the hydraulic storage scale.
      warm_storage_exit_threshold = max(jam_stor_eps, 0.01 * ch_vol_cap)

      !Ice-jam trigger reference flow from centralized reach scale.

      q_ratio = q_in_rate_raw / q_jam_ref_rate
      sd_ch(ich)%icejam_qratio = q_ratio

      ! q-rise diagnostics and moving-average lag triggers were removed; q_prev stores raw inflow only.

      !thermal forcing and memory indices
      tw_ice = sd_ch(ich)%tmp_prx
      if (tw_ice < -20. .or. tw_ice > 40.) then
        tw_ice = t_air
      endif

      !ice growth responds to cold water/air conditions.
      t_ice_growth = min(tw_ice, t_air)

      !daily ice melt uses a combined water-air thermal proxy.
      t_ice_decay = 0.5 * tw_ice + 0.5 * t_air

      t_freeze = max(0., prm%ice_frz_tmp - t_ice_growth)
      sd_ch(ich)%ice_freeze_dd = prm%freeze_memory * sd_ch(ich)%ice_freeze_dd + t_freeze
      if (sd_ch(ich)%ice_freeze_dd < 1.e-6) sd_ch(ich)%ice_freeze_dd = 0.

      !! ------------------------------------------------------------------
      !! 1. Ice-cover growth during cold periods.
      !! Use a Stefan-type temperature-index relationship:
      !! target ice storage ratio increases with sqrt(accumulated freezing index).
      !! ------------------------------------------------------------------
      if (t_ice_growth < prm%ice_frz_tmp) then
          ice_target_thick = min(prm%ice_maturity_ref_thick, &
                  prm%ice_growth_coeff * sqrt(max(0., sd_ch(ich)%ice_freeze_dd)))
          ice_target = ice_target_thick * ice_area
          ice_growth_pot = max(0., ice_target - sd_ch(ich)%ice)
          ice_growth_pot_thick = ice_growth_pot / max(ice_area, 1.e-6)

          !! The Stefan-type target is an equilibrium/diagnostic target based
          !! on accumulated freezing memory.  Without a daily cap, the model
          !! can unrealistically "catch up" several centimeters in one day
          !! after intermittent freeze-thaw periods.  Limit growth by an
          !! explicit reach-mean daily thickness cap.
          ice_growth_cap_thick = max(0., prm%max_daily_ice_growth_thick)
          ice_growth_cap = ice_growth_cap_thick * ice_area
          ice_growth_pot = min(ice_growth_pot, ice_growth_cap)

          ice_avail = max(0., ch_stor(ich)%flo) + prm%ice_freeze_inflow_frac * max(0., ht1%flo)
          ice_growth = min(ice_growth_pot, ice_avail)
          ice_growth_thick = ice_growth / max(ice_area, 1.e-6)

          freeze_from_chstor = min(max(0., ch_stor(ich)%flo), ice_growth)
          ch_stor(ich)%flo = ch_stor(ich)%flo - freeze_from_chstor

          freeze_remain = ice_growth - freeze_from_chstor
          freeze_from_ht1 = min(prm%ice_freeze_inflow_frac * max(0., ht1%flo), freeze_remain)
          ht1%flo = ht1%flo - freeze_from_ht1

          sd_ch(ich)%ice = sd_ch(ich)%ice + freeze_from_chstor + freeze_from_ht1
      endif

      !! ------------------------------------------------------------------
      !! 2. Ice-cover melt / deterioration.
      !! Ice melt returns water to ht1%flo; event release uses thaw memory.
      !! ------------------------------------------------------------------
      if (t_ice_decay > prm%ice_melt_tmp .and. sd_ch(ich)%ice > 0.) then
          !! Thickness-based melt: reduce reach-mean ice thickness by a daily
          !! degree-day melt depth rather than by a fraction of current ice volume.
          !! This avoids the nonphysical behavior that thick ice melts faster only
          !! because more ice volume is present.
          ice_melt_thick = prm%ice_decay_coeff * (t_ice_decay - prm%ice_melt_tmp)
          ice_melt_thick = max(0., ice_melt_thick)
          ice_decay = min(sd_ch(ich)%ice, ice_melt_thick * ice_area)

          sd_ch(ich)%ice = sd_ch(ich)%ice - ice_decay
          ht1%flo = ht1%flo + ice_decay
      endif

      sd_ch(ich)%ice = max(0., sd_ch(ich)%ice)
      ice_stor_eps = max(1.0, 1.e-6 * ice_cover_max)

      if (sd_ch(ich)%ice <= ice_stor_eps .and. t_air > prm%warm_flush_tmp) then
          ht1%flo = ht1%flo + sd_ch(ich)%ice
          sd_ch(ich)%ice = 0.
      endif

      !the capacity limit is ice_cover_max locally, but ice transported from upstream can
      !make sd_ch(ich)%ice greater than ice_cover_max, this is reasonable and necessary.
      !sd_ch(ich)%ice = max(0., min(sd_ch(ich)%ice, ice_cover_max))
      sd_ch(ich)%ice = max(0., sd_ch(ich)%ice)
      sim_ice_thick = sd_ch(ich)%ice / max(ice_area, 1.e-6)
      ice_maturity = sim_ice_thick / max(prm%ice_maturity_ref_thick, 1.e-6)
      ice_ratio = max(0., min(1., ice_maturity))
      ice_depth_ratio = sim_ice_thick / max(sd_ch(ich)%chd, 1.0e-6)
      ice_depth_ratio = max(0., min(1., ice_depth_ratio))
      !! Numerical state cleanup before any phase/trigger logic.  These
      !! thresholds remove tiny residual ice/mobile pools that otherwise keep
      !! ratios nonzero through the warm season.  They are volume-scale relative
      !! thresholds, not calibration controls.
      if (sd_ch(ich)%ice <= ice_state_eps .and. (t_ice_decay > prm%ice_melt_tmp .or. &
          sd_ch(ich)%ice_thaw_dd >= prm%storage_cleanup_thaw_dd)) then
          ht1%flo = ht1%flo + sd_ch(ich)%ice
          sd_ch(ich)%ice = 0.
          sim_ice_thick = 0.
          ice_maturity = 0.
          ice_ratio = 0.
          ice_depth_ratio = 0.
      endif

      if (sd_ch(ich)%ice_mobile <= mobile_state_eps) sd_ch(ich)%ice_mobile = 0.
      if (sd_ch(ich)%ice_mobile_pass <= mobile_state_eps) sd_ch(ich)%ice_mobile_pass = 0.
      mobile_jam_ratio = sd_ch(ich)%ice_mobile / ice_cover_max
      mobile_pass_jam_ratio = sd_ch(ich)%ice_mobile_pass / max(ice_cover_max, 1.e-6)
      if (mobile_jam_ratio < 1.e-4) then
          sd_ch(ich)%ice_mobile = 0.
          mobile_jam_ratio = 0.
      endif
      mobile_jam_ratio = max(0., min(1., mobile_jam_ratio))
      if (mobile_pass_jam_ratio < 1.e-4) mobile_pass_jam_ratio = 0.
      mobile_pass_jam_ratio = max(0., min(1., mobile_pass_jam_ratio))

      ice_absent = (sim_ice_thick <= prm%warm_ice_thick .and. sd_ch(ich)%ice <= ice_state_eps)
      mobile_absent = (sd_ch(ich)%ice_mobile <= mobile_state_eps .and. sd_ch(ich)%ice_mobile_pass <= mobile_state_eps)

      !! ------------------------------------------------------------------
      !! Channel-scale rain-on-snow/ice diagnosis.
      !! sd_ch(ich)%ros and sd_ch(ich)%snow_melt are supplied by the HRU snowmelt
      !! routines and aggregated to each channel before routing.  Therefore this
      !! module no longer re-partitions precipitation or recomputes ROS locally.
      !! The only additional filter here is a minimum channel snowmelt depth, so
      !! trace ROS flags do not trigger breakup by themselves.
      !! ------------------------------------------------------------------
      snow_melt_mm = max(0., sd_ch(ich)%snow_melt)
      ros_day = sd_ch(ich)%ros .and. snow_melt_mm >= prm%ros_min_melt_mm

      !! ------------------------------------------------------------------
      !! Thaw-memory index for breakup and force-flush logic.
      !! Based on positive degree-days above a Tmax threshold.
      !! Rain-on-snow/ice lowers the threshold, consistent with field guidance.
      !! ------------------------------------------------------------------
      if (ros_day) then
          t_thaw = max(0., tmax - prm%thaw_tmax_base_ros)
      else
          t_thaw = max(0., tmax - prm%thaw_tmax_base)
      endif

      !! If the daily mean remains below freezing, account for nighttime refreezing
      !! by reducing the effective thaw forcing.
      if (t_air < prm%thaw_tave_base) then
          t_thaw = 0.5 * t_thaw
      endif

      sd_ch(ich)%ice_thaw_dd = prm%thaw_memory * sd_ch(ich)%ice_thaw_dd + t_thaw
      if (t_freeze > 0.) sd_ch(ich)%ice_thaw_dd = 0.7 * sd_ch(ich)%ice_thaw_dd
      if (sd_ch(ich)%ice_thaw_dd < 1.e-6) sd_ch(ich)%ice_thaw_dd = 0.
      thaw_weakening_index = sd_ch(ich)%ice_thaw_dd / &
              max(sd_ch(ich)%ice_freeze_dd + sd_ch(ich)%ice_thaw_dd, 1.e-6)
      thaw_weakening_index = max(0., min(1., thaw_weakening_index))
      drift_weak_zone = thaw_weakening_index > prm%freezeup_strong_index .and. &
                        thaw_weakening_index < prm%jam_release_weakening_index

      !! ------------------------------------------------------------------
      !! Trigger variables.
      !! qratio/qrise and minor/major event logic no longer control BREAKUP.
      !! These diagnostics are retained only where needed for non-BREAKUP
      !! disturbance/mobilization bookkeeping.
      !! ------------------------------------------------------------------
      thermal_trigger = (t_thaw > 0. .and. &
              sd_ch(ich)%ice_thaw_dd >= prm%storage_cleanup_thaw_dd)

      !! ------------------------------------------------------------------
      !! Phase state transition.
      !! Seasonal phase transitions use absolute reach-mean ice thickness and
      !! follow a one-way sequence:
      !!   WARM -> FREEZEUP -> DEEPWINTER -> BREAKUP -> WARM.
      !! Normal seasonal BREAKUP is allowed only from DEEPWINTER.  In FREEZEUP,
      !! thaw / rising-flow disturbances may trigger jam disturbance events, but they
      !! do not change the seasonal phase to BREAKUP.
      !! ------------------------------------------------------------------
      old_phase = sd_ch(ich)%ice_phase
      new_phase = old_phase

      select case (old_phase)

      case (ICE_WARM)
          !! Single-sequence seasonal logic.  WARM -> FREEZEUP is not simply
          !! "a little ice exists today"; after a spring BREAKUP, a warm-season
          !! lockout prevents a brief April cold snap from starting a new ice season.
          !! A cold-start simulation with phase_days <= 1 is still allowed to enter
          !! FREEZEUP if the freeze signal is already strong.
          warm_to_freezeup_ready = &
              (sd_ch(ich)%ice_phase_days >= prm%warm_min_days_before_freezeup) .or. &
              (sd_ch(ich)%ice_phase_days <= 1 .and. &
               sd_ch(ich)%ice_freeze_dd >= 2.0 * prm%freezeup_freeze_dd)

          if (warm_to_freezeup_ready .and. &
              sd_ch(ich)%ice_freeze_dd >= prm%freezeup_freeze_dd .and. &
              sim_ice_thick >= prm%freezeup_ice_thick .and. &
              thaw_weakening_index <= prm%freezeup_strong_index) then
              new_phase = ICE_FREEZEUP
          endif

      case (ICE_FREEZEUP)
          !! FREEZEUP matures into DEEPWINTER once a stable cover exists.
          !! The max-days guard preserves the one-way seasonal sequence while
          !! preventing a weak/variable FREEZEUP from persisting for months.
          freezeup_max_ready = sd_ch(ich)%ice_phase_days >= prm%freezeup_max_days
          if ((sd_ch(ich)%ice_phase_days >= prm%freezeup_min_days .and. &
              sim_ice_thick >= prm%deepwinter_ice_thick .and. &
               thaw_weakening_index <= prm%freezeup_strong_index) .or. &
              (freezeup_max_ready .and. &
               sim_ice_thick >= prm%freezeup_ice_thick .and. &
               thaw_weakening_index <= prm%breakup_onset_weakening_index)) then
              new_phase = ICE_DEEPWINTER
          endif

      case (ICE_DEEPWINTER)
          !! DEEPWINTER -> BREAKUP is controlled by phase duration and thermal
          !! weakening only.  ROS/snowmelt no longer bypasses the minimum
          !! DEEPWINTER duration; it can affect thaw_weak and ice melt, but not
          !! directly trigger phase change.  The max-days guard prevents the
          !! model from being trapped in DEEPWINTER after the seasonal thaw window.
          breakup_material_ready = &
              sim_ice_thick >= prm%warm_ice_thick .or. &
              sd_ch(ich)%ice > ice_state_eps .or. &
              sd_ch(ich)%ice_cover_stor > jam_stor_eps .or. &
              sd_ch(ich)%ice_mobile > mobile_state_eps .or. &
              sd_ch(ich)%ice_mobile_pass > mobile_state_eps

              breakup_thermal_ready = thaw_weakening_index >= prm%breakup_onset_weakening_index
          deepwinter_age_ready = &
              sd_ch(ich)%ice_phase_days >= prm%deepwinter_min_days_before_breakup
          deepwinter_max_ready = &
              sd_ch(ich)%ice_phase_days >= prm%deepwinter_max_days_before_breakup .and. &
              thaw_weakening_index >= prm%freezeup_strong_index

          if (breakup_material_ready .and. &
              ((deepwinter_age_ready .and. breakup_thermal_ready) .or. deepwinter_max_ready)) then
                  new_phase = ICE_BREAKUP
              endif

      case (ICE_BREAKUP)
          !! BREAKUP -> WARM is a confirmed seasonal exit, not a temporary warm
          !! day or momentary loss of local ice thickness.  A max-days guard
          !! allows cleanup once ice is small and the breakup phase has persisted
          !! for a long period, even if tiny storage residuals remain.
          breakup_long_enough = &
              sd_ch(ich)%ice_phase_days >= prm%breakup_min_days_before_warm
          breakup_max_ready = &
              sd_ch(ich)%ice_phase_days >= prm%breakup_max_days_before_warm
          thermal_warm_ready = &
              thaw_weakening_index >= prm%warm_season_weakening_index .and. &
              sd_ch(ich)%ice_thaw_dd >= prm%flush_thaw_dd
          ice_small_enough = sim_ice_thick <= prm%warm_ice_thick
          storage_not_active = &
              sd_ch(ich)%ice_jam_stor <= prm%warm_storage_exit_ratio * max(jam_stor_max, 1.e-6) .and. &
              sd_ch(ich)%ice_cover_stor <= prm%warm_storage_exit_ratio * max(cover_stor_max, 1.e-6)

          if ((breakup_long_enough .and. thermal_warm_ready .and. &
               ice_small_enough .and. storage_not_active) .or. &
              (breakup_max_ready .and. ice_small_enough .and. &
               thaw_weakening_index >= prm%jam_release_weakening_index)) then
              new_phase = ICE_WARM
          endif

      case default
          new_phase = ICE_WARM

      end select

      phase_changed_today = (new_phase /= old_phase)
      breakup_onset_today = (old_phase == ICE_DEEPWINTER .and. new_phase == ICE_BREAKUP)
      seasonal_breakup_reset = (old_phase == ICE_BREAKUP .and. new_phase == ICE_WARM)

      if (new_phase == old_phase) then
          sd_ch(ich)%ice_phase_days = sd_ch(ich)%ice_phase_days + 1
      else
          sd_ch(ich)%ice_phase = new_phase
          sd_ch(ich)%ice_phase_days = 1
      endif
      !! Seasonal BREAKUP onset is a transition day: high cover_to_jam and
      !! ice-to-mobile conversion are allowed, but mature jam block/release starts
      !! only from following BREAKUP days.

      !! When BREAKUP has ended, the previous ice season is over.  Return
      !! residual ice mass and ice-related water storage to the routing water
      !! balance, then reset all seasonal ice-jam states so the next freeze-up
      !! starts from zero.  This preserves mass balance because these pools were
      !! originally formed by removing water from ht1%flo/ch_stor.
      if (seasonal_breakup_reset) then
          warm_cleanup_return = sd_ch(ich)%ice + sd_ch(ich)%ice_mobile + &
                                sd_ch(ich)%ice_mobile_pass + &
                                sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor
          warm_cleanup_return = max(0., warm_cleanup_return)
          ht1%flo = ht1%flo + warm_cleanup_return

          sd_ch(ich)%ice = 0.
          sd_ch(ich)%ice_mobile = 0.
          sd_ch(ich)%ice_mobile_pass = 0.
          sd_ch(ich)%ice_cover_stor = 0.
          sd_ch(ich)%ice_jam_stor = 0.
          sd_ch(ich)%ice_release_active = JAM_NONE
          sd_ch(ich)%ice_release_days = 0
          sd_ch(ich)%ice_block_days = 0
          total_ice_stor = 0.
          sim_ice_thick = 0.
          ice_maturity = 0.
          ice_ratio = 0.
          mobile_jam_ratio = 0.
      endif

      !! WARM handling / cleanup.
      !! Only the explicit BREAKUP->WARM seasonal reset clears the previous ice
      !! season.  There is no daily WARM flushing here, so late-WARM freeze-up
      !! ice can accumulate before the next FREEZEUP.

      !! ------------------------------------------------------------------
      !! Stable ice-cover retention.
      !! This process is active only during FREEZEUP and DEEPWINTER.  BREAKUP
      !! is handled by explicit jam block/release; allowing cover retention in
      !! BREAKUP caused cross-season storage accumulation and phase lock.
      !! ------------------------------------------------------------------
      if (sd_ch(ich)%ice_phase == ICE_FREEZEUP .or. &
          sd_ch(ich)%ice_phase == ICE_DEEPWINTER) then
          if (sim_ice_thick >= prm%retention_ice_thick) then
              q_underice_cap = q_jam_ref_rate * max(prm%underice_cap_min, &
                      prm%underice_cap_open * (1.0 - ice_ratio)**prm%underice_cap_exp)

              underice_excess = max(0., ht1%flo - q_underice_cap * 86400.)

              block_capacity = prm%ice_cover_ret_cap_coeff * sd_ch(ich)%ice * reach%jam_block_modifier
              block_capacity = max(0., block_capacity)

              !! Stable-cover retention has its own conceptual storage capacity.
              !! Do not reuse jam_stor_max/remaining_jam_capacity here; those
              !! are reserved for the BREAKUP ice-jam reservoir.
              cover_stor_max = max(0., prm%ice_cover_ret_stor_frac * &
                      reach%jam_storage_modifier * ch_vol_cap)
              cover_remaining_capacity = max(0., cover_stor_max - sd_ch(ich)%ice_cover_stor)

              !! Stable ice-cover retention is separated from breakup-jam
              !! blocking for debugging.  The effective maximum retention
              !! fraction is damped as q_ratio increases, so stable cover mainly
              !! suppresses small winter pulses rather than unrealistically
              !! blocking most of a strong hydraulic event.
              select case (sd_ch(ich)%ice_phase)
              case (ICE_FREEZEUP)
                  phase_ret_mult = prm%freezeup_ret_mult
              case (ICE_DEEPWINTER)
                  phase_ret_mult = prm%deepwinter_ret_mult
              case (ICE_BREAKUP)
                  phase_ret_mult = prm%breakup_ret_mult
              case default
                  phase_ret_mult = 0.0
              end select

              retention_frac_eff = phase_ret_mult * prm%ice_cover_ret_frac_max * ice_ratio * &
                      reach%jam_form_modifier / (1.0 + prm%ice_cover_ret_q_damp * max(0., q_ratio))
              if (sd_ch(ich)%ice_phase == ICE_DEEPWINTER) then
                  retention_frac_eff = max(retention_frac_eff, prm%deepwinter_ret_frac_min)
                  retention_frac_eff = min(retention_frac_eff, prm%deepwinter_ret_frac_max)
              endif
              retention_frac_eff = max(0., min(prm%ice_cover_ret_frac_max, retention_frac_eff))

              blocked_cover = min(block_capacity, underice_excess)
              blocked_cover = min(blocked_cover, retention_frac_eff * ht1%flo)
              blocked_cover = min(blocked_cover, cover_remaining_capacity)
              blocked_cover = max(0., blocked_cover)

              ht1%flo = ht1%flo - blocked_cover
              sd_ch(ich)%ice_cover_stor = sd_ch(ich)%ice_cover_stor + blocked_cover
              sd_ch(ich)%icejam_block = sd_ch(ich)%icejam_block + blocked_cover
              if (blocked_cover > 0.) then
                  cover_stor_ratio = sd_ch(ich)%ice_cover_stor / cover_stor_max
                  cover_stor_ratio = max(0., min(1., cover_stor_ratio))
                  total_ice_stor = sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor
              endif
              cover_q_cap = q_underice_cap
              cover_underice_excess = underice_excess
              cover_block_capacity = block_capacity
              cover_stor_max_dbg = cover_stor_max
              ! cover_remaining_capacity already set above
        endif
      endif

      !! ------------------------------------------------------------------
      !! Mechanical trigger based on ice load and support capacity.
      !! ------------------------------------------------------------------
      select case (sd_ch(ich)%ice_phase)
      case (ICE_FREEZEUP)
          ice_strength_factor = 0.7
      case (ICE_DEEPWINTER)
          ice_strength_factor = 1.5
      case (ICE_BREAKUP)
          ice_strength_factor = 0.8
      case default
          ice_strength_factor = 1.0
      end select

      ice_load = sd_ch(ich)%ice + sd_ch(ich)%ice_mobile
      ice_support_capacity = prm%ice_support_frac * ice_cap_vol * ice_strength_factor / max(reach%mechanical_weakness_modifier, 0.2)

      if (ice_load > ice_support_capacity .and. q_in_rate_raw >= 0.0) then
          mechanical_trigger = .true.
      endif

      !! ------------------------------------------------------------------
      !! Pre-event mobile ice generation.
      !!
      !! ice_mobile is not only produced by explicit jam-breakup release.
      !! It also represents broken/floating ice that can be generated by
      !! progressive thaw, rising flow, or mechanical weakening before a
      !! mature jam release occurs.  This allows upstream broken ice to be
      !! advected downstream and captured in susceptible reaches, which is
      !! essential for representing large breakup-jam development.
      !! ------------------------------------------------------------------
      ice_mobilized_drift = 0.
      ice_mobilized_dynamic = 0.
      ice_mobilized_pre = 0.
      drift_frac_eff = 0.
      dynamic_frac_eff = 0.

      if (sd_ch(ich)%ice > ice_stor_eps .and. &
              q_in_rate_raw >= prm%mobile_q_min .and. sim_ice_thick >= prm%mobile_ice_thick) then

          !! Background drift: weak, progressive conversion of local cover ice
          !! into transportable broken ice during FREEZEUP/BREAKUP progression.
          select case (sd_ch(ich)%ice_phase)
          case (ICE_BREAKUP)
              if (thermal_trigger .or. ros_day .or. &
                  sd_ch(ich)%ice_thaw_dd >= prm%mobile_thaw_dd) then
                  drift_frac_eff = prm%drift_mobilization_frac * prm%mobile_breakup_drift_multiplier
              endif
          case (ICE_FREEZEUP)
              if (sd_ch(ich)%ice_thaw_dd >= prm%mobile_thaw_dd .and. &
                  q_ratio >= 0.0) then
                  drift_frac_eff = prm%drift_mobilization_frac * prm%mobile_freezeup_drift_multiplier
              endif
          case (ICE_DEEPWINTER)
              !! Deep-winter mobile ice generation is strongly damped and only
              !! responds to meaningful thaw/ROS disturbance.
              if (thermal_trigger .and. ros_day) then
                  drift_frac_eff = prm%drift_mobilization_frac * prm%mobile_deepwinter_dynamic_weight
              endif
          end select

          drift_frac_eff = drift_frac_eff * mobile_order_mult
          ice_mobilized_drift = drift_frac_eff * sd_ch(ich)%ice

          !! Dynamic/mechanical mobilization: stronger conversion caused by
          !! ROS forcing or load-induced ice-cover failure.  qratio/qrise no
          !! longer control this pathway.
          if (mechanical_trigger .and. &
              (thermal_trigger .or. ros_day .or. sd_ch(ich)%ice_phase == ICE_BREAKUP)) then
              dynamic_frac_eff = prm%dynamic_mobilization_frac
          endif

          if (sd_ch(ich)%ice_phase == ICE_DEEPWINTER) then
              dynamic_frac_eff = dynamic_frac_eff * prm%mobile_deepwinter_dynamic_weight
          endif

          dynamic_frac_eff = dynamic_frac_eff * mobile_order_mult
          ice_mobilized_dynamic = dynamic_frac_eff * sd_ch(ich)%ice

          ice_mobilized_pre = ice_mobilized_drift + ice_mobilized_dynamic
          ice_mobilized_pre = min(ice_mobilized_pre, prm%mobile_max_daily_frac * sd_ch(ich)%ice)
          ice_mobilized_pre = max(0., min(ice_mobilized_pre, sd_ch(ich)%ice))

          if (ice_mobilized_pre > 0.) then
              sd_ch(ich)%ice = sd_ch(ich)%ice - ice_mobilized_pre
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile + ice_mobilized_pre
              sd_ch(ich)%ice = max(0., sd_ch(ich)%ice)
              sim_ice_thick = sd_ch(ich)%ice / max(ice_area, 1.e-6)
              ice_maturity = sim_ice_thick / max(prm%ice_maturity_ref_thick, 1.e-6)
              ice_ratio = max(0., min(1., ice_maturity))
          endif
          mobile_jam_ratio = sd_ch(ich)%ice_mobile / ice_cover_max
          mobile_pass_jam_ratio = sd_ch(ich)%ice_mobile_pass / max(ice_cover_max, 1.e-6)
          mobile_jam_ratio = max(0., min(1., mobile_jam_ratio))
          mobile_pass_jam_ratio = max(0., min(1., mobile_pass_jam_ratio))
      endif

      !! ------------------------------------------------------------------
      !! Release existing ice-related stored water.
      !!
      !! Two storage pools are kept separate:
      !!   ice_cover_stor : water retained by stable under-ice/ice-cover conveyance reduction
      !!   ice_jam_stor   : water blocked by explicit freeze-up/breakup ice-jam formation
      !!
      !! Release thresholds therefore use source-specific maturity ratios:
      !!   cover_stor_ratio      = ice_cover_stor / stable-cover storage capacity
      !!   jam_stor_ratio  = ice_jam_stor / unified jam storage capacity
      !! ------------------------------------------------------------------
      if (sd_ch(ich)%ice_cover_stor <= jam_stor_eps) sd_ch(ich)%ice_cover_stor = 0.
      if (sd_ch(ich)%ice_jam_stor <= jam_stor_eps) sd_ch(ich)%ice_jam_stor = 0.

      cover_stor_ratio = sd_ch(ich)%ice_cover_stor / cover_stor_max
      cover_stor_ratio = max(0., min(1., cover_stor_ratio))
      jam_stor_ratio = sd_ch(ich)%ice_jam_stor / max(jam_stor_max, 1.e-6)
      jam_stor_ratio = max(0., min(1., jam_stor_ratio))
      total_ice_stor = sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor

      !! v17: save the start-of-day jam storage before any release or blocking.
      !! New jam additions later in the day (cover_to_jam + block_jam) are
      !! constrained by this start-of-day free capacity, so storage freed by
      !! same-day release cannot be immediately refilled.  This adds a daily
      !! routing lag and prevents artificial release/block cancellation.
      jam_stor_start = sd_ch(ich)%ice_jam_stor
      jam_remain_capacity_start = max(0., jam_stor_max - jam_stor_start)
      jam_capacity_used_today = 0.

      !! BREAKUP internal state machine.  The seasonal onset day is a DRIFT
      !! transition day: it can move cover_stor into jam_stor and mobilize ice,
      !! but does not perform mature jam blocking or release.  Subsequent days
      !! are mutually exclusive DRIFT / BLOCK / RELEASE states.
      breakup_release_gate = .false.
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
          breakup_release_gate = (thaw_weakening_index >= prm%jam_release_weakening_index)

          if (sd_ch(ich)%ice_release_active /= JAM_NONE) then
          if (sd_ch(ich)%ice_jam_stor <= jam_stor_eps) then
              sd_ch(ich)%ice_release_active = JAM_NONE
              sd_ch(ich)%ice_release_days = 0
              else
                  sd_ch(ich)%ice_release_days = sd_ch(ich)%ice_release_days + 1
                  sd_ch(ich)%ice_block_days = 0
                  if (sd_ch(ich)%ice_release_days > prm%release_auto_hold_days .and. &
                      thaw_weakening_index < 0.8 * prm%jam_release_weakening_index) then
                      sd_ch(ich)%ice_release_active = JAM_NONE
                      sd_ch(ich)%ice_release_days = 0
                  endif
          endif
      endif

          if (sd_ch(ich)%ice_release_active /= JAM_NONE) then
              breakup_day_type = BRK_DAY_RELEASE
          else if (breakup_onset_today) then
              breakup_day_type = BRK_DAY_DRIFT
              sd_ch(ich)%ice_block_days = 0
              sd_ch(ich)%ice_release_days = 0
          else if (breakup_release_gate .and. sd_ch(ich)%ice_jam_stor > jam_stor_eps) then
              sd_ch(ich)%ice_release_active = JAM_RELEASE
              sd_ch(ich)%ice_release_days = 1
              sd_ch(ich)%ice_block_days = 0
              breakup_day_type = BRK_DAY_RELEASE
          else if (thaw_weakening_index < prm%jam_release_weakening_index .and. &
               block_flow_ready .and. &
               jam_remain_capacity_start > jam_stor_eps .and. &
               (mobile_jam_ratio >= prm%jam_mobile_trigger_ratio .or. &
                mobile_pass_jam_ratio >= prm%jam_mobile_trigger_ratio .or. &
                (sd_ch(ich)%ice_cover_stor > jam_stor_eps .and. &
                 (mobile_jam_ratio + mobile_pass_jam_ratio) >= 0.5 * prm%jam_mobile_trigger_ratio))) then
              breakup_day_type = BRK_DAY_BLOCK
              sd_ch(ich)%ice_block_days = sd_ch(ich)%ice_block_days + 1
              sd_ch(ich)%ice_release_days = 0
          else
              breakup_day_type = BRK_DAY_DRIFT
              sd_ch(ich)%ice_block_days = 0
              sd_ch(ich)%ice_release_days = 0
          endif
      else
          sd_ch(ich)%ice_release_active = JAM_NONE
          sd_ch(ich)%ice_release_days = 0
          sd_ch(ich)%ice_block_days = 0
      endif

      if (total_ice_stor <= jam_stor_eps) then
          sd_ch(ich)%ice_cover_stor = 0.
          sd_ch(ich)%ice_jam_stor = 0.
          sd_ch(ich)%ice_release_active = JAM_NONE
          sd_ch(ich)%ice_release_days = 0
      endif
      
      if (total_ice_stor > jam_stor_eps) then
          released_event = 0.
          released_jam_event = 0.
          released_leak = 0.
          released_cover_leak = 0.
          released_jam_leak = 0.
          jam_release_ratio = 0.
          mobilization_ratio = 0.
          cover_stor_before = sd_ch(ich)%ice_cover_stor
          jam_stor_before = sd_ch(ich)%ice_jam_stor
          stor_before = total_ice_stor

          !! BREAKUP release is controlled only by thaw weakening and episode
          !! memory.  All older minor/major, qratio/qrise, force-flush, and
          !! coupled-cover release paths have been removed.
          if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
              if (sd_ch(ich)%ice_release_active /= JAM_NONE .and. &
                  sd_ch(ich)%ice_jam_stor > jam_stor_eps) then

                  release_recession_day = .true.
                  breakup_day_type = BRK_DAY_RELEASE

                  release_weak_eff = max(0., min(1., thaw_weakening_index)) ** prm%jam_release_weak_exp
                  release_frac_eff = prm%jam_leak_frac + &
                      release_weak_eff * (prm%jam_release_frac_max - prm%jam_leak_frac)
                  release_frac_eff = max(prm%jam_leak_frac, &
                      min(prm%jam_release_frac_max, release_frac_eff))

                  if (sd_ch(ich)%ice_release_days <= 1) then
                      release_ramp_factor = prm%release_ramp_day1_frac
                  else if (sd_ch(ich)%ice_release_days == 2) then
                      release_ramp_factor = prm%release_ramp_day2_frac
                  else
                      release_ramp_factor = 1.0
                  endif
                  release_ramp_factor = max(0., min(1., release_ramp_factor))
                  !! Strong HRU-diagnosed ROS/snowmelt events can produce rapid
                  !! mechanical breakup.  In those cases, avoid strongly damping
                  !! the first release day with the generic episode ramp.
                  if (ros_day .and. snow_melt_mm >= prm%ros_min_melt_mm) then
                      release_ramp_factor = max(release_ramp_factor, 0.80)
                  endif
                  release_frac_eff = release_frac_eff * release_ramp_factor
                  recession_frac = release_frac_eff

                  released_jam_event = release_frac_eff * sd_ch(ich)%ice_jam_stor
                  released_jam_event = max(0., min(released_jam_event, sd_ch(ich)%ice_jam_stor))
                  if (released_jam_event > active_release_min) sd_ch(ich)%ice_jam_flag = JAM_RELEASE
                  endif
              endif

          !! No direct or coupled cover-storage event release is allowed.  Cover
          !! storage can become jam storage, leak slowly, or be cleared only by the
          !! explicit seasonal BREAKUP->WARM reset.
    
          released_jam_event = max(0., min(released_jam_event, sd_ch(ich)%ice_jam_stor))
          if (jam_stor_before > jam_stor_eps) then
              jam_release_ratio = released_jam_event / jam_stor_before
          else
              jam_release_ratio = 0.
          endif
          jam_release_ratio = max(0., min(1., jam_release_ratio))

          released_event = released_jam_event
          if (stor_before > jam_stor_eps) then
              release_ratio = released_event / stor_before
          else
              release_ratio = 0.
          endif
              release_ratio = max(0., min(1., release_ratio))

          !! Jam release can mobilize local cover ice into transportable floes.
          mobilization_ratio = jam_release_ratio
          mobilization_ratio = max(0., min(1., mobilization_ratio))

          ice_mobilized_event = 0.
          if (released_event > active_release_min) then
              ice_mobilized_event = mobilization_ratio * sd_ch(ich)%ice
              ice_mobilized_event = max(0., min(ice_mobilized_event, sd_ch(ich)%ice))
          endif

          sd_ch(ich)%ice = sd_ch(ich)%ice - ice_mobilized_event
          sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile + ice_mobilized_event
          sd_ch(ich)%ice = max(0., sd_ch(ich)%ice)
          sim_ice_thick = sd_ch(ich)%ice / max(ice_area, 1.e-6)
          ice_maturity = sim_ice_thick / max(prm%ice_maturity_ref_thick, 1.e-6)
          ice_ratio = max(0., min(1., ice_maturity))

          sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor - released_jam_event

          cover_stor_ratio = sd_ch(ich)%ice_cover_stor / max(cover_stor_max, 1.e-6)
          cover_stor_ratio = max(0., min(1., cover_stor_ratio))
          jam_stor_ratio = sd_ch(ich)%ice_jam_stor / max(jam_stor_max, 1.e-6)
          jam_stor_ratio = max(0., min(1., jam_stor_ratio))

          !! Background leakage is not an event and does not alter release mode.
          if (released_event <= active_release_min) then
              if (sd_ch(ich)%ice_phase == ICE_DEEPWINTER .and. .not. ros_day) then
                  released_cover_leak = max(released_cover_leak, &
                          prm%deepwinter_cover_leak_frac * sd_ch(ich)%ice_cover_stor)
              endif

              if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
                  released_jam_leak = max(released_jam_leak, &
                          prm%jam_leak_frac * sd_ch(ich)%ice_jam_stor)

                  if (breakup_day_type == BRK_DAY_DRIFT) then
                      !! DRIFT is not a strong release episode, but it should
                      !! represent open-flow / pass-through recession rather than
                      !! near-zero leakage.  Use a derived intermediate fraction
                      !! between background leak and maximum release without adding
                      !! a new calibration parameter.
                      released_jam_leak = max(released_jam_leak, &
                              sqrt(prm%jam_leak_frac * prm%jam_release_frac_max) * &
                              sd_ch(ich)%ice_jam_stor)
                  released_cover_leak = max(released_cover_leak, &
                              sqrt(prm%jam_leak_frac * prm%jam_release_frac_max) * &
                              sd_ch(ich)%ice_cover_stor)
                  else if (breakup_release_gate .and. breakup_day_type == BRK_DAY_RELEASE) then
                      released_cover_leak = max(released_cover_leak, &
                              prm%jam_leak_frac * sd_ch(ich)%ice_cover_stor)
                  endif
              endif
          endif

          released_cover_leak = max(0., min(released_cover_leak, sd_ch(ich)%ice_cover_stor))
          released_jam_leak = max(0., min(released_jam_leak, sd_ch(ich)%ice_jam_stor))
          released_leak = released_cover_leak + released_jam_leak

          sd_ch(ich)%ice_cover_stor = sd_ch(ich)%ice_cover_stor - released_cover_leak
          sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor - released_jam_leak

          !! Storage bounds check.  Blocking days are exempt so jam storage can
          !! accumulate.  WARM does not trigger daily flushing here; the seasonal
          !! BREAKUP->WARM reset has already cleared the previous ice season.
          cover_stor_cap = cover_stor_max
          jam_stor_cap = jam_stor_max
          cover_overflow_release = 0.
          jam_overflow_release = 0.

          if (sd_ch(ich)%ice_phase /= ICE_WARM .and. &
              breakup_day_type == BRK_DAY_RELEASE) then
                  cover_overflow_release = max(0., sd_ch(ich)%ice_cover_stor - cover_stor_cap)
                  jam_overflow_release = max(0., sd_ch(ich)%ice_jam_stor - jam_stor_cap)
              endif

          cover_overflow_release = max(0., min(cover_overflow_release, sd_ch(ich)%ice_cover_stor))
          jam_overflow_release = max(0., min(jam_overflow_release, sd_ch(ich)%ice_jam_stor))

          if (cover_overflow_release > 0.) then
              sd_ch(ich)%ice_cover_stor = sd_ch(ich)%ice_cover_stor - cover_overflow_release
              released_cover_leak = released_cover_leak + cover_overflow_release
          endif
          if (jam_overflow_release > 0.) then
              sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor - jam_overflow_release
              released_jam_leak = released_jam_leak + jam_overflow_release
          endif
          released_leak = released_cover_leak + released_jam_leak

          if (sd_ch(ich)%ice_cover_stor < 1.e-6) sd_ch(ich)%ice_cover_stor = 0.
          if (sd_ch(ich)%ice_jam_stor < 1.e-6) sd_ch(ich)%ice_jam_stor = 0.

          cover_stor_ratio = sd_ch(ich)%ice_cover_stor / max(cover_stor_max, 1.e-6)
          cover_stor_ratio = max(0., min(1., cover_stor_ratio))
          jam_stor_ratio = sd_ch(ich)%ice_jam_stor / max(jam_stor_max, 1.e-6)
          jam_stor_ratio = max(0., min(1., jam_stor_ratio))
          total_ice_stor = sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor

          if (ice_absent .and. mobile_absent .and. total_ice_stor <= warm_storage_exit_threshold) then
              sd_ch(ich)%ice_cover_stor = 0.
              sd_ch(ich)%ice_jam_stor = 0.
              total_ice_stor = 0.
              sd_ch(ich)%ice_release_active = JAM_NONE
              sd_ch(ich)%ice_release_days = 0
          endif

          released = released_event + released_leak
          ht1%flo = ht1%flo + released
          sd_ch(ich)%icejam_release = released

          if (total_ice_stor <= jam_stor_eps) then
              sd_ch(ich)%ice_cover_stor = 0.
              sd_ch(ich)%ice_jam_stor = 0.
              sd_ch(ich)%ice_release_active = JAM_NONE
          endif

          endif

      ice_mobilized_total = ice_mobilized_pre + ice_mobilized_event + ice_mobilized_cover_break
      ice_mobilized = ice_mobilized_total

      !! Split newly generated local mobile ice into a locally captured jam-material
      !! pool (ice_mobile) and a pass-through pool (ice_mobile_pass) that will be
      !! routed on the next daily step.  Incoming pass-through ice is still handled
      !! exclusively by sd_channel_ice_advect.  This prevents all newly generated
      !! local mobile ice from immediately leaving the reach, while still allowing
      !! a pass-through component to move downstream with a one-day delay.
      if (ice_mobilized_total > mobile_state_eps) then
          local_mobile_capture_frac = prm%mobile_capture_base + &
                  prm%mobile_capture_susc_weight * reach%ice_capture_modifier + &
                  prm%mobile_capture_ice_weight * ice_ratio + &
                  prm%mobile_capture_depth_weight * ice_depth_ratio
          select case (sd_ch(ich)%ice_phase)
          case (ICE_FREEZEUP)
              local_mobile_capture_frac = max(local_mobile_capture_frac, prm%freezeup_capture_min)
          case (ICE_DEEPWINTER)
              local_mobile_capture_frac = max(local_mobile_capture_frac, prm%deepwinter_capture_min)
          case (ICE_BREAKUP)
              local_mobile_capture_frac = max(local_mobile_capture_frac, prm%breakup_capture_min)
          case default
              if (ice_ratio <= 0.05) local_mobile_capture_frac = min(local_mobile_capture_frac, prm%warm_capture_max)
          end select
          local_mobile_capture_frac = max(prm%mobile_capture_min, min(prm%mobile_capture_max, local_mobile_capture_frac))
          local_mobile_pass_frac = max(0., 1.0 - local_mobile_capture_frac)
          ice_mobile_generated_pass = min(sd_ch(ich)%ice_mobile, local_mobile_pass_frac * ice_mobilized_total)
          if (ice_mobile_generated_pass > 0.) then
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile - ice_mobile_generated_pass
              sd_ch(ich)%ice_mobile_pass = sd_ch(ich)%ice_mobile_pass + ice_mobile_generated_pass
          endif
      endif

      !! Existing mobile ice should not remain as a static local stock until
      !! WARM cleanup.  Jam release flushes a fraction of the local mobile pool
      !! into the pass-through pool, using jam_release_ratio as the release
      !! strength proxy.  DRIFT days also pass a small fraction downstream.
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP .and. sd_ch(ich)%ice_mobile > mobile_state_eps) then
          if (breakup_day_type == BRK_DAY_RELEASE .and. released_jam_event > active_release_min) then
              mobile_flush_frac = max(0., min(1., jam_release_ratio))
              mobile_flushed_by_release = mobile_flush_frac * sd_ch(ich)%ice_mobile
              mobile_flushed_by_release = max(0., min(mobile_flushed_by_release, sd_ch(ich)%ice_mobile))
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile - mobile_flushed_by_release
              sd_ch(ich)%ice_mobile_pass = sd_ch(ich)%ice_mobile_pass + mobile_flushed_by_release
          else if (breakup_day_type == BRK_DAY_DRIFT) then
              mobile_drift_pass_frac = prm%drift_mobilization_frac * prm%mobile_breakup_drift_multiplier
              mobile_drift_pass_frac = max(0., min(1., mobile_drift_pass_frac))
              mobile_drift_pass = mobile_drift_pass_frac * sd_ch(ich)%ice_mobile
              mobile_drift_pass = max(0., min(mobile_drift_pass, sd_ch(ich)%ice_mobile))
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile - mobile_drift_pass
              sd_ch(ich)%ice_mobile_pass = sd_ch(ich)%ice_mobile_pass + mobile_drift_pass
          endif
          sd_ch(ich)%ice_mobile = max(0., sd_ch(ich)%ice_mobile)
      endif

      !! Recompute the mobile-ice ratio after event mobilization.
      !! This is essential when a cover-only release generates mobile ice:
      !! the newly mobilized ice can be captured and form a jam on the same day.
      mobile_jam_ratio = sd_ch(ich)%ice_mobile / ice_cover_max
      mobile_pass_jam_ratio = sd_ch(ich)%ice_mobile_pass / max(ice_cover_max, 1.e-6)
      mobile_jam_ratio = max(0., min(1., mobile_jam_ratio))
      mobile_pass_jam_ratio = max(0., min(1., mobile_pass_jam_ratio))

      !! ------------------------------------------------------------------
      !! 4. Determine whether same-day cover-to-jam or mature jam blocking can
      !! occur.  BREAKUP states are mutually exclusive:
      !!   DRIFT   : cover_to_jam may occur on onset day, but no block/release.
      !!   BLOCK   : cover_to_jam + block_jam can occur.
      !!   RELEASE : jam_stor release can occur; no cover_to_jam/block_jam.
      !! ------------------------------------------------------------------
      allow_new_jam_today = .false.
      do_jam_formation_today = .false.
      do_cover_to_jam_today = .false.
      do_onset_cover_block_today = .false.
      jam_formation_ready = .false.

          select case (sd_ch(ich)%ice_phase)
          case (ICE_FREEZEUP)
          allow_new_jam_today = (released_event <= active_release_min .and. &
                                 released_jam_event <= active_release_min)
          if (allow_new_jam_today) then
              jam_formation_ready = (sim_ice_thick >= prm%jam_material_ice_thick .or. &
                   mobile_jam_ratio >= prm%jam_mobile_trigger_ratio) .and. &
                   (thermal_trigger .or. ros_day .or. mechanical_trigger)
          endif

          case (ICE_BREAKUP)
          do_cover_to_jam_today = (breakup_onset_today .or. breakup_day_type == BRK_DAY_BLOCK)
          do_onset_cover_block_today = breakup_onset_today
          do_jam_formation_today = (breakup_day_type == BRK_DAY_BLOCK .or. do_onset_cover_block_today)
          if (do_jam_formation_today) then
              !! v38: Mature BLOCK requires actual ice material.  Residual
              !! ice_jam_stor is water storage, not a source of new jam
              !! material, and therefore must not by itself trigger new
              !! blocking.  This allows mobile-ice pass-through / DRIFT days
              !! when floating ice is present but insufficient to form an
              !! effective jam.
              jam_formation_ready = breakup_onset_today .or. &
                  (block_flow_ready .and. &
                   (mobile_jam_ratio >= prm%jam_mobile_trigger_ratio .or. &
                  mobile_pass_jam_ratio >= prm%jam_mobile_trigger_ratio .or. &
                  (sd_ch(ich)%ice_cover_stor > jam_stor_eps .and. &
                     (mobile_jam_ratio + mobile_pass_jam_ratio) >= 0.5 * prm%jam_mobile_trigger_ratio)))
          else
              jam_formation_ready = do_cover_to_jam_today
          endif

          case default
              jam_formation_ready = .false.
          end select

      if (sd_ch(ich)%ice_phase == ICE_FREEZEUP) then
          do_jam_formation_today = allow_new_jam_today .and. jam_formation_ready
          do_cover_to_jam_today = .false.
      endif

      !! ------------------------------------------------------------------
      !! Mobile ice remains mobile jam material.  Do not convert captured or
      !! locally present mobile ice into sd_ch%ice, because sd_ch%ice is used to
      !! diagnose stable ice-cover thickness and phase.  Jam blocking capacity
      !! below uses sd_ch%ice + sd_ch%ice_mobile as the available ice load.
      !! ------------------------------------------------------------------

      !! ------------------------------------------------------------------
      !! Set blocking parameters and form an explicit ice jam.
      !! On breakup onset, reclassify mature cover storage into jam-controlled
      !! storage before blocking new incoming water.  This order represents
      !! stable-cover backwater becoming localized jam backwater before the
      !! newly arriving flow is impounded behind the jam.
      !! ------------------------------------------------------------------
      block_capacity = 0.
      block_cap_coeff = 0.
      block_frac_max = 0.
      jam_stor_max_frac_eff = prm%jam_form_stor_max_frac * reach%jam_storage_modifier
      jam_stor_max = max(0., jam_stor_max_frac_eff * ch_vol_cap)
      jam_remaining_capacity_step = 0.
      cover_to_jam = 0.
      cover_to_jam_frac_eff = 0.

      if (do_jam_formation_today) then
          !! Mature BREAKUP block uses the unified jam-blocking rule.  On the
          !! seasonal onset day, apply only weak residual cover-controlled
          !! obstruction: the reach is transitioning from stable cover to mobile
          !! ice, but a mature ice jam is not yet assumed.
          block_cap_coeff = prm%jam_form_block_cap_coeff
          block_frac_max = prm%jam_form_block_frac_max
          if (do_onset_cover_block_today) then
              !! v38: derive weak residual cover obstruction from the onset
              !! cover-to-jam transfer fraction instead of exposing another
              !! parameter.  A larger cover_to_jam transfer leaves less
              !! residual cover-controlled obstruction on the onset day.
              onset_block_mult = max(0., min(1., 1.0 - prm%breakup_onset_cover_to_jam_frac))
              block_cap_coeff = block_cap_coeff * onset_block_mult
              block_frac_max = block_frac_max * onset_block_mult
          endif
          jam_stor_max_frac_eff = prm%jam_form_stor_max_frac * reach%jam_storage_modifier
          jam_stor_max = max(0., jam_stor_max_frac_eff * ch_vol_cap)
      endif

      !! 6a. Existing cover-controlled storage becomes jam-controlled storage first.
      !! v12 uses one cover_breakup_frac because the model should not know at
      !! jam-buildup time how large the eventual release will be.
      !! The same fraction also mobilizes local cover ice into mobile ice,
      !! representing synchronous breakup of the ice cover.
      if (sd_ch(ich)%ice_cover_stor > jam_stor_eps .and. jam_stor_max > jam_stor_eps) then
          if (do_cover_to_jam_today .and. sd_ch(ich)%ice_phase == ICE_BREAKUP) then
              !! v25: distinguish the initial DEEPWINTER->BREAKUP transition
              !! from ordinary BREAKUP rebuilding.  At seasonal breakup onset,
              !! the reach commonly has abundant stable-cover storage and ice;
              !! a large fraction can be reclassified as jam-controlled storage.
              !! Later rebuilding uses the ordinary cover_breakup_frac.
              if (breakup_onset_today) then
                  cover_to_jam_frac_eff = prm%breakup_onset_cover_to_jam_frac
              else
              cover_to_jam_frac_eff = prm%cover_breakup_frac
          endif
          endif

          if (cover_to_jam_frac_eff > 0.) then
              !! v17: cover-to-jam also uses only start-of-day free jam capacity.
              !! This keeps same-day release from creating immediately reusable
              !! capacity and preserves a daily lag in jam buildup.
              cover_to_jam_capacity = max(0., jam_remain_capacity_start - jam_capacity_used_today)
              cover_to_jam = min(cover_to_jam_frac_eff * sd_ch(ich)%ice_cover_stor, cover_to_jam_capacity)
              cover_to_jam = max(0., min(cover_to_jam, sd_ch(ich)%ice_cover_stor))

          sd_ch(ich)%ice_cover_stor = sd_ch(ich)%ice_cover_stor - cover_to_jam
          sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor + cover_to_jam
          sd_ch(ich)%icejam_block = sd_ch(ich)%icejam_block + cover_to_jam
              jam_capacity_used_today = jam_capacity_used_today + cover_to_jam

              ice_mobilized_cover_break = min(cover_to_jam_frac_eff * sd_ch(ich)%ice, &
                      prm%mobile_max_daily_frac * sd_ch(ich)%ice)
              ice_mobilized_cover_break = max(0., min(ice_mobilized_cover_break, sd_ch(ich)%ice))
              sd_ch(ich)%ice = sd_ch(ich)%ice - ice_mobilized_cover_break
              sd_ch(ich)%ice_mobile = sd_ch(ich)%ice_mobile + ice_mobilized_cover_break
              sim_ice_thick = sd_ch(ich)%ice / max(ice_area, 1.e-6)
              ice_maturity = sim_ice_thick / max(prm%ice_maturity_ref_thick, 1.e-6)
              ice_ratio = max(0., min(1., ice_maturity))
              mobile_jam_ratio = sd_ch(ich)%ice_mobile / max(ice_cover_max, 1.e-6)
              mobile_pass_jam_ratio = sd_ch(ich)%ice_mobile_pass / max(ice_cover_max, 1.e-6)
              mobile_jam_ratio = max(0., min(1., mobile_jam_ratio))
              mobile_pass_jam_ratio = max(0., min(1., mobile_pass_jam_ratio))
      endif
      endif

      !! v17: block_jam uses only the start-of-day capacity that has not
      !! already been used by cover_to_jam.  It cannot occupy storage capacity
      !! made available by release earlier in the same daily step.
      jam_remaining_capacity_step = max(0., jam_remain_capacity_start - jam_capacity_used_today)
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
          ice_block_material = sd_ch(ich)%ice_mobile + sd_ch(ich)%ice_mobile_pass
      else
          ice_block_material = sd_ch(ich)%ice + sd_ch(ich)%ice_mobile
      endif
      ice_load_block_capacity = block_cap_coeff * ice_block_material * reach%jam_block_modifier
      if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
          ice_load_block_capacity = ice_load_block_capacity * flow_supply_factor
      endif
      ice_load_block_capacity = max(0., ice_load_block_capacity)

      !! Once a jam is present or building, its ability to reduce conveyance is
      !! not only proportional to the daily local ice load.  It also reflects
      !! the presence of an existing jam-storage pool and reach susceptibility.
      !! v21: do NOT scale jam maturity by ice_jam_stor / jam_stor_max, because
      !! jam_stor_max is now a loose safety capacity.  Otherwise increasing
      !! jam_form_stor_max_frac would unintentionally reduce blocking capacity.
      jam_presence_factor = sd_ch(ich)%ice_jam_stor / max(reach%hyd_storage_scale, 1.e-6)
      jam_presence_factor = max(0., min(1., jam_presence_factor))
      jam_material_factor = max(mobile_jam_ratio, mobile_pass_jam_ratio)
      jam_material_factor = max(0., min(1., jam_material_factor))
      !! v38: existing jam storage can enhance constriction only when there
      !! is sufficient ice material to support a jam.  This prevents small
      !! residual jam_stor from maintaining strong blocking during mobile-ice
      !! DRIFT / pass-through periods.
      jam_maturity_factor = max(mobile_jam_ratio, jam_presence_factor * jam_material_factor)
      jam_maturity_factor = max(0., min(1., jam_maturity_factor))

      !! v40: qraw is not used as a hard event trigger, but it should control
      !! how effective a mature jam can be at blocking flow.  At very low
      !! qraw, mobile ice is more likely to drift/pass through without building
      !! strong backwater.  At qraw near or above the reach reference flow, the
      !! full constriction capacity is available.
      !! flow_supply_factor was computed near qraw initialization and is also
      !! used in the BREAKUP day-type decision.
      jam_constriction_capacity = block_frac_max * ht1%flo * reach%jam_susc * &
                                  jam_maturity_factor * flow_supply_factor
      jam_constriction_capacity = max(0., jam_constriction_capacity)

      if (sd_ch(ich)%ice_phase == ICE_BREAKUP) then
          block_capacity = max(ice_load_block_capacity, jam_constriction_capacity)
      else
          block_capacity = ice_load_block_capacity
      endif

      !! 6b. Then block today's incoming water behind the remaining jam capacity.
      if (block_capacity > 1.e-6 .and. jam_remaining_capacity_step > 1.e-6 .and. &
              block_frac_max > 1.e-6 .and. &
              (sd_ch(ich)%ice_phase /= ICE_BREAKUP .or. breakup_day_type == BRK_DAY_BLOCK .or. do_onset_cover_block_today)) then
          blocked_jam = min(block_capacity, block_frac_max * ht1%flo)
          blocked_jam = min(blocked_jam, jam_remaining_capacity_step)
          blocked_jam = max(0., blocked_jam)

          ht1%flo = ht1%flo - blocked_jam
          sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor + blocked_jam
          sd_ch(ich)%icejam_block = sd_ch(ich)%icejam_block + blocked_jam
          jam_capacity_used_today = jam_capacity_used_today + blocked_jam
      endif

      jam_block_capacity = block_capacity
      jam_stor_max_dbg = jam_stor_max
      jam_remaining_capacity = max(0., jam_stor_max - sd_ch(ich)%ice_jam_stor)

      blocked = blocked_cover + blocked_jam

      cover_stor_ratio = sd_ch(ich)%ice_cover_stor / cover_stor_max
      cover_stor_ratio = max(0., min(1., cover_stor_ratio))
      jam_stor_ratio = sd_ch(ich)%ice_jam_stor / max(jam_stor_max, 1.e-6)
      jam_stor_ratio = max(0., min(1., jam_stor_ratio))
      total_ice_stor = sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor

      !! ------------------------------------------------------------------
      !! 7. Final synchronization with the hydrograph used by ch_rtmusk.
      !! ch_rtmusk uses ob(icmd)%tsin(irtstep), so this must be updated
      !! after all release/blocking operations.
      !! ------------------------------------------------------------------
      ht1%flo = max(0., ht1%flo)

      if (time%step == 1) then
          ob(icmd)%tsin(1) = ht1%flo
      else
          tsin_sum = sum(ob(icmd)%tsin)
          raw_flo = max(1.e-6, q_in_rate_raw * 86400.)
          if (tsin_sum <= 1.e-6 .and. ht1%flo > 1.e-6) then
              !! If ice-jam release creates water on a day with zero original
              !! subdaily inflow, scaling a zero hydrograph would still give
              !! zero. Distribute the daily adjusted volume uniformly.
              ob(icmd)%tsin(:) = ht1%flo / real(size(ob(icmd)%tsin))
          else
            adj_ratio = ht1%flo / raw_flo
            adj_ratio = max(0., adj_ratio)
            ob(icmd)%tsin(:) = ob(icmd)%tsin(:) * adj_ratio
          endif
      endif

      !! Final numerical cleanup and diagnostics.
      if (sd_ch(ich)%ice <= ice_state_eps .and. (t_ice_decay > prm%ice_melt_tmp .or. &
          sd_ch(ich)%ice_phase == ICE_WARM)) then
          ht1%flo = ht1%flo + sd_ch(ich)%ice
          sd_ch(ich)%ice = 0.
      endif
      if (sd_ch(ich)%ice_mobile <= mobile_state_eps) sd_ch(ich)%ice_mobile = 0.
      if (sd_ch(ich)%ice_mobile_pass <= mobile_state_eps) sd_ch(ich)%ice_mobile_pass = 0.

      sim_ice_thick = sd_ch(ich)%ice / max(ice_area, 1.e-6)
      ice_maturity = sim_ice_thick / max(prm%ice_maturity_ref_thick, 1.e-6)
      ice_ratio = max(0., min(1., ice_maturity))
      ice_depth_ratio = sim_ice_thick / max(sd_ch(ich)%chd, 1.0e-6)
      ice_depth_ratio = max(0., min(1., ice_depth_ratio))
      mobile_jam_ratio = sd_ch(ich)%ice_mobile / ice_cover_max
      mobile_pass_jam_ratio = sd_ch(ich)%ice_mobile_pass / max(ice_cover_max, 1.e-6)
      if (mobile_jam_ratio < 1.e-4) then
          sd_ch(ich)%ice_mobile = 0.
          mobile_jam_ratio = 0.
      endif
      mobile_jam_ratio = max(0., min(1., mobile_jam_ratio))
      if (mobile_pass_jam_ratio < 1.e-4) mobile_pass_jam_ratio = 0.
      mobile_pass_jam_ratio = max(0., min(1., mobile_pass_jam_ratio))

      total_ice_stor = sd_ch(ich)%ice_cover_stor + sd_ch(ich)%ice_jam_stor

      sd_ch(ich)%icejam_qadj = ht1%flo / 86400.

      !Important: q_prev stores raw inflow before ice-jam adjustment.
      sd_ch(ich)%q_prev = q_in_rate_raw
      
      if (ich == 68) then
          write(9003,*) time%yrc, time%day, ich, &
            "tmax", tmax, "tave", t_air, &
            "phase", sd_ch(ich)%ice_phase, &
            "phase_d", sd_ch(ich)%ice_phase_days, &
            "phase_changed", phase_changed_today, &
            "breakup_onset", breakup_onset_today, &
            "breakup_day_type", breakup_day_type, &
            "release_state", sd_ch(ich)%ice_release_active, &
            "release_days", sd_ch(ich)%ice_release_days, &
            "block_days", sd_ch(ich)%ice_block_days, &
            "form_today", do_jam_formation_today, &
            "qraw", q_in_rate_raw, "qadj", sd_ch(ich)%icejam_qadj, &
            "sim_ice_thick", sim_ice_thick, &
            "ice", sd_ch(ich)%ice, &
            "mobile", sd_ch(ich)%ice_mobile, &
            "mobile_pass", sd_ch(ich)%ice_mobile_pass, &
            "mobile_pass_ratio", mobile_pass_jam_ratio, &
            "cover_stor", sd_ch(ich)%ice_cover_stor, &
            "jam_stor", sd_ch(ich)%ice_jam_stor, &
            "cover_stor_max", cover_stor_max_dbg, &
            "cover_remain_cap", cover_remaining_capacity, &
            "jam_stor_max", jam_stor_max_dbg, &
            "jam_stor_ratio", jam_stor_ratio, &
            "jam_remain_cap_start", jam_remain_capacity_start, &
            "jam_cap_used_today", jam_capacity_used_today, &
            "block_cover", blocked_cover, &
            "cover_to_jam", cover_to_jam, &
            "cover_to_jam_frac", cover_to_jam_frac_eff, &
            "block_jam", blocked_jam, &
            "rel_jam_event", released_jam_event, &
            "rel_jam_leak", released_jam_leak, &
            "release_frac_eff", release_frac_eff, &
            "release_weak_eff", release_weak_eff, &
            "release_ramp", release_ramp_factor, &
            "warm_cleanup_return", warm_cleanup_return, &
            "thaw_weak", thaw_weakening_index, &
            "drift_zone", merge(1, 0, drift_weak_zone), &
            "breakup_onset_thr", prm%breakup_onset_weakening_index, &
            "jam_release_thr", prm%jam_release_weakening_index, &
            "breakup_release_gate", merge(1, 0, breakup_release_gate), &
            "frz_dd", sd_ch(ich)%ice_freeze_dd, &
            "thaw_dd", sd_ch(ich)%ice_thaw_dd, &
            "ros", ros_day, &
            "snow_melt", snow_melt_mm, &
            "susc", jam_susc, &
            "jam_presence_fac", jam_presence_factor, &
            "jam_material_fac", jam_material_factor, &
            "jam_maturity_fac", jam_maturity_factor, &
            "flow_supply_fac", flow_supply_factor, &
            "block_flow_ready", merge(1, 0, block_flow_ready), &
            "mobile_flush", mobile_flushed_by_release, &
            "mobile_drift_pass", mobile_drift_pass, &
            "onset_block_mult", onset_block_mult
      end if

      return 

end subroutine sd_channel_icejam

subroutine sd_channel_ice_advect(j)

      !!    ~ ~ ~ PURPOSE ~ ~ ~
      !!    Advect mobile/broken channel ice to downstream channel objects.
      !!
      !!    Conceptual interpretation:
      !!      sd_ch%ice_mobile is transportable broken ice/floes generated by
      !!      progressive cover deterioration, dynamic/mechanical failure, or
      !!      event-breakup release.  It is an ice-mass state variable, not liquid
      !!      water, and therefore this routine does not modify ht1, ht2, ob%hd,
      !!      ob%tsin, or the liquid-water balance.
      !!
      !!    Updated routing rules:
      !!      1. Mobile ice is transferred along explicit downstream chandeg links.
      !!      2. Downstream reaches capture only part of the incoming mobile ice.
      !!         Capture increases with ice-jam susceptibility, local ice condition,
      !!         and winter/freezing phase.
      !!      3. Captured incoming ice is stored in downstream ice_mobile and can
      !!         participate in local jam formation on the next daily step.
      !!      4. Uncaptured incoming ice is stored in downstream ice_mobile_pass,
      !!         a pass-through pool that is routed on the next daily step without
      !!         directly contributing to local jam formation.
      !!      4. If a downstream object exists but none is a chandeg, unsent mobile
      !!         ice is returned to local stationary ice.  If no downstream object
      !!         exists, mobile ice leaves the represented channel-ice domain.

      use hydrograph_module
      use sd_channel_module
      use channel_module
      use sd_channel_icejam_module

      implicit none

      integer, intent(in) :: j

      integer, parameter :: ICE_WARM       = 0
      integer, parameter :: ICE_FREEZEUP   = 1
      integer, parameter :: ICE_DEEPWINTER = 2
      integer, parameter :: ICE_BREAKUP    = 3

      type(icejam_param_type), save :: prm
      type(icejam_reach_scale_type) :: reach_dn
      logical, save :: prm_initialized = .false.
      real, parameter :: ice_eps = 1.e-4

      integer :: iout
      integer :: iob_dn
      integer :: ich_dn
      integer :: ord_dn

      real :: ice_out
      real :: ice_to_dn
      real :: ice_sent
      real :: ice_unsent
      real :: ice_capture
      real :: ice_pass
      real :: frac_dn
      real :: capture_frac
      real :: capture_capacity
      real :: jam_susc_dn
      real :: ch_vol_cap_dn
      real :: ice_cover_max_dn
      real :: ice_ratio_dn
      real :: sim_ice_thick_dn
      real :: ice_depth_ratio_dn

      logical :: has_downstream_channel
      logical :: has_downstream_object

      if (j <= 0) return

      !! Only pass-through mobile ice is routed by this routine.  Captured/local
      !! ice_mobile remains in the reach as local jam material and does not
      !! participate in same-day downstream routing.  This keeps ice_mobile_pass
      !! equal to previous-day incoming pass-through ice plus any pass-through
      !! ice assigned to this reach by upstream routing.
      if (sd_ch(j)%ice_mobile_pass <= ice_eps) then
            sd_ch(j)%ice_mobile_pass = 0.
            return
      endif

      if (.not. prm_initialized) then
            call icejam_default_params(prm)
            call icejam_validate_params(prm)
            prm_initialized = .true.
      endif

      ice_out = sd_ch(j)%ice_mobile_pass
      ice_sent = 0.
      has_downstream_channel = .false.
      has_downstream_object = .false.

      do iout = 1, ob(icmd)%src_tot

            iob_dn = ob(icmd)%obj_out(iout)
            if (iob_dn <= 0) cycle

            has_downstream_object = .true.

            frac_dn = ob(icmd)%frac_out(iout)
            frac_dn = max(0., min(1., frac_dn))
            if (frac_dn <= ice_eps) cycle

            if (trim(ob(iob_dn)%typ) == "chandeg") then
                  has_downstream_channel = .true.

                  ich_dn = ob(iob_dn)%num
                  if (ich_dn > 0) then
                        ice_to_dn = frac_dn * ice_out

                        !! Avoid sending more than available if frac_out sums
                        !! slightly greater than one.
                        ice_to_dn = min(ice_to_dn, max(0., ice_out - ice_sent))

                        if (ice_to_dn > ice_eps) then
                              !! Downstream ice-jam susceptibility and local ice capacity.
                              !! This uses the same continuous reach-scale function as sd_channel_icejam.
                              ord_dn = sd_ch(ich_dn)%order
                              call icejam_compute_reach_scale(prm, sd_ch(ich_dn)%chw, sd_ch(ich_dn)%chl, &
                                      sd_ch(ich_dn)%chd, sd_ch(ich_dn)%chs, sd_ch(ich_dn)%sinu, &
                                      ch_rcurv(ich_dn)%elev(1)%flo_rate, ch_rcurv(ich_dn)%elev(2)%flo_rate, reach_dn)

                              jam_susc_dn = reach_dn%jam_susc
                              ch_vol_cap_dn = reach_dn%hyd_storage_scale
                              ice_cover_max_dn = max(1.e-6, reach_dn%ice_cap_vol)
                              sim_ice_thick_dn = sd_ch(ich_dn)%ice / max(reach_dn%ice_area, 1.e-6)
                              ice_ratio_dn = sim_ice_thick_dn / max(prm%ice_maturity_ref_thick, 1.e-6)
                              ice_ratio_dn = max(0., min(1., ice_ratio_dn))
                              ice_depth_ratio_dn = sim_ice_thick_dn / max(sd_ch(ich_dn)%chd, 1.e-6)
                              ice_depth_ratio_dn = max(0., min(1., ice_depth_ratio_dn))

                              !! Capture fraction.  Susceptible, ice-covered, and still-
                              !! cold reaches capture more incoming floes.  Warm/open
                              !! reaches pass most of the mobile ice downstream.
                              capture_frac = prm%mobile_capture_base + &
                                      prm%mobile_capture_susc_weight * reach_dn%ice_capture_modifier + &
                                      prm%mobile_capture_ice_weight * ice_ratio_dn + &
                                      prm%mobile_capture_depth_weight * ice_depth_ratio_dn

                              select case (sd_ch(ich_dn)%ice_phase)
                              case (ICE_FREEZEUP)
                                    capture_frac = max(capture_frac, prm%freezeup_capture_min)
                              case (ICE_DEEPWINTER)
                                    capture_frac = max(capture_frac, prm%deepwinter_capture_min)
                              case (ICE_BREAKUP)
                                    capture_frac = max(capture_frac, prm%breakup_capture_min)
                              case default
                                    if (ice_ratio_dn <= 0.05) then
                                          capture_frac = min(capture_frac, prm%warm_capture_max)
                                    endif
                              end select

                              capture_frac = max(prm%mobile_capture_min, min(prm%mobile_capture_max, capture_frac))

                              !! Incoming broken ice remains mobile jam material in the
                              !! downstream reach.  The capture fraction is retained as a
                              !! diagnostic / future residence-time indicator, but we do
                              !! not add captured mobile ice to sd_ch%ice because that would
                              !! inflate stable-cover thickness and distort phase logic.
                              capture_capacity = prm%mobile_capture_capacity_mult * ice_cover_max_dn
                              capture_capacity = max(0., capture_capacity)

                              ice_capture = min(capture_frac * ice_to_dn, capture_capacity)
                              ice_capture = max(0., min(ice_capture, ice_to_dn))
                              ice_pass = ice_to_dn - ice_capture
                              ice_pass = max(0., ice_pass)

                              sd_ch(ich_dn)%ice_mobile = sd_ch(ich_dn)%ice_mobile + ice_capture
                              sd_ch(ich_dn)%ice_mobile_pass = sd_ch(ich_dn)%ice_mobile_pass + ice_pass

                              ice_sent = ice_sent + ice_to_dn
                        endif
                  endif
            endif

            if (ice_sent >= ice_out - ice_eps) exit

      end do

      ice_unsent = max(0., ice_out - ice_sent)

      if (has_downstream_channel) then
            !! Some ice was routed to downstream channel(s).  All routed local
            !! mobile material has left this reach.  Any unsent remainder stays
            !! in the pass-through pool and can try to move again on the next day.
            sd_ch(j)%ice_mobile_pass = ice_unsent
      else
            if (has_downstream_object) then
                  !! There is a downstream object, but it is not an explicit
                  !! chandeg. Keep unsent ice mobile locally rather than returning
                  !! it to stable cover ice, so cover-thickness diagnostics remain clean.
                  sd_ch(j)%ice_mobile_pass = ice_unsent
            else
                  !! Outlet of the represented channel-ice domain. Mobile ice exits
                  !! the modelled ice domain and is not converted to liquid water.
                  sd_ch(j)%ice_mobile_pass = 0.
            endif
      endif

      return
end subroutine sd_channel_ice_advect
