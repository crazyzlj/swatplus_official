module sd_channel_icejam_module

!!    Conceptual river-ice / ice-jam parameter module.
!!
!!    Design principle:
!!      1. Keep the liquid-water ice-jam effect as one conceptual backwater
!!         wedge storage: ht1%flo <-> sd_ch%ice_wedge_stor.
!!      2. Keep solid ice as a separate water-equivalent state variable:
!!         sd_ch%ice_vol.  Melt returns liquid water; mobile ice does not.
!!      3. Separate the seasonal ice regime from dynamic jam episodes:
!!           ice_phase     : OPEN, FREEZEUP, STABLE_COVER, BREAKUP
!!           is_jamming    : current breakup-tail jam build-up episode
!!           is_releasing  : current jam/wedge release episode
!!      4. Jam formation is controlled by mobile-ice supply exceeding the local
!!         transport capacity at a susceptible reach.  This is closer to real
!!         breakup-jam mechanics than triggering a jam directly from temperature.
!!      5. Keep parameter families compact.  Derived reach-scale modifiers are
!!         computed here from channel geometry, slope, sinuosity, and reference flow.

      implicit none

      public :: icejam_param_type
      public :: icejam_reach_scale_type
      public :: icejam_default_params
      public :: icejam_validate_params
      public :: icejam_compute_reach_scale
      public :: icejam_clamp
      public :: icejam_sigmoid

      type :: icejam_param_type

          ! --- Thermal ice growth / decay ---
          real :: ice_maturity_ref_thick = 0.30 ! m; larger -> slower maturity, less capture/major-bg; with stable_ice_thick
          real :: ice_frz_tmp = -1.0 ! degC; higher -> easier ice growth; with ice_growth_coeff/freezing DD
          real :: ice_melt_tmp = 0.0 ! degC; lower -> earlier melt; with ice_decay_coeff/max_melt_frac_ice
          real :: ice_growth_coeff = 0.07 ! m/sqrt(DD); higher -> faster ice_vol, stronger storage; with max growth caps
          real :: max_daily_ice_growth_thick = 0.025 ! m/day; caps daily growth; limits cold-spell ice buildup
          real :: ice_decay_coeff = 0.020 ! m/degC/day; higher -> faster ice loss, earlier OPEN/BREAKUP tail
          real :: ice_freeze_inflow_frac = 0.08 ! frac; caps water converted to ice from daily inflow; with max_freeze_frac_stor
          real :: max_freeze_frac_stor = 0.20 ! frac; caps freezing by channel storage; with ice_freeze_inflow_frac
          real :: max_melt_frac_ice = 0.30 ! frac/day; caps ice meltwater addition and ice_vol decay
          real :: ice_min_vol = 1.0 ! m3; numerical cutoff below which ice/mobile states are cleared
          real :: ice_tail_decay = 0.80 ! frac/day; larger -> faster warm-tail ice/mobile removal; affects OPEN cleanup

          ! --- Freeze/thaw memory and phase gates ---
          real :: freeze_memory = 0.97 ! 0-1; higher -> longer cold memory, later thaw/breakup; with thaw_memory
          real :: thaw_memory = 0.90 ! 0-1; higher -> longer warm memory, earlier weakening/breakup
          real :: freezeup_freeze_dd = 6.0 ! DD; higher -> harder FREEZEUP entry; with freezeup_strong_index
          real :: freezeup_strong_index = 0.35 ! 0-1; higher -> stricter freezeup signal; with freezeup_freeze_dd
          real :: breakup_onset_weakening_index = 0.35 ! 0-1; higher -> later BREAKUP onset; with mechanical_breakup_warm_min
          real :: jam_release_weakening_index = 0.50 ! 0-1; higher -> stricter release weakening; with release_force_ratio
          real :: warm_season_weakening_index = 0.70 ! 0-1; higher -> later OPEN cleanup under warm season
          real :: thaw_tmax_base = 0.0 ! degC; lower -> easier thaw weakening; with thaw_tave_base/ROS threshold
          real :: thaw_tmax_base_ros = -1.0 ! degC; lower -> ROS can weaken ice at colder Tmax
          real :: thaw_tave_base = -2.0 ! degC; lower -> easier thaw memory growth; with thaw_tmax_base
          integer :: freezeup_min_days = 14 ! days; minimum FREEZEUP duration; prevents immediate STABLE switch
          integer :: freezeup_max_days = 45 ! days; maximum FREEZEUP duration; forces STABLE after enough cold season
          integer :: stable_min_days_before_breakup = 45 ! days; minimum STABLE age before BREAKUP; with onset gates
          integer :: stable_max_days_before_breakup = 140 ! days; upper STABLE age guard; prevents frozen phase persistence
          integer :: breakup_min_days_before_open = 3 ! days; minimum BREAKUP tail duration before OPEN
          integer :: breakup_max_days_before_open = 70 ! days; upper BREAKUP duration guard; with force_open_day
          integer :: warm_min_days_before_freezeup = 20 ! days; minimum OPEN age before new FREEZEUP
          integer :: new_ice_year_start_day = 270 ! DOY; resets ice-year memory; affects peak/integrity histories
          integer :: cold_start_freezeup_end_day = 120 ! DOY; early-year cold-start window for FREEZEUP logic
          integer :: breakup_onset_start_day = 45 ! DOY; earliest BREAKUP window only; not stable-cover leakage cutoff
          integer :: breakup_release_start_day = 60 ! DOY; earliest release window; with release_force_ratio
          integer :: breakup_warm_exit_start_day = 120 ! DOY; earliest warm OPEN exit window
          integer :: breakup_tail_end_day = 170 ! DOY; late-season tail cutoff; with force_open_day
          integer :: breakup_force_open_day = 180 ! DOY; hard OPEN cleanup after this day if residual ice persists
          real :: warm_ice_thick = 0.01 ! m; below this ice is warm-tail/open-like; affects cleanup
          real :: freezeup_ice_thick = 0.03 ! m; ice thickness for FREEZEUP readiness
          real :: stable_ice_thick = 0.08 ! m; intact-cover threshold; with maturity/integrity for stable protection
          real :: warm_storage_exit_ratio = 0.03 ! 0-1; lower -> release episode exits only when wedge nearly empty
          real :: flush_thaw_dd = 10.0 ! DD; higher -> harder warm-flush memory; with warm_flush_release_days
          real :: force_open_ice_thick = 0.005 ! m; below this force OPEN after warm/late season
          real :: force_open_integrity = 0.05 ! 0-1; below this force OPEN after warm/late season

          ! --- Structural resistance and hydraulic forcing ---
          real :: integrity_gain_freeze = 0.04 ! per DD; higher -> faster structural ice recovery; capped by max gain
          real :: integrity_loss_thaw = 0.08 ! per DD; higher -> faster structural weakening during thaw
          real :: integrity_loss_ros = 0.015 ! per mm; higher -> more ROS damage to structural integrity
          real :: structural_max_gain = 0.12 ! 1/day; caps structural integrity gain; with integrity_gain_freeze
          real :: structural_max_loss = 0.12 ! 1/day; caps structural integrity loss outside deep winter
          real :: structural_max_loss_deepwinter = 0.06 ! 1/day; lower -> stronger deep-winter integrity persistence
          real :: deepwinter_integrity_floor = 0.70 ! 0-1; higher -> stronger stable-cover protection/major-bg
          real :: surface_weak_loss_thaw = 0.12 ! per DD; higher -> stronger fast surface weakening; affects BREAKUP
          real :: surface_weak_loss_ros = 0.035 ! per mm; higher -> stronger ROS surface weakening; affects winter_drain
          real :: surface_weak_recovery_freeze = 0.10 ! per DD; higher -> faster Iweak recovery in cold spells
          real :: surface_weak_memory = 0.75 ! 0-1; higher -> longer surface-weakening memory
          real :: shock_lambda = 1.0 ! scale; higher -> sharper hydraulic shock response; with force weights
          real :: thermal_force_weight = 0.08 ! weight; higher -> thermal contribution to damage/force index
          real :: damage_force_weight = 1.0 ! weight; higher -> hydraulic damage contribution to force index
          real :: breakup_tmax_base = 0.0 ! degC; lower -> easier warm forcing for breakup
          real :: breakup_tave_base = 0.0 ! degC; lower -> easier mean-temp forcing for breakup
          real :: tave_weight = 0.35 ! 0-1; higher -> more tave weight vs tmax in warm forcing
          real :: breakup_tmax_min = 2.0 ! degC; higher -> stricter warm signal for breakup/warm flush
          real :: breakup_tave_min = 0.0 ! degC; higher -> stricter mean-temp warm signal
          real :: breakup_ros_min = 2.0 ! mm; higher -> stricter ROS contribution to breakup forcing
          real :: alpha_min = 0.05 ! bankfull frac; lower -> less under-ice capacity under weak ice
          real :: alpha_max = 0.45 ! bankfull frac; higher -> more under-ice capacity when ice is weak/open
          real :: breakup_force_ratio_min = 0.85 ! F/R; higher -> harder STABLE mechanical breakup
          real :: release_force_ratio = 1.10 ! F/R; higher -> harder ordinary release trigger
          real :: release_integrity_max = 0.55 ! 0-1; lower -> release only after stronger structural weakening

          ! --- Mobile ice and jam-episode mechanics ---
          real :: thermal_mobile_frac = 0.10 ! frac; higher -> more thermally broken ice to mobile pool
          real :: mechanical_breakup_fr = 0.90 ! F/R; higher -> harder mechanical ice breakup/mobile generation
          real :: mechanical_breakup_fr_scale = 0.70 ! scale; larger -> smoother F/R response for mechanical breakup
          real :: mechanical_mobile_frac = 0.16 ! frac; higher -> more ice_vol converted to mobile ice by breakup
          real :: jam_mobile_min_vol = 100.0 ! m3; higher -> harder jam formation and major-bg memory
          real :: jam_material_min_vol = 500.0 ! m3; higher -> more total ice material required for jam
          real :: ice_transport_cap_base = 5000.0 ! m3/day; higher -> less mobile excess, fewer jams
          real :: ice_transport_cap_min = 100.0 ! m3/day; lower bound on mobile-ice transport capacity
          real :: ice_transport_q_exp = 1.20 ! exp; higher -> transport capacity rises faster with flow
          real :: jam_susc_transport_weight = 0.70 ! 0-1; higher -> susceptible reaches have lower transport capacity
          real :: jam_mobile_excess_capture_frac = 0.70 ! frac; higher -> more mobile excess captured into jam/wedge
          integer :: min_jam_days = 1 ! days; minimum age for major release; ordinary uses stricter >=3
          integer :: max_jam_days = 5 ! days; age for aged-jam leakage; should exceed ordinary min age
          integer :: release_duration_days = 3 ! days; duration of major release burst; ordinary forced to 1 day
          integer :: post_release_lock_days = 5 ! days; ordinary release lockout; with major_post_release_lock_days
          integer :: jam_inactive_max_days = 2 ! days; clears stale jam if no activity persists
          real :: post_release_ice_retention = 0.60 ! frac; higher -> more ice remains after release for tail effects
          real :: post_release_capture_frac = 0.10 ! frac; higher -> more capture allowed after ordinary release
          real :: post_release_leak_mult = 4.00 ! mult; major post-release drainage strength; ordinary uses separate mult
          real :: ordinary_post_release_leak_mult = 1.50 ! mult; higher -> more leakage after ordinary release
          real :: aged_jam_leak_mult = 1.25 ! mult; higher -> more leakage from old active jam; gated by F/R>=1
          real :: freezeup_release_frac = 0.45 ! frac; higher -> more release if freezeup jam fails/breaks
          real :: release_ice_to_mobile_frac = 0.45 ! frac; higher -> more jam ice returns to mobile pool at release
          real :: mobile_resistance_weight = 0.30 ! weight; higher -> mobile ice adds more hydraulic resistance
          real :: reference_ice_vol = 10000.0 ! m3; scaling volume for ice-material nondimensional factors

          ! --- Major jam / warm-flush release gate ---
          real :: major_freeze_dd_min = 20.0 ! DD; higher -> stricter deep-winter prerequisite for major release
          real :: major_ice_maturity_min = 0.45 ! 0-1; higher -> stronger ice storage needed for major background
          real :: major_integrity_peak_min = 0.55 ! 0-1; higher -> requires stronger prior intact cover
          real :: major_snowpack_ante_min = 25.0 ! mm SWE; higher -> requires larger antecedent snow storage
          real :: major_snowpack_peak_min = 40.0 ! mm SWE; higher -> requires larger seasonal snowpack peak
          real :: major_frz_surf_min = 0.35 ! 0-1; higher -> requires stronger frozen-soil state for major event
          real :: major_frz_area_min = 0.35 ! frac; higher -> requires larger frozen area for major event
          real :: major_warm_tave_min = 1.0 ! degC; higher -> stricter warm-flush mean-temp gate
          real :: major_warm_tmax_min = 3.0 ! degC; higher -> stricter warm-flush Tmax gate
          real :: major_snomelt_min = 8.0 ! mm/d; higher -> stricter snowmelt flush gate
          real :: major_ros_min = 2.0 ! mm/d; higher -> stricter ROS flush/drain gate
          real :: major_qrise_min = 0.03 ! bankfull/d; higher -> stricter rising-flow gate
          real :: major_fr_min = 0.80 ! F/R; lower gate in major forcing factor; with start_fr_min
          real :: major_release_max = 0.65 ! frac/day; higher -> stronger first-day major release, less tail
          real :: major_capture_boost = 0.15 ! frac; higher -> more jam capture under major gate before release
          real :: major_release_storage_boost = 0.00 ! frac; higher -> more same-day storage available; kept 0 to avoid full flush
                                                        !V3.13: no same-day full flush
          integer :: warm_flush_release_days = 3 ! days; memory of warm-flush trigger for release window
          real :: warm_flush_memory_base_min = 0.75 ! 0-1; lower -> easier warm-flush memory from background state
          real :: mechanical_breakup_base_min = 0.75 ! 0-1; lower -> easier STABLE mechanical BREAKUP transition
          real :: mechanical_breakup_surface_weak_min = 0.70 ! 0-1; lower -> easier surface-weak breakup trigger
          real :: mechanical_breakup_warm_min = 0.25 ! 0-1; lower -> earlier BREAKUP; major gate not directly changed
          real :: major_release_start_fr_min = 2.50 ! F/R; higher -> delays major release until stronger hydraulics
          real :: major_wedge_ratio_min = 0.15 ! 0-1; higher -> fewer major candidates; filters small-wedge years
          integer :: major_release_pending_days = 3 ! days; major background/trigger memory length
          integer :: major_post_release_lock_days = 10 ! days; lockout after major release; suppresses repeated releases
          real :: major_post_release_capture_frac = 0.03 ! frac; lower -> less post-major recapture/tail release
          real :: ordinary_release_max_frac = 0.30 ! frac; higher -> stronger non-major release; with capacity cap
          real :: ordinary_release_capacity_frac = 0.05 ! frac cap; higher -> larger non-major release scaled by wedge capacity

          ! --- One conceptual liquid backwater-wedge storage ---
          real :: bankfull_storage_min = 1000.0 ! m3; minimum storage scale for wedge capacity in small channels
          real :: wedge_capacity_bankfull_mult = 2.0 ! mult; higher -> larger max wedge storage; affects major volume
          real :: wedge_capture_cover_frac = 0.08 ! frac; higher -> more STABLE ice-cover impoundment
          real :: deepwinter_cover_capture_frac = 0.40 ! frac; higher -> more intact deep-winter capture/storage
          real :: deepwinter_cover_flow_boost = 0.35 ! frac; higher -> more high-flow capture under intact cover
          real :: deepwinter_cover_q_ref = 0.10 ! bankfull frac; lower -> high-flow capture boost saturates earlier
          real :: deepwinter_cover_q_damp = 0.25 ! damp; higher -> suppresses deep-winter capture at high flow
          real :: deepwinter_capture_q_damp_frac = 0.25 ! frac; lower -> weaker second-stage q damping in deep winter
          real :: deepwinter_leak_mult = 0.05 ! mult; lower -> more deep-winter storage retention
          integer :: winter_pulse_drain_days = 3 ! days; length of stable-cover controlled drainage memory
          real :: winter_drain_excess_frac = 0.30 ! frac; higher -> stronger non-event winter drainage of excess flow
          real :: winter_drain_storage_frac = 0.015 ! frac/day; higher -> larger additional winter drainage from wedge
          real :: breakup_jam_k_frac = 0.75 ! frac; higher -> stronger dynamic K/X effect in BREAKUP jam state
          real :: wedge_capture_jam_frac = 0.75 ! frac; higher -> more capture during active jam episode
          real :: wedge_capture_tail_frac = 0.15 ! frac; higher -> more breakup-tail capture without active jam
          real :: wedge_base_capture_frac = 0.02 ! frac; higher -> more background capture in ice season
          real :: wedge_capture_q_damp = 1.5 ! damp; higher -> less capture at high q_rel
          real :: underice_alpha_min = 0.05 ! bankfull frac; lower -> less intact under-ice conveyance, more capture
          real :: underice_alpha_max = 0.35 ! bankfull frac; higher -> more conveyance under weak/open ice, less capture
          real :: wedge_release_storage_frac = 0.25 ! frac; higher -> larger storage-based release tendency
          real :: wedge_release_min = 0.35 ! frac; minimum release fraction for active release episode
          real :: wedge_release_max = 0.85 ! frac; maximum release fraction; bounds major/non-major release
          real :: open_wedge_leak_frac = 0.25 ! frac/day; higher -> faster OPEN cleanup of residual wedge
          real :: tail_wedge_leak_frac = 0.05 ! frac/day; base leakage rate used by stable/breakup multipliers
          real :: stable_unprotected_leak_max_mult = 1.00 ! mult; caps STABLE non-event leakage when cover not intact
          real :: breakup_background_leak_mult = 0.50 ! mult; lower -> less non-event BREAKUP tail leakage
          real :: mobile_wedge_capacity_weight = 0.50 ! weight; higher -> mobile ice increases wedge capacity more
          real :: mobile_wedge_capture_weight = 0.35 ! weight; higher -> mobile ice increases capture efficiency more

          ! --- Muskingum hydraulic modifiers used by ch_rtmusk ---
          integer :: icejam_msk_dynamic = 0 ! flag; 1 updates K/X, 0 volume-only; compare sensitivity carefully
          real :: k_cover_mult = 1.25 ! mult; higher -> slower routing under cover if dynamic K/X enabled
          real :: k_jam_mult = 5.00 ! mult; higher -> stronger jam routing delay if dynamic K/X enabled
          real :: k_release_mult = 1.50 ! mult; higher -> slower release-wave routing if dynamic K/X enabled
          real :: k_tail_max = 1.00 ! mult; upper tail K boost after jam/release
          real :: k_min_mult = 0.50 ! mult; lower bound for dynamic K multiplier
          real :: k_max_mult = 8.00 ! mult; upper bound for dynamic K multiplier
          real :: x_cover = 0.15 ! 0-0.5; Muskingum X under cover; lower -> more attenuation
          real :: x_jam = 0.001 ! 0-0.5; Muskingum X during jam; near zero -> strong attenuation
          real :: x_release = 0.05 ! 0-0.5; Muskingum X during release; affects peak timing/attenuation

          ! --- Downstream mobile-ice capture ---
          real :: mobile_capture_base = 0.10 ! frac; baseline downstream mobile-ice capture
          real :: mobile_capture_susc_weight = 0.55 ! weight; higher -> susceptibility controls mobile capture more
          real :: mobile_capture_ice_weight = 0.20 ! weight; higher -> local ice amount controls capture more
          real :: mobile_capture_depth_weight = 0.15 ! weight; higher -> shallow/depth effect controls capture more
          real :: mobile_capture_min = 0.05 ! frac; lower bound of downstream mobile-ice capture
          real :: mobile_capture_max = 0.90 ! frac; upper bound of downstream mobile-ice capture
          real :: mobile_capture_capacity_mult = 1.50 ! mult; higher -> larger local capacity to retain incoming ice
          real :: warm_capture_max = 0.15 ! frac; max capture in warm/open season; lower -> faster mobile passage
          real :: freezeup_capture_min = 0.35 ! frac; min capture during FREEZEUP; higher -> more ice retention
          real :: stable_capture_min = 0.45 ! frac; min capture during STABLE; higher -> more cover retention
          real :: breakup_capture_min = 0.40 ! frac; min capture during BREAKUP; higher -> more tail jams

          ! --- Reach susceptibility transfer function ---
          real :: slope_ref = 0.001 ! m/m; reference slope for susceptibility; with slope_pow/slope_min
          real :: slope_pow = 0.5 ! exp; higher -> slope contrast affects jam susceptibility more
          real :: slope_min = 1.0e-6 ! m/m; avoids zero-slope division; affects very flat channels
          real :: sinu_alpha = 0.70 ! weight; higher -> sinuosity contributes more to susceptibility
          real :: length_ref = 5.0 ! km; larger -> weaker length effect on susceptibility
          real :: w_slope = 0.35 ! weight; slope component in jam_susc; normalize with other weights
          real :: w_sinu = 0.25 ! weight; sinuosity component in jam_susc
          real :: w_inter = 0.25 ! weight; hydraulic interaction/width-depth component in jam_susc
          real :: w_len = 0.15 ! weight; reach-length component in jam_susc
          real :: width_min = 1.0 ! m; lower bound for geometry scaling; affects small channels
          real :: depth_min = 0.2 ! m; lower bound for depth scaling; affects shallow channels
          real :: jam_susc_min = 0.10 ! 0-1; lower bound for reach jam susceptibility
          real :: jam_susc_max = 0.95 ! 0-1; upper bound for reach jam susceptibility
          real :: jam_storage_base = 0.50 ! base; minimum storage modifier from susceptibility
          real :: jam_storage_weight = 0.50 ! weight; higher -> susceptibility increases wedge capacity more
          real :: jam_block_base = 0.30 ! base; minimum blockage modifier from susceptibility
          real :: jam_block_weight = 0.70 ! weight; higher -> susceptibility reduces transport/increases blockage
          real :: jam_capture_base = 0.20 ! base; minimum capture modifier from susceptibility
          real :: jam_capture_weight = 0.80 ! weight; higher -> susceptibility increases jam capture more
          real :: jam_mech_min = 0.20 ! 0-1; lower bound mechanical weakness modifier
          real :: jam_ref_frac = 0.50 ! frac; reference flow fraction for jam scaling/transport

      end type icejam_param_type

      type :: icejam_reach_scale_type
          real :: ice_area = 1.0
          real :: ice_cap_vol = 1.0
          real :: hyd_storage_scale = 1.0
          real :: q_jam_ref_rate = 0.05
          real :: jam_susc = 0.5
          real :: jam_form_modifier = 0.5
          real :: jam_storage_modifier = 0.75
          real :: jam_block_modifier = 0.65
          real :: ice_capture_modifier = 0.60
          real :: mechanical_weakness_modifier = 0.50
          real :: mobile_order_mult = 1.0
      end type icejam_reach_scale_type

contains

      real function icejam_clamp(x, xmin, xmax)
          real, intent(in) :: x, xmin, xmax
          icejam_clamp = max(xmin, min(xmax, x))
      end function icejam_clamp

      real function icejam_sigmoid(x)
          real, intent(in) :: x
          if (x > 50.0) then
              icejam_sigmoid = 1.0
          else if (x < -50.0) then
              icejam_sigmoid = 0.0
          else
              icejam_sigmoid = 1.0 / (1.0 + exp(-x))
          endif
      end function icejam_sigmoid

      subroutine icejam_default_params(p)
          type(icejam_param_type), intent(out) :: p
          p = icejam_param_type()
      end subroutine icejam_default_params

      subroutine icejam_validate_params(p)
          type(icejam_param_type), intent(in) :: p

          if (p%ice_maturity_ref_thick <= 0.) stop "icejam V3 parameter error: ice_maturity_ref_thick <= 0"
          if (p%ice_growth_coeff < 0. .or. p%ice_decay_coeff < 0.) stop "icejam V3 parameter error: invalid ice growth/decay"
          if (p%max_daily_ice_growth_thick < 0. .or. p%max_freeze_frac_stor < 0.) &
              stop "icejam V3 parameter error: invalid freeze limits"
          if (p%max_melt_frac_ice < 0. .or. p%max_melt_frac_ice > 1.) stop "icejam V3 parameter error: invalid melt limit"
          if (p%freezeup_max_days < p%freezeup_min_days) stop "icejam V3 parameter error: freezeup_max_days < freezeup_min_days"
          if (p%stable_max_days_before_breakup < p%stable_min_days_before_breakup) &
              stop "icejam V3 parameter error: stable max < min"
          if (p%breakup_max_days_before_open < p%breakup_min_days_before_open) stop "icejam V3 parameter error: breakup max < min"
          if (p%breakup_tail_end_day < p%breakup_onset_start_day) stop "icejam V3 parameter error: invalid breakup_tail_end_day"
          if (p%breakup_force_open_day < p%breakup_warm_exit_start_day) &
              stop "icejam V3 parameter error: invalid breakup_force_open_day"
          if (.not. (p%warm_ice_thick <= p%freezeup_ice_thick .and. p%freezeup_ice_thick <= p%stable_ice_thick)) &
              stop "icejam V3 parameter error: inconsistent ice-thickness phase thresholds"
          if (p%alpha_min < 0. .or. p%alpha_max <= p%alpha_min) stop "icejam V3 parameter error: invalid resistance alpha"
          if (p%underice_alpha_min < 0. .or. p%underice_alpha_max <= p%underice_alpha_min) &
              stop "icejam V3 parameter error: invalid under-ice capacity alpha"
          if (p%ice_transport_cap_base <= 0. .or. p%ice_transport_cap_min < 0. .or. p%ice_transport_q_exp <= 0.) &
              stop "icejam V3 parameter error: invalid mobile-ice transport parameters"
          if (p%jam_susc_transport_weight < 0. .or. p%jam_susc_transport_weight > 0.95) &
              stop "icejam V3 parameter error: invalid jam_susc_transport_weight"
          if (p%jam_mobile_excess_capture_frac < 0. .or. p%jam_mobile_excess_capture_frac > 1.) &
              stop "icejam V3 parameter error: invalid jam_mobile_excess_capture_frac"
          if (p%wedge_release_min < 0. .or. p%wedge_release_max > 1. .or. p%wedge_release_min > p%wedge_release_max) &
              stop "icejam V3 parameter error: invalid wedge release fractions"
          if (p%icejam_msk_dynamic /= 0 .and. p%icejam_msk_dynamic /= 1) &
              stop "icejam V3 parameter error: icejam_msk_dynamic must be 0 or 1"
          if (p%major_release_max < p%wedge_release_min .or. p%major_release_max > 0.99) &
              stop "icejam V3 parameter error: invalid major_release_max"
          if (p%warm_flush_release_days < 1) &
              stop "icejam V3 parameter error: warm_flush_release_days < 1"
          if (p%warm_flush_memory_base_min < 0.0 .or. p%warm_flush_memory_base_min > 1.0) &
              stop "icejam V3 parameter error: invalid warm_flush_memory_base_min"
          if (p%mechanical_breakup_base_min < 0.0 .or. p%mechanical_breakup_base_min > 1.0) &
              stop "icejam V3 parameter error: invalid mechanical_breakup_base_min"
          if (p%mechanical_breakup_surface_weak_min < 0.0 .or. &
              p%mechanical_breakup_surface_weak_min > 1.0) &
              stop "icejam V3 parameter error: invalid mechanical_breakup_surface_weak_min"
          if (p%mechanical_breakup_warm_min < 0.0 .or. p%mechanical_breakup_warm_min > 1.0) &
              stop "icejam V3 parameter error: invalid mechanical_breakup_warm_min"
          if (p%major_release_start_fr_min < p%major_fr_min) &
              stop "icejam V3 parameter error: invalid major_release_start_fr_min"
          if (p%major_wedge_ratio_min < 0.0 .or. p%major_wedge_ratio_min > 1.0) &
              stop "icejam V3 parameter error: invalid major_wedge_ratio_min"
          if (p%major_release_pending_days < 0) &
              stop "icejam V3 parameter error: invalid major_release_pending_days"
          if (p%major_post_release_lock_days < p%post_release_lock_days) &
              stop "icejam V3 parameter error: major_post_release_lock_days < post_release_lock_days"
          if (p%major_post_release_capture_frac < 0. .or. p%major_post_release_capture_frac > p%post_release_capture_frac) &
              stop "icejam V3 parameter error: invalid major_post_release_capture_frac"
          if (p%ordinary_release_max_frac <= 0.0 .or. p%ordinary_release_max_frac > p%wedge_release_max) &
              stop "icejam V3 parameter error: invalid ordinary_release_max_frac"
          if (p%ordinary_release_capacity_frac <= 0.0 .or. p%ordinary_release_capacity_frac > 1.0) &
              stop "icejam V3 parameter error: invalid ordinary_release_capacity_frac"
          if (p%major_ice_maturity_min <= 0.0 .or. p%major_ice_maturity_min > 1.0) &
              stop "icejam V3 parameter error: invalid major_ice_maturity_min"
          if (p%deepwinter_cover_q_ref <= 0.0) &
              stop "icejam V3 parameter error: invalid deepwinter_cover_q_ref"
          if (p%deepwinter_capture_q_damp_frac < 0.0 .or. p%deepwinter_capture_q_damp_frac > 1.0) &
              stop "icejam V3 parameter error: invalid deepwinter_capture_q_damp_frac"
          if (p%deepwinter_leak_mult < 0.0 .or. p%deepwinter_leak_mult > 1.0) &
              stop "icejam V3 parameter error: invalid deepwinter_leak_mult"
          if (p%winter_drain_excess_frac < 0.0 .or. p%winter_drain_excess_frac > 1.0) &
              stop "icejam V3 parameter error: invalid winter_drain_excess_frac"
          if (p%winter_drain_storage_frac < 0.0 .or. p%winter_drain_storage_frac > 1.0) &
              stop "icejam V3 parameter error: invalid winter_drain_storage_frac"
          if (p%winter_pulse_drain_days < 1) &
              stop "icejam V3 parameter error: winter_pulse_drain_days < 1"
          if (p%breakup_jam_k_frac < 0.0 .or. p%breakup_jam_k_frac > 1.0) &
              stop "icejam V3 parameter error: invalid breakup_jam_k_frac"
          if (p%post_release_capture_frac < 0. .or. p%post_release_capture_frac > 1.) &
              stop "icejam V3 parameter error: invalid post_release_capture_frac"
          if (p%post_release_leak_mult < 1. .or. p%ordinary_post_release_leak_mult < 1. .or. &
              p%aged_jam_leak_mult < 1.) &
              stop "icejam V3 parameter error: invalid post-release/aged-jam leak multipliers"
          if (p%stable_unprotected_leak_max_mult < p%deepwinter_leak_mult .or. &
              p%breakup_background_leak_mult < 0.0) &
              stop "icejam V3 parameter error: invalid non-event leakage multipliers"
          if (p%freezeup_release_frac < 0. .or. p%freezeup_release_frac > 1.) &
              stop "icejam V3 parameter error: invalid freezeup_release_frac"
          if (p%k_min_mult <= 0. .or. p%k_max_mult < p%k_min_mult) stop "icejam V3 parameter error: invalid K multipliers"
          if (p%x_cover < 0. .or. p%x_cover >= 0.5 .or. p%x_jam < 0. .or. p%x_jam >= 0.5) &
              stop "icejam V3 parameter error: invalid Muskingum X values"
          if (p%slope_ref <= 0. .or. p%slope_pow <= 0. .or. p%length_ref <= 0.) &
              stop "icejam V3 parameter error: invalid susceptibility reference values"
          if (p%jam_susc_min < 0. .or. p%jam_susc_max > 1. .or. p%jam_susc_min >= p%jam_susc_max) &
              stop "icejam V3 parameter error: invalid susceptibility bounds"
          if (p%w_slope < 0. .or. p%w_sinu < 0. .or. p%w_inter < 0. .or. p%w_len < 0.) &
              stop "icejam V3 parameter error: susceptibility weights must be non-negative"
          if (p%w_slope + p%w_sinu + p%w_inter + p%w_len <= 0.) &
              stop "icejam V3 parameter error: at least one susceptibility weight must be positive"
      end subroutine icejam_validate_params

      subroutine icejam_compute_reach_scale(p, ch_width, ch_length_km, ch_depth, ch_slope, &
              ch_sinu, q_rate_low, q_rate_high, scale)
          type(icejam_param_type), intent(in) :: p
          real, intent(in) :: ch_width, ch_length_km, ch_depth, ch_slope, ch_sinu
          real, intent(in) :: q_rate_low, q_rate_high
          type(icejam_reach_scale_type), intent(out) :: scale

          real :: width_eff
          real :: depth_eff
          real :: length_eff
          real :: slope_eff
          real :: sinu_eff
          real :: slope_score
          real :: sinu_score
          real :: inter_score
          real :: len_score
          real :: i_raw
          real :: i_den

          width_eff = max(p%width_min, ch_width)
          depth_eff = max(p%depth_min, ch_depth)
          length_eff = max(1.0, 1000.0 * max(0.001, ch_length_km))
          slope_eff = max(p%slope_min, ch_slope)
          sinu_eff = max(1.0, ch_sinu)

          scale%ice_area = max(1.0, width_eff * length_eff)
          scale%ice_cap_vol = max(1.0e-6, p%ice_maturity_ref_thick * scale%ice_area)
          scale%hyd_storage_scale = max(p%bankfull_storage_min, width_eff * depth_eff * length_eff)

          ! Low slope, high sinuosity, long reaches are more susceptible to
          ! frazil/floe accumulation and ice-jam backwater.  The interaction term
          ! emphasizes low-gradient sinuous reaches.
          slope_score = (p%slope_ref / max(slope_eff, p%slope_min)) ** p%slope_pow
          slope_score = icejam_clamp(slope_score, 0.0, 1.0)
          sinu_score = 1.0 - exp(-p%sinu_alpha * max(0.0, sinu_eff - 1.0))
          sinu_score = icejam_clamp(sinu_score, 0.0, 1.0)
          inter_score = slope_score * sinu_score
          len_score = 1.0 - exp(-max(0.0, ch_length_km) / max(p%length_ref, 1.0e-6))
          len_score = icejam_clamp(len_score, 0.0, 1.0)

          i_raw = p%w_slope * slope_score + p%w_sinu * sinu_score + &
                  p%w_inter * inter_score + p%w_len * len_score
          i_den = max(1.0e-6, p%w_slope + p%w_sinu + p%w_inter + p%w_len)
          scale%jam_susc = p%jam_susc_min + (p%jam_susc_max - p%jam_susc_min) * i_raw / i_den
          scale%jam_susc = icejam_clamp(scale%jam_susc, p%jam_susc_min, p%jam_susc_max)

          scale%jam_form_modifier = scale%jam_susc
          scale%jam_storage_modifier = icejam_clamp(p%jam_storage_base + p%jam_storage_weight * scale%jam_susc, 0.0, 2.0)
          scale%jam_block_modifier = icejam_clamp(p%jam_block_base + p%jam_block_weight * scale%jam_susc, 0.0, 2.0)
          scale%ice_capture_modifier = icejam_clamp(p%jam_capture_base + p%jam_capture_weight * scale%jam_susc, 0.0, 2.0)
          scale%mechanical_weakness_modifier = max(p%jam_mech_min, scale%jam_susc)
          scale%q_jam_ref_rate = p%jam_ref_frac * sqrt(max(q_rate_low, 1.0e-6) * max(q_rate_high, 1.0e-6))
          scale%q_jam_ref_rate = max(scale%q_jam_ref_rate, 0.05)
          scale%mobile_order_mult = icejam_clamp(0.70 + 0.15 * log(max(width_eff, 1.0)), 0.70, 1.40)

      end subroutine icejam_compute_reach_scale

end module sd_channel_icejam_module
