      subroutine sq_snom
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine predicts daily snow accumulation, snowmelt,
!!    rain-on-snow heat exchange, snowpack liquid-water retention/release,
!!    and explicit diurnal refreezing in the snowpack.
!!
!!    Core formulation:
!!    1) precipitation phase, rain-on-snow heat, temperature-index melt,
!!       refreezing, liquid retention, and liquid release are integrated in
!!       the same sub-daily sinusoidal temperature loop;
!!    2) daily precipitation is conserved and may be shifted toward warm
!!       sub-steps only on mixed-phase transition days;
!!    3) TIMP is converted from a daily lag factor to a consistent sub-daily
!!       relaxation coefficient;
!!    4) degree-day melt is an empirical potential driven by the combined
!!       air-snow temperature index; explicit cold-content accounting is used
!!       only for rain-on-snow advective heat;
!!    5) liquid drainage is locked when sub-step air temperature is at or below
!!       freezing, while liquid retention depends on snowpack ripening;
!!    6) observed solar radiation controls daylight degree-day melt
!!       enhancement, while humidity and wind control rain-on-snow heat;
!!    7) daily rain-on-snow intensity can open a bounded preferential-flow
!!       pathway, allowing a fraction of new liquid input to bypass matrix
!!       storage without releasing historical liquid water;
!!    8) snomlt stores generated solid-snow melt, while precip_eff receives
!!       snow-free rainfall plus liquid water released from the snowpack.
!!
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru         |none          |HRU number
!!    snocov1      |none          |1st shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 50% snow cover
!!    snocov2      |none          |2nd shape parameter for snow cover equation
!!                                |This parameter is determined by solving the
!!                                |equation for 95% snow cover
!!    snotmp       |deg C         |temperature of snow pack in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    wst(:)%weat%ts(:)  |mm H2O  |precipitation for the time step during day
!!    snofall      |mm H2O        |amount of precipitation falling as snow on day
!!    snomlt       |mm H2O        |solid snow melted during the day
!!    precip_eff   |mm H2O        |rain bypass plus liquid water released from snowpack
!!    precip_eff_liq |mm H2O       |liquid rainfall directly reaching land surface
!!                                |after canopy interception and snow-cover bypass;
!!                                |excludes snowmelt and snowpack liquid release
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sin, Cos, Exp, Atan, Max, Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use time_module
      use hydrograph_module
      use basin_module, only : bsn_cc
      use hru_module, only : hru, ihru, precip_eff, precip_eff_liq, snofall, snomlt
      use sd_channel_module
      use climate_module, only:  w
      use output_landscape_module
      
      implicit none

      integer :: j = 0                 !none      |HRU number
      integer :: iob = 0               !none      |object number
      integer :: iru = 0               !none      |routing unit containing the HRU
      integer :: icha = 0              !none      |channel index
      integer :: k = 0                 !none      |sub-daily integration counter
      integer, parameter :: n_sub = 24 !none      |number of sub-daily integration points

      real :: smfac = 0.               !mm/degC/d |daily snowmelt factor
      real :: rto_sno = 0.             !none      |snow water / full-cover snow water
      real :: snocov = 0.              !none      |fraction of HRU area covered with snow
      real :: snotmp = 0.              !deg C     |original SWAT+ snowpack temperature proxy
      real :: old_sno_mm = 0.          !mm        |solid snow SWE before adding new snow
      real :: snowpack_swe = 0.        !mm        |total snowpack SWE, solid snow plus retained liquid
      real :: new_snow_temp = 0.       !deg C     |temperature assigned to new snowfall
      real :: sub_timp = 0.            !none      |sub-daily snow temperature lag factor
      real :: sno_before = 0.          !mm        |temporary snow SWE before melt/refreeze
      real :: sno_after = 0.           !mm        |temporary snow SWE after process
      
      real :: precip_day = 0.          !mm        |daily precipitation entering the snow routine
      real :: precip_sub = 0.          !mm        |sub-daily precipitation amount
      real :: rain_sub = 0.            !mm        |sub-daily liquid precipitation
      real :: snow_sub = 0.            !mm        |sub-daily solid precipitation
      real :: rain_mm = 0.             !mm        |daily liquid precipitation after phase partitioning
      real :: snow_mm = 0.             !mm        |daily solid precipitation after phase partitioning
      real :: rain_snow = 0.           !mm        |daily rainfall falling on snow-covered area
      real :: rain_bypass = 0.         !mm        |daily rainfall falling on snow-free area
      real :: rain_snow_sub = 0.       !mm        |sub-daily rainfall falling on snow-covered area
      real :: rain_bypass_sub = 0.     !mm        |sub-daily rainfall falling on snow-free area
      real :: frac_rain = 0.           !none      |sub-daily liquid precipitation weight
      real :: precip_weight = 0.       !none      |sub-daily precipitation distribution weight
      real :: precip_weight_sum = 0.   !none      |normalizing sum of precipitation weights
      real :: precip_bias_strength = 0.!none      |warm-period precipitation bias strength
      real :: warm_factor = 0.         !none      |relative warmth for precipitation weighting
      real :: t_weight = 0.            !deg C     |temporary sub-daily temperature for weighting
      real :: phase_shift = 0.          !rad       |temperature phase angle with delayed daily maximum
      real :: temp_phase_delay = 0.     !hours     |delay of Tmax relative to noon

      real :: melt_dd = 0.             !mm        |daily degree-day melt generated
      real :: melt_ros = 0.            !mm        |daily rain-on-snow heat melt generated
      real :: melt_dd_sub = 0.         !mm        |sub-daily degree-day melt generated
      real :: melt_ros_sub = 0.        !mm        |sub-daily rain-on-snow heat melt generated
      real :: melt_act_sub = 0.        !mm        |sub-daily actual melt limited by solid SWE
      real :: melt_index_sub = 0.      !deg C     |sub-daily positive melt-temperature index
      real :: melt_pot_sub = 0.        !mm        |sub-daily potential melt energy before cold-content correction
      real :: melt_energy_sub = 0.     !mm        |melt-equivalent energy remaining after cold-content correction

      real :: liq_max_ratio = 0.       !none      |liquid-water holding capacity fraction
      real :: liq_capacity = 0.        !mm        |liquid water retained in snowpack
      real :: liq_release = 0.         !mm        |daily liquid water released from snowpack
      real :: liq_release_sub = 0.     !mm        |sub-daily liquid water release
      real :: liq_excess_sub = 0.      !mm        |liquid water above current holding capacity
      real :: liq_theta_min = 0.       !none      |minimum liquid holding ratio for cold/dry snow
      real :: liq_theta_max = 0.       !none      |maximum liquid holding ratio for ripe/wet snow
      real :: liq_ripen_k = 0.         !1/degC    |snow-temperature ripening coefficient
      real :: liq_temp_factor = 0.     !none      |temperature-based ripening factor
      real :: liq_wet_factor = 0.      !none      |liquid-water-based ripening factor
      real :: liq_ripen_state = 0.     !none      |instantaneous ripening state
      real :: liq_ripen_memory = 0.    !none      |one-day nondecreasing ripening memory
      real :: long_term_ripe = 0.      !none      |cross-day snow-temperature ripening proxy
      real :: current_ripe_drive = 0.  !none      |combined short- and long-term ripening driver
      logical :: drain_ready = .false. !true      |warm/ripe snowpack can drain excess liquid

      real :: refreeze = 0.            !mm        |daily liquid water refrozen into snowpack
      real :: refreeze_sub = 0.        !mm        |sub-daily liquid water refrozen
      real :: refreeze_excess_sub = 0. !mm        |excess liquid refrozen during cold-lock retention
      real :: refreeze_pot = 0.        !mm        |sub-daily refreezing potential
      real :: refreeze_factor = 0.     !none      |fraction of melt factor used for refreezing

      real :: tsnow = 0.               !deg C     |all-snow precipitation threshold
      real :: train = 0.               !deg C     |all-rain precipitation threshold
      real :: transition_width = 0.    !deg C     |mixed rain-snow transition width
      real :: train_ros = 0.           !deg C     |sub-daily rain temperature proxy

      real :: pi = 0.                  !none      |pi constant
      real :: t_sub = 0.               !deg C     |sub-daily air temperature
      real :: tave_day = 0.            !deg C     |mean of daily maximum and minimum temperature
      real :: amp_day = 0.             !deg C     |daily sinusoidal temperature amplitude

      real :: cpw_lf = 0.              !mm/degC   |specific heat of water / latent heat of fusion
      real :: cis_lf = 0.              !mm/degC   |specific heat of ice / latent heat of fusion
      real :: rain_heat_melt_eq = 0.   !mm        |rain advective heat as melt-equivalent depth
      real :: cold_content = 0.        !mm        |snowpack cold content as melt-equivalent depth
      real :: heat_excess = 0.         !mm        |rain heat remaining after satisfying cold content
      real :: ros_heat_mult = 0.       !none      |base rain heat multiplier
      real :: ros_heat_mult_sub = 0.   !none      |sub-daily rain heat multiplier
      real :: ros_heat_rain_scale = 0. !mm/substep|liquid rain scale for heat multiplier
      real :: ros_heat_mult_cap = 0.   !none      |weather-regime upper bound of rain heat multiplier
      real :: rain_scale_factor = 0.   !none      |rain amount activation for ROS heat multiplier
      real :: slr_ref = 0.             !MJ/m2/d   |reference spring shortwave radiation
      real :: solrad_eff = 0.          !MJ/m2/d   |non-negative daily solar radiation
      real :: solradmx_eff = 0.        !MJ/m2/d   |non-negative maximum/clear-sky solar radiation
      real :: rad_ratio = 0.           !none      |daily solar radiation ratio
      real :: rad_modifier = 1.        !none      |daily shortwave melt modifier
      real :: rad_sub_modifier = 1.    !none      |daylight-only sub-step shortwave modifier
      real :: rad_min_mult = 0.        !none      |minimum daytime solar melt multiplier
      real :: rad_max_mult = 0.        !none      |maximum daytime solar melt multiplier
      real :: solar_shape = 0.         !none      |diurnal daylight mask for shortwave radiation
      real :: rhum_frac = 0.           !none      |relative humidity as fraction, 0-1
      real :: rhum_ros_thr = 0.        !none      |humidity threshold for turbulent ROS heat
      real :: humid_factor = 0.        !none      |humidity activation for turbulent ROS heat
      real :: wind_func = 0.           !none      |wind activation for turbulent ROS heat
      real :: wind_ros_scale = 0.      !m/s       |wind speed scale for turbulent ROS heat
      real :: latent_multiplier = 1.   !none      |humidity-wind rain-on-snow heat multiplier
      real :: ros_turb_melt_eq = 0.   !mm        |unresolved turbulent/latent ROS melt energy
      real :: ros_turb_melt_per_mm = 0.!mm/mm    |turbulent melt-equivalent per mm ROS rain
      real :: ros_turb_temp_min = 0.  !deg C     |lower air-temperature bound for turbulent ROS heat
      real :: ros_turb_temp_ref = 0.  !deg C     |air-temperature scale for turbulent ROS heat
      real :: ros_turb_temp_factor = 0.!none     |temperature activation for turbulent ROS heat
      real :: liq_input_sub = 0.      !mm        |new liquid input to the snowpack in this sub-step
      real :: pref_flow_frac = 0.     !none      |preferential-flow fraction of new liquid input
      real :: pref_flow_release_sub = 0.!mm      |sub-daily preferential-flow release
      real :: pref_flow_max_frac = 0. !none      |maximum preferential-flow bypass fraction
      real :: pref_flow_input_scale = 0.!mm      |liquid input scale for preferential-flow activation
      real :: pref_flow_rain_scale = 0.!mm       |sub-step ROS rain scale for preferential-flow activation
      real :: pref_flow_rain_day_thr = 0.!mm/day |daily ROS rainfall threshold for preferential-flow opening
      real :: pref_flow_rain_day_scale = 0.!mm/day|daily ROS rainfall scale for preferential-flow opening
      real :: pref_flow_wet_ref = 0.  !none      |liquid-water ratio for preferential-flow wetness activation
      real :: pref_flow_temp_min = 0. !deg C     |minimum air temperature for preferential drainage
      real :: pref_flow_drive = 0.    !none      |liquid-input activation for preferential flow
      real :: pref_flow_day_drive = 0.!none      |daily ROS rainfall activation for preferential flow
      real :: pref_flow_weather = 0.  !none      |humid/windy ROS activation for preferential flow
      real :: pref_flow_wetness = 0. !none      |snowpack liquid-water activation for preferential flow
      real :: pref_flow_ripe = 0.     !none      |ripening activation for preferential flow
      real :: pref_old_liq_frac = 0.  !none      |fraction of above-capacity old liquid released by preferential flow
      real :: pref_old_liq_release = 0.!mm      |old above-capacity liquid released by preferential flow
      real :: liq_excess_pref = 0.    !mm        |above-capacity liquid available before preferential flow
      real :: cold_lock_ros_rain_day_thr = 0.!mm/day |daily ROS rain needed to relax cold-lock
      real :: cold_lock_wet_min = 0.  !none      |minimum wetness/ripeness to relax cold-lock
      real :: cold_lock_temp_min = 0. !deg C     |lower air-temperature bound for ROS cold-lock relaxation
      logical :: pref_flow_event = .false. !true |daily ROS rainfall can open preferential paths
      logical :: ros_unlock_cold_lock = .false.!true |wet ROS can partially relax cold-lock

      logical :: snow_process = .false.      !true if snowpack process occurs during the day
      logical :: ros_event = .false.         !true when rain falls on snow-covered area
      logical :: ros_melt_event = .false.    !true when rain heat produces additional melt
      logical :: cold_lock_release = .false.  !true when cold conditions prevent liquid drainage
      
      j = ihru

      !! Use the original SWAT+ degree-day method when bsn_cc%snom is not enabled.
      !! This branch is intentionally kept close to the official sq_snom.f90
      !! implementation to minimize conflicts when synchronizing with upstream code.
      if (bsn_cc%snom /= 1) then
        hru(j)%sno_liq = 0.
        snotmp = 0.
        precip_eff_liq = precip_eff

        !! estimate snow pack temperature
        snotmp = snotmp * (1. - hru(j)%sno%timp) + w%tave * hru(j)%sno%timp

        if (w%tave <= hru(j)%sno%falltmp) then
          !! calculate snow fall
          hru(j)%sno_mm = hru(j)%sno_mm + precip_eff
          snofall = precip_eff
          precip_eff = 0.
          precip_eff_liq = 0.
          !! set subdaily effective precip to zero
          if (time%step > 1) w%ts = 0.
        endif
 
        if (w%tmax > hru(j)%sno%melttmp .and. hru(j)%sno_mm > 0.) then
          !! adjust melt factor for time of year
          smfac = (hru(j)%sno%meltmx + hru(j)%sno%meltmn) / 2. + Sin((time%day - 81) / 58.09) *     &
                        (hru(j)%sno%meltmx - hru(j)%sno%meltmn) / 2.        !! 365/2pi = 58.09
          snomlt = smfac * (((snotmp + w%tmax)/2.) - hru(j)%sno%melttmp)

          !! adjust for areal extent of snow cover
          if (hru(j)%sno_mm < hru(j)%sno%covmx) then
            rto_sno = hru(j)%sno_mm / hru(j)%sno%covmx
            snocov = rto_sno / (rto_sno + Exp(hru(j)%snocov1 - hru(j)%snocov2 * rto_sno))
          else
            snocov = 1.
          endif
          snomlt = snomlt * snocov
          if (snomlt < 0.) snomlt = 0.
          if (snomlt > hru(j)%sno_mm) snomlt = hru(j)%sno_mm
          hru(j)%sno_mm = hru(j)%sno_mm - snomlt
          precip_eff = precip_eff + snomlt
          if (time%step > 1) then
            w%ts(:) = w%ts(:) + snomlt / time%step
          end if
          if (precip_eff < 0.) precip_eff = 0.
        else
          snomlt = 0.
        end if

        return
      endif

      !!reset daily variables
      snofall = 0.
      snomlt = 0.
      precip_day = precip_eff
      precip_eff = 0.
      precip_eff_liq = 0.

      rain_mm = 0.
      snow_mm = 0.
      rain_snow = 0.
      rain_bypass = 0.
      melt_dd = 0.
      melt_ros = 0.
      liq_release = 0.
      refreeze = 0.
      snow_process = .false.
      ros_event = .false.
      ros_melt_event = .false.
      liq_ripen_memory = 0.

      !! physical conversion factors expressed as melt-equivalent water depth.
      !! cpw_lf = 4.186 / 334 = 0.0125; cis_lf = 2.10 / 334 = 0.0063.
      cpw_lf = 0.0125
      cis_lf = 0.0063

      !! Compact parameter set for the enhanced snow module.  The remaining
      !! coefficients below are structural constants, while these fields are
      !! exposed through snow_parameters and may be updated by calibration.
      transition_width = max(0.1, hru(j)%sno%tband)
      liq_theta_max = max(0.0, hru(j)%sno%liqmx)
      liq_theta_min = min(liq_theta_max, 0.40 * liq_theta_max)
      liq_ripen_k = 0.6
      refreeze_factor = max(0.0, min(1.0, hru(j)%sno%refz))
      temp_phase_delay = 2.0

      !! Hybrid energy modifiers.  Solar radiation is applied only during
      !! daylight sub-steps; humidity and wind activate unresolved turbulent /
      !! latent heat during rain-on-snow events. These internal parameters
      !! should be promoted to formal snow parameters after sensitivity tests.
      slr_ref = 18.0             !! reference spring shortwave radiation, MJ/m2/day
      rad_min_mult = 0.65       !! minimum daylight solar melt multiplier
      rad_max_mult = 1.55       !! maximum daylight solar melt multiplier

      ros_heat_mult = 1.00       !! base multiplier for sensible rain heat
      ros_heat_mult_cap = 3.8   !! upper multiplier under humid/windy ROS conditions
      ros_heat_rain_scale = 0.5  !! sub-step rain scale activating ROS heat multiplier
      rhum_ros_thr = 0.6        !! RH threshold above which turbulent ROS heat grows
      wind_ros_scale = 4.0       !! wind scale for turbulent ROS heat activation, m/s

      !! Independent turbulent/latent ROS melt. This term is not multiplied by
      !! the degree-day temperature index; it represents unresolved condensation,
      !! sensible turbulent exchange, and wind-driven heat transfer during warm,
      !! humid rainfall on snow. It is still constrained by snow cold content.
      ros_turb_melt_per_mm = max(0.0, hru(j)%sno%rosk) !! melt-equivalent mm per mm ROS rain at full activation
      ros_turb_temp_min = -2.0    !! allow near-freezing humid/windy ROS to add energy
      ros_turb_temp_ref = 2.0     !! full temperature activation near 2 C and above

      !! Preferential/macropore flow through ripe or rain-loaded snow. This does
      !! not create water; it allows a bounded fraction of new liquid input to
      !! bypass matrix storage on intense ROS/melt sub-steps.
      pref_flow_max_frac = max(0.0, min(0.95, hru(j)%sno%pfmax))
      pref_flow_input_scale = 0.25
      pref_flow_rain_scale = 1.0
      pref_flow_rain_day_thr = 0.5
      pref_flow_rain_day_scale = 5.0
      pref_flow_wet_ref = 0.025
      pref_flow_temp_min = -0.5
      pref_old_liq_frac = 0.16
      cold_lock_ros_rain_day_thr = max(0.0, hru(j)%sno%clrain)
      cold_lock_wet_min = 0.5
      cold_lock_temp_min = -1.8

      pi = 4.0 * atan(1.0)
      tave_day = 0.5 * (w%tmax + w%tmin)
      amp_day = 0.5 * max(0.0, w%tmax - w%tmin)

      !! Daily solar modifier.  It is converted to a daylight-only sub-step
      !! multiplier inside the integration loop so nighttime warm advection is
      !! not incorrectly amplified by daytime shortwave radiation.
      solrad_eff = max(0.0, w%solrad)
      solradmx_eff = max(0.0, w%solradmx)
      if (solrad_eff > 1.e-6) then
        if (solradmx_eff > 1.e-6) then
          rad_ratio = solrad_eff / max(1.e-6, min(slr_ref, solradmx_eff))
        else
          rad_ratio = solrad_eff / max(1.e-6, slr_ref)
        endif
        rad_modifier = max(rad_min_mult, min(rad_max_mult, rad_ratio))
      else
        rad_modifier = 1.0
      endif

      !! Relative humidity is confirmed to be stored as a 0-1 fraction.
      rhum_frac = max(0.0, min(1.0, w%rhum))
      humid_factor = max(0.0, min(1.0, (rhum_frac - rhum_ros_thr) / max(1.e-6, 1.0 - rhum_ros_thr)))
      wind_func = 1.0 - exp(-max(0.0, w%windsp) / max(1.e-6, wind_ros_scale))
      latent_multiplier = ros_heat_mult + (ros_heat_mult_cap - ros_heat_mult) * humid_factor * wind_func
      latent_multiplier = max(ros_heat_mult, min(ros_heat_mult_cap, latent_multiplier))

      !! Daily rain-on-snow intensity opens preferential-flow pathways only
      !! during meaningful liquid-rain events. The actual bypass remains
      !! bounded primarily by sub-step new liquid input, with only a small
      !! above-capacity stored-liquid release during same-day ROS forcing.
      pref_flow_event = .false.
      pref_flow_day_drive = 0.0

      tsnow = hru(j)%sno%falltmp
      train = hru(j)%sno%falltmp + transition_width

      !! Daily precipitation timing is unknown at daily resolution. A uniform
      !! drizzle assumption creates false cold-hour snowfall on high-DTR spring
      !! storm days. Use a normalized warm-period bias only when the daily
      !! temperature range crosses the rain-snow transition band. This preserves
      !! the exact daily precipitation total while assigning more precipitation
      !! to the warmer part of mixed-phase days.
      if (precip_day > 0. .and. w%tmax > tsnow .and. w%tmin < train) then
        precip_bias_strength = 1.50
      else
        precip_bias_strength = 0.0
      endif

      precip_weight_sum = 0.0
      do k = 1, n_sub
        if (amp_day > 1.e-6) then
          !! v5 phase correction: minimum temperature near 02:00 and maximum
          !! near 14:00, matching the common spring diurnal thermal lag.
          phase_shift = 2.0 * pi * (real(k) - temp_phase_delay) / real(n_sub)
          t_weight = w%tmin + (w%tmax - w%tmin) * 0.5 * (1.0 - cos(phase_shift))
          warm_factor = max(0.0, (t_weight - tave_day) / max(amp_day, 1.e-6))
        else
          warm_factor = 0.0
        endif
        precip_weight_sum = precip_weight_sum + 1.0 + precip_bias_strength * warm_factor
      end do
      if (precip_weight_sum <= 1.e-9) precip_weight_sum = real(n_sub)

      !! Convert the daily TIMP factor into a sub-daily relaxation coefficient.
      !! This preserves the intended daily lag while allowing daytime warming and
      !! nighttime cooling inside the same daily routine.
      if (hru(j)%sno%timp <= 0.) then
        sub_timp = 0.
      else if (hru(j)%sno%timp >= 1.) then
        sub_timp = 1.
      else
        sub_timp = 1.0 - (1.0 - hru(j)%sno%timp) ** (1.0 / real(n_sub))
      endif

      !! Initialize ripening memory from snow temperature and retained liquid.
      !! This preserves metamorphic memory after nocturnal refreezing.
      if (hru(j)%sno_mm > 1.e-6) then
        long_term_ripe = exp(liq_ripen_k * hru(j)%sno_tmp)
        long_term_ripe = max(0.0, min(1.0, long_term_ripe))
        liq_wet_factor = hru(j)%sno_liq / max(1.e-6, hru(j)%sno_mm * liq_theta_max)
        liq_wet_factor = max(0.0, min(1.0, liq_wet_factor))
        liq_ripen_memory = max(long_term_ripe, liq_wet_factor)
      else
        liq_ripen_memory = 0.0
      endif

      !! Seasonal melt factor. The original SWAT+ seasonal structure is retained,
      !! but it is applied to the integrated sub-daily positive-temperature area.
      smfac = (hru(j)%sno%meltmx + hru(j)%sno%meltmn) / 2. +              &
              sin((time%day - 81) / 58.09) *                             &
              (hru(j)%sno%meltmx - hru(j)%sno%meltmn) / 2.
      smfac = max(0.0, smfac)

      !! ================================================================
      !! Full diurnal snow energy and water integration
      !! ================================================================
      do k = 1, n_sub

        if (amp_day > 1.e-6) then
          !! Delay the warmest sub-step by about two hours relative to solar
          !! noon; use the same curve for precipitation and snow processes.
          phase_shift = 2.0 * pi * (real(k) - temp_phase_delay) / real(n_sub)
          t_sub = w%tmin + (w%tmax - w%tmin) * 0.5 * (1.0 - cos(phase_shift))
          warm_factor = max(0.0, (t_sub - tave_day) / max(amp_day, 1.e-6))
        else
          t_sub = tave_day
          warm_factor = 0.0
        endif

        if (precip_day > 0.) then
          precip_weight = 1.0 + precip_bias_strength * warm_factor
          precip_sub = precip_day * precip_weight / precip_weight_sum
        else
          precip_sub = 0.0
        endif

        !! Reset sub-step process variables to avoid carryover between
        !! warm/rain and cold/dry sub-steps.
        melt_ros_sub = 0.
        melt_dd_sub = 0.
        melt_pot_sub = 0.
        melt_energy_sub = 0.
        melt_act_sub = 0.
        refreeze_sub = 0.
        refreeze_excess_sub = 0.
        liq_release_sub = 0.
        pref_flow_release_sub = 0.
        pref_old_liq_release = 0.
        liq_excess_pref = 0.
        liq_input_sub = 0.
        ros_turb_melt_eq = 0.
        drain_ready = .false.

        !! ---- 1. Sub-daily precipitation phase partitioning ----
        if (precip_sub > 0.) then
          if (train <= tsnow + 1.e-6) then
            if (t_sub <= tsnow) then
              frac_rain = 0.0
            else
              frac_rain = 1.0
            endif
          else if (t_sub <= tsnow) then
            frac_rain = 0.0
          else if (t_sub >= train) then
            frac_rain = 1.0
          else
            frac_rain = (t_sub - tsnow) / max(1.e-6, train - tsnow)
          endif
          rain_sub = precip_sub * frac_rain
          snow_sub = precip_sub - rain_sub
        else
          rain_sub = 0.
          snow_sub = 0.
        endif

        rain_mm = rain_mm + rain_sub
        snow_mm = snow_mm + snow_sub

        !! ---- 2. Add snowfall to the solid snowpack ----
        if (snow_sub > 0.) then
          snow_process = .true.
          old_sno_mm = hru(j)%sno_mm
          new_snow_temp = min(t_sub, 0.0)

          if (old_sno_mm > 1.e-6) then
            hru(j)%sno_tmp = (old_sno_mm * hru(j)%sno_tmp +             &
                              snow_sub * new_snow_temp) /               &
                             (old_sno_mm + snow_sub)
          else
            hru(j)%sno_tmp = new_snow_temp
          endif

          hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
          hru(j)%sno_mm = hru(j)%sno_mm + snow_sub
          snofall = snofall + snow_sub
        endif

        !! ---- 3. If no snowpack exists, liquid rain bypasses directly ----
        if (hru(j)%sno_mm <= 1.e-6) then
          if (hru(j)%sno_liq > 0.) then
            precip_eff = precip_eff + hru(j)%sno_liq
            hru(j)%sno_liq = 0.
          endif
          precip_eff = precip_eff + rain_sub
          precip_eff_liq = precip_eff_liq + rain_sub
          cycle
        endif

        snow_process = .true.

        !! ---- 4. Update snow temperature with sub-daily lag ----
        !! Air temperature above freezing can only warm the snowpack toward 0 C;
        !! it cannot increase snow temperature above 0 C.
        hru(j)%sno_tmp = hru(j)%sno_tmp * (1.0 - sub_timp) +             &
                         min(t_sub, 0.0) * sub_timp
        hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)

        !! ---- 5. Snow cover fraction for this sub-step ----
        if (hru(j)%sno%covmx > 1.e-6 .and. hru(j)%sno_mm < hru(j)%sno%covmx) then
          rto_sno = hru(j)%sno_mm / hru(j)%sno%covmx
          rto_sno = max(0.0, rto_sno)
          snocov = rto_sno / (rto_sno + exp(hru(j)%snocov1 -             &
                                            hru(j)%snocov2 * rto_sno))
        else
          snocov = 1.
        endif
        snocov = max(0.0, min(1.0, snocov))

        !! ---- 6. Rainfall split between snow-covered and snow-free areas ----
        rain_snow_sub = rain_sub * snocov
        rain_bypass_sub = rain_sub - rain_snow_sub
        rain_snow = rain_snow + rain_snow_sub
        rain_bypass = rain_bypass + rain_bypass_sub
        precip_eff = precip_eff + rain_bypass_sub
        precip_eff_liq = precip_eff_liq + rain_bypass_sub

        !! ---- 7. Rain-on-snow heat and liquid input ----
        if (rain_snow_sub > 0.) then
          ros_event = .true.
          hru(j)%sno_liq = hru(j)%sno_liq + rain_snow_sub
          liq_input_sub = liq_input_sub + rain_snow_sub

          !! Rain heat is driven only by liquid rain actually falling on the
          !! snow-covered fraction. Do not use total daily precipitation here,
          !! because snowfall should not amplify rain advective heat.
          train_ros = max(0.0, t_sub)

          !! ROS heat multiplier is activated only by liquid rain falling on
          !! snow-covered area.  Humid and windy conditions represent unresolved
          !! turbulent sensible/latent heat exchange; snowfall cannot amplify it.
          rain_scale_factor = 1.0 - exp(-rain_snow_sub / max(1.e-6, ros_heat_rain_scale))
          ros_heat_mult_sub = ros_heat_mult + (latent_multiplier - ros_heat_mult) * &
                              rain_scale_factor
          ros_heat_mult_sub = max(ros_heat_mult, min(ros_heat_mult_cap, ros_heat_mult_sub))
          !! Sensible heat of the rain itself is proportional to rain
          !! temperature, but turbulent/latent ROS heat is an independent energy
          !! term driven by humid, windy conditions and activated by liquid rain.
          ros_turb_temp_factor = (t_sub - ros_turb_temp_min) / max(1.e-6, ros_turb_temp_ref - ros_turb_temp_min)
          ros_turb_temp_factor = max(0.0, min(1.0, ros_turb_temp_factor))
          ros_turb_melt_eq = rain_snow_sub * ros_turb_melt_per_mm * &
                             humid_factor * wind_func * ros_turb_temp_factor * rain_scale_factor

          rain_heat_melt_eq = rain_snow_sub * cpw_lf * train_ros * ros_heat_mult_sub + &
                              ros_turb_melt_eq

          if (rain_heat_melt_eq > 0. .and. hru(j)%sno_mm > 1.e-6) then
            cold_content = max(0.0, -hru(j)%sno_tmp) * hru(j)%sno_mm * cis_lf
            heat_excess = rain_heat_melt_eq - cold_content

            if (heat_excess >= 0.) then
              hru(j)%sno_tmp = 0.0
              melt_ros_sub = min(heat_excess, hru(j)%sno_mm)
              if (melt_ros_sub > 0.) then
                hru(j)%sno_mm = hru(j)%sno_mm - melt_ros_sub
                hru(j)%sno_liq = hru(j)%sno_liq + melt_ros_sub
                liq_input_sub = liq_input_sub + melt_ros_sub
                melt_ros = melt_ros + melt_ros_sub
                snomlt = snomlt + melt_ros_sub
                ros_melt_event = .true.
              endif
            else
              !! Rain heat only partially offsets cold content.
              if (hru(j)%sno_mm > 1.e-6) then
                hru(j)%sno_tmp = - (cold_content - rain_heat_melt_eq) /  &
                                 (hru(j)%sno_mm * cis_lf)
                hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
              endif
            endif
          endif
        endif

        !! ---- 8. Temperature-index melt from positive sub-daily energy ----
        !! Degree-day melt is an empirical potential and is not reduced again
        !! by explicit cold-content subtraction. Cold snow is represented in
        !! the driver through (t_sub + sno_tmp) / 2.
        if (hru(j)%sno_mm > 1.e-6) then
          melt_index_sub = max(0.0, 0.5 * (t_sub + hru(j)%sno_tmp) -      &
                                      hru(j)%sno%melttmp)

          !! Apply shortwave enhancement only during daylight. At night the
          !! modifier returns to 1.0, so warm-air advection can melt snow but is
          !! not suppressed or amplified by daily shortwave radiation. During
          !! daylight, the modifier moves from 1.0 toward the daily radiation
          !! modifier and remains bounded by the daylight solar limits.
          solar_shape = max(0.0, sin(pi * (real(k) - 6.0) / 12.0))
          if (solar_shape > 1.e-6) then
            rad_sub_modifier = 1.0 + (rad_modifier - 1.0) * solar_shape
            rad_sub_modifier = max(rad_min_mult, min(rad_max_mult, rad_sub_modifier))
          else
            rad_sub_modifier = 1.0
          endif

          melt_dd_sub = smfac * melt_index_sub / real(n_sub) * snocov * rad_sub_modifier
          melt_act_sub = min(max(0.0, melt_dd_sub), hru(j)%sno_mm)

          if (melt_act_sub > 0.) then
            hru(j)%sno_mm = hru(j)%sno_mm - melt_act_sub
            hru(j)%sno_liq = hru(j)%sno_liq + melt_act_sub
            liq_input_sub = liq_input_sub + melt_act_sub
            melt_dd = melt_dd + melt_act_sub
            snomlt = snomlt + melt_act_sub
            hru(j)%sno_tmp = 0.0
          endif
        endif

        !! ---- 9. Explicit sub-daily refreezing during below-freezing periods ----
        !! This is the key DTR correction: meltwater generated in warm hours can
        !! refreeze during cold hours even when the daily mean temperature is above 0 C.
        if (hru(j)%sno_liq > 0. .and. hru(j)%sno_mm > 1.e-6 .and. t_sub < 0.) then
          refreeze_pot = refreeze_factor * smfac * max(0.0, -t_sub) / real(n_sub)
          refreeze_sub = min(hru(j)%sno_liq, max(0.0, refreeze_pot))

          if (refreeze_sub > 0.) then
            hru(j)%sno_liq = hru(j)%sno_liq - refreeze_sub
            hru(j)%sno_mm = hru(j)%sno_mm + refreeze_sub
            refreeze = refreeze + refreeze_sub

            !! Latent heat released by refreezing warms the snowpack toward 0 C.
            if (hru(j)%sno_mm > 1.e-6) then
              hru(j)%sno_tmp = hru(j)%sno_tmp +                           &
                               refreeze_sub / (hru(j)%sno_mm * cis_lf)
              hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
            endif
          endif
        endif

        !! ---- 10. Dynamic liquid-water retention and release ----
        if (hru(j)%sno_mm > 1.e-6) then
          !! Update ripening from liquid-water connectivity and snow temperature.
          !! Release still requires warm, wet, or ripe conditions.
          long_term_ripe = exp(liq_ripen_k * hru(j)%sno_tmp)
          long_term_ripe = max(0.0, min(1.0, long_term_ripe))

          liq_wet_factor = hru(j)%sno_liq / max(1.e-6, hru(j)%sno_mm * liq_theta_max)
          liq_wet_factor = max(0.0, min(1.0, liq_wet_factor))

          current_ripe_drive = max(liq_wet_factor, long_term_ripe)
          current_ripe_drive = max(0.0, min(1.0, current_ripe_drive))
          liq_ripen_memory = max(liq_ripen_memory, current_ripe_drive)

          liq_max_ratio = liq_theta_min + (liq_theta_max - liq_theta_min) * &
                          liq_ripen_memory
          liq_max_ratio = max(liq_theta_min, min(liq_theta_max, liq_max_ratio))

          liq_capacity = hru(j)%sno_mm * liq_max_ratio

          !! Preferential/macropore flow. A meaningful daily ROS event can
          !! open connected flow paths through wet or ripe snow. The primary
          !! bypass is limited to a bounded fraction of same-sub-step new liquid
          !! input. A small additional fraction of above-capacity stored liquid
          !! can drain during mature wet ROS events; this targets large-rain-on-
          !! snowpack events without allowing full snowpack dump.
          pref_flow_event = (rain_snow >= pref_flow_rain_day_thr)
          if (liq_input_sub > 1.e-9 .and. hru(j)%sno_liq > 1.e-9 .and. &
              pref_flow_event .and. &
              (t_sub > pref_flow_temp_min .or. rain_snow_sub > 0.0)) then
            pref_flow_drive = 1.0 - exp(-liq_input_sub / max(1.e-6, pref_flow_input_scale))
            pref_flow_day_drive = 1.0 - exp(-rain_snow / max(1.e-6, pref_flow_rain_day_scale))
            pref_flow_weather = humid_factor * wind_func * &
                                 (1.0 - exp(-rain_snow_sub / max(1.e-6, pref_flow_rain_scale)))
            pref_flow_wetness = hru(j)%sno_liq / max(1.e-6, hru(j)%sno_mm)
            pref_flow_wetness = max(0.0, min(1.0, pref_flow_wetness / max(1.e-6, pref_flow_wet_ref)))
            pref_flow_ripe = max(liq_ripen_memory, max(pref_flow_wetness, pref_flow_weather))
            pref_flow_frac = pref_flow_max_frac * pref_flow_drive * pref_flow_day_drive * pref_flow_ripe
            pref_flow_frac = max(0.0, min(pref_flow_max_frac, pref_flow_frac))
            liq_excess_pref = max(0.0, hru(j)%sno_liq - liq_capacity)
            pref_old_liq_release = pref_old_liq_frac * pref_flow_day_drive * &
                                   pref_flow_ripe * liq_excess_pref
            pref_flow_release_sub = min(hru(j)%sno_liq, pref_flow_frac * liq_input_sub + &
                                        pref_old_liq_release)

            if (pref_flow_release_sub > 0.) then
              hru(j)%sno_liq = hru(j)%sno_liq - pref_flow_release_sub
              liq_release = liq_release + pref_flow_release_sub
              precip_eff = precip_eff + pref_flow_release_sub
            endif
          endif

          liq_excess_sub = hru(j)%sno_liq - liq_capacity

          if (liq_excess_sub > 0.) then

            !! Cold-lock and drainage-gate rules:
            !! 1) sub-freezing air normally locks drainage and favors refreezing;
            !! 2) mature wet snow under active ROS forcing can partially relax
            !!    the cold-lock even when t_sub is slightly below 0 C. This is
            !!    needed for near-freezing large-rain-on-snowpack events;
            !! 3) drainage still requires warm, wet, or ripe forcing.
            ros_unlock_cold_lock = (rain_snow >= cold_lock_ros_rain_day_thr .and. &
                                    rain_snow_sub > 0.0 .and.                    &
                                    t_sub > cold_lock_temp_min .and.              &
                                    max(liq_wet_factor, liq_ripen_memory) >= cold_lock_wet_min)
            cold_lock_release = (t_sub <= 0.0 .and. .not. ros_unlock_cold_lock)
            drain_ready = (.not. cold_lock_release) .and.                 &
                          (t_sub > 0.0 .or. rain_snow_sub > 0.0 .or.       &
                           melt_act_sub > 0.0 .or. liq_wet_factor > 0.75 .or. &
                           ros_unlock_cold_lock)

            if (cold_lock_release) then
              refreeze_pot = max(0.0, refreeze_factor * smfac * max(0.0, -t_sub) / real(n_sub)) + &
                             max(0.0, -hru(j)%sno_tmp) * hru(j)%sno_mm * cis_lf
              refreeze_excess_sub = min(liq_excess_sub, max(0.0, refreeze_pot))

              if (refreeze_excess_sub > 0.) then
                hru(j)%sno_liq = hru(j)%sno_liq - refreeze_excess_sub
                hru(j)%sno_mm = hru(j)%sno_mm + refreeze_excess_sub
                refreeze = refreeze + refreeze_excess_sub

                if (hru(j)%sno_mm > 1.e-6) then
                  hru(j)%sno_tmp = hru(j)%sno_tmp +                       &
                                   refreeze_excess_sub / (hru(j)%sno_mm * cis_lf)
                  hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
                endif
              endif

              liq_release_sub = 0.

            else if (drain_ready) then
              liq_release_sub = liq_excess_sub
              hru(j)%sno_liq = liq_capacity
            else
              liq_release_sub = 0.
            endif

          else
            liq_release_sub = 0.
          endif

          if (liq_release_sub > 0.) then
            liq_release = liq_release + liq_release_sub
            precip_eff = precip_eff + liq_release_sub
          endif
        else
          if (hru(j)%sno_liq > 0.) then
            precip_eff = precip_eff + hru(j)%sno_liq
            liq_release = liq_release + hru(j)%sno_liq
            hru(j)%sno_liq = 0.
          endif
          hru(j)%sno_mm = 0.
          hru(j)%sno_tmp = 0.
        endif

      end do

      !! Clear negligible snowpack and release remaining liquid water.
      if (hru(j)%sno_mm < 1.e-6) then
        if (hru(j)%sno_liq > 0.) then
          precip_eff = precip_eff + hru(j)%sno_liq
          liq_release = liq_release + hru(j)%sno_liq
        endif
        hru(j)%sno_mm = 0.
        hru(j)%sno_tmp = 0.
        hru(j)%sno_liq = 0.
      else
        hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
      endif

      if (precip_eff < 0.) precip_eff = 0.
      if (precip_eff_liq < 0.) precip_eff_liq = 0.
      precip_eff_liq = min(precip_eff_liq, precip_eff)

      if (bsn_cc%snom == 1 .and. bsn_cc%icejam == 1) then

      !! Update rain-on-snow and snowmelt diagnostics for the downstream channel.
      !! snow_melt_m3 stores generated solid-snow melt volume. The channel-scale
      !! mm/day value should be computed later as:
      !! snow_melt = snow_melt_m3 / (upstream_area_ha * 10).
      !! ros_water_m3 stores liquid water relevant to ice-jam forcing:
      !! rain falling on snow-covered area plus snowpack liquid release.
      iob = j + sp_ob1%hru - 1
      if (ob(iob)%ru_tot > 0) then
        iru = ob(iob)%ru(1)                     ! lsu number
        iru = iru + sp_ob1%ru - 1
        if (ob(iru)%src_tot > 0) then
          if (ob(iru)%obtyp_out(1) == 'sdc') then
            icha = ob(iru)%obj_out(1)          ! channel object index
            icha = icha - sp_ob1%chandeg + 1  ! sequential channel index
            if (icha > 0 .and. icha <= size(sd_ch)) then
              snowpack_swe = hru(j)%sno_mm + hru(j)%sno_liq
              if (snowpack_swe > 0.) then
                sd_ch(icha)%snowpack_m3 = sd_ch(icha)%snowpack_m3 + snowpack_swe * hru(j)%area_ha * 10.
                sd_ch(icha)%snowpack_area_ha = sd_ch(icha)%snowpack_area_ha + hru(j)%area_ha
              endif
              !! Phase-change snowmelt volume: mm * ha * 10 = m3.
              !! snomlt is solid snow converted to liquid water; it is not
              !! necessarily snowpack liquid outflow.
              if (snomlt > 0.) then
                sd_ch(icha)%snow_melt_area_ha = sd_ch(icha)%snow_melt_area_ha + hru(j)%area_ha
                sd_ch(icha)%snow_melt_m3 = sd_ch(icha)%snow_melt_m3 + snomlt * hru(j)%area_ha * 10.
              endif

              !! Liquid-water forcing for ice-jam weakening and hydraulic shock:
              !! rain_snow is liquid rain intercepted by snow cover;
              !! liq_release is liquid water released from the snowpack.
              if ((rain_snow + liq_release) > 0.) then
                sd_ch(icha)%ros_water_m3 = sd_ch(icha)%ros_water_m3 + &
                                           (rain_snow + liq_release) * hru(j)%area_ha * 10.
              endif

              if (ros_event) then
                sd_ch(icha)%ros_area_ha = sd_ch(icha)%ros_area_ha + hru(j)%area_ha
              endif
            endif
          endif
        endif
      endif

      endif

      if (time%step > 1 .and. snow_process) then
        w%ts(:) = precip_eff / time%step
      end if

      return
      end subroutine sq_snom