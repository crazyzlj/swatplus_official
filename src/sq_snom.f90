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
!!    6) snomlt stores generated solid-snow melt, while precip_eff receives
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
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sin, Exp, Atan, Max, Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use time_module
      use hydrograph_module
      use hru_module, only : hru, ihru, precip_eff, snofall, snomlt
      use sd_channel_module
      use climate_module, only:  w
      use output_landscape_module
      
      implicit none

      integer :: j = 0                 !none      |HRU number
      integer :: iru = 0               !none      |routing unit containing the HRU
      integer :: icha = 0              !none      |channel index
      integer :: k = 0                 !none      |sub-daily integration counter
      integer, parameter :: n_sub = 24 !none      |number of sub-daily integration points

      real :: smfac = 0.               !mm/degC/d |daily snowmelt factor
      real :: rto_sno = 0.             !none      |snow water / full-cover snow water
      real :: snocov = 0.              !none      |fraction of HRU area covered with snow
      real :: old_sno_mm = 0.          !mm        |snow SWE before adding new snow
      real :: new_snow_temp = 0.       !deg C     |temperature assigned to new snowfall
      real :: sub_timp = 0.            !none      |sub-daily snow temperature lag factor
      
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

      real :: liq_max_ratio = 0.       !none      |liquid-water holding capacity fraction
      real :: liq_capacity = 0.        !mm        |liquid water retained in snowpack
      real :: liq_release = 0.         !mm        |daily liquid water released from snowpack
      real :: liq_release_sub = 0.     !mm        |sub-daily liquid water release
      real :: liq_excess_sub = 0.      !mm        |liquid water above current holding capacity
      real :: liq_theta_min = 0.       !none      |minimum liquid holding ratio for cold/dry snow
      real :: liq_theta_max = 0.       !none      |maximum liquid holding ratio for ripe/wet snow
      real :: liq_ripen_k = 0.         !1/degC    |snow-temperature ripening coefficient
      real :: liq_wet_factor = 0.      !none      |liquid-water-based ripening factor
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
      real :: ros_heat_mult_max = 0.   !none      |upper bound of rain heat multiplier
      real :: ros_heat_rain_scale = 0. !mm/substep|liquid rain scale for heat multiplier

      logical :: snow_process = .false.      !true if snowpack process occurs during the day
      logical :: ros_event = .false.         !true when rain falls on snow-covered area
      logical :: ros_melt_event = .false.    !true when rain heat produces additional melt
      logical :: cold_lock_release = .false.  !true when cold conditions prevent liquid drainage
      
      j = ihru

      !!reset daily variables
      snofall = 0.
      snomlt = 0.
      precip_day = precip_eff
      precip_eff = 0.

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

      !! Internal parameters introduced in this routine. These should ideally be
      !! moved to basin/HRU snow parameters after sensitivity testing.
      transition_width = 2.0
      liq_theta_min = 0.03
      liq_theta_max = 0.10
      liq_ripen_k = 0.5
      refreeze_factor = 0.50
      temp_phase_delay = 2.0

      !! Rain advective heat must be driven by actual liquid rain falling on
      !! snow, not by total precipitation. Keep the multiplier conservative;
      !! it represents unresolved condensation/turbulent heat effects only.
      ros_heat_mult = 1.00
      ros_heat_mult_max = 1.25
      ros_heat_rain_scale = 0.50

      pi = 4.0 * atan(1.0)
      tave_day = 0.5 * (w%tmax + w%tmin)
      amp_day = 0.5 * max(0.0, w%tmax - w%tmin)
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
          !! Thermal phase correction: minimum temperature near 02:00 and
          !! maximum temperature near 14:00.
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
        melt_act_sub = 0.
        refreeze_sub = 0.
        refreeze_excess_sub = 0.
        liq_release_sub = 0.
        drain_ready = .false.

        !! ---- 1. Sub-daily precipitation phase partitioning ----
        if (precip_sub > 0.) then
          if (t_sub <= tsnow) then
            frac_rain = 0.0
          else if (t_sub >= train) then
            frac_rain = 1.0
          else if (train > tsnow) then
            frac_rain = (t_sub - tsnow) / (train - tsnow)
          else
            if (t_sub > tsnow) then
              frac_rain = 1.0
            else
              frac_rain = 0.0
            endif
          endif
          frac_rain = max(0.0, min(1.0, frac_rain))
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

        !! ---- 7. Rain-on-snow heat and liquid input ----
        if (rain_snow_sub > 0.) then
          ros_event = .true.
          hru(j)%sno_liq = hru(j)%sno_liq + rain_snow_sub

          !! Rain heat is driven only by liquid rain actually falling on the
          !! snow-covered fraction. Do not use total daily precipitation here,
          !! because snowfall should not amplify rain advective heat.
          train_ros = max(0.0, t_sub)
          ros_heat_mult_sub = ros_heat_mult + (ros_heat_mult_max - ros_heat_mult) * &
                              (1.0 - exp(-rain_snow_sub / max(1.e-6, ros_heat_rain_scale)))
          ros_heat_mult_sub = max(ros_heat_mult, min(ros_heat_mult_max, ros_heat_mult_sub))
          rain_heat_melt_eq = rain_snow_sub * cpw_lf * train_ros * ros_heat_mult_sub

          if (rain_heat_melt_eq > 0. .and. hru(j)%sno_mm > 1.e-6) then
            cold_content = max(0.0, -hru(j)%sno_tmp) * hru(j)%sno_mm * cis_lf
            heat_excess = rain_heat_melt_eq - cold_content

            if (heat_excess >= 0.) then
              hru(j)%sno_tmp = 0.0
              melt_ros_sub = min(heat_excess, hru(j)%sno_mm)
              if (melt_ros_sub > 0.) then
                hru(j)%sno_mm = hru(j)%sno_mm - melt_ros_sub
                hru(j)%sno_liq = hru(j)%sno_liq + melt_ros_sub
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
          melt_dd_sub = smfac * melt_index_sub / real(n_sub) * snocov
          melt_act_sub = min(max(0.0, melt_dd_sub), hru(j)%sno_mm)

          if (melt_act_sub > 0.) then
            hru(j)%sno_mm = hru(j)%sno_mm - melt_act_sub
            hru(j)%sno_liq = hru(j)%sno_liq + melt_act_sub
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

          liq_excess_sub = hru(j)%sno_liq - liq_capacity

          if (liq_excess_sub > 0.) then

            !! Cold-lock and drainage-gate rules:
            !! 1) if the current air temperature is at or below freezing,
            !!    surface crusting immediately locks drainage; liquid above
            !!    capacity is first refrozen and otherwise retained;
            !! 2) even under non-cold conditions, drainage requires a warm or
            !!    wet/ripe forcing. This avoids releasing water only because a
            !!    diagnostic capacity changed.
            cold_lock_release = (t_sub <= 0.0)
            drain_ready = (.not. cold_lock_release) .and.                 &
                          (t_sub > 0.0 .or. rain_snow_sub > 0.0 .or.       &
                           melt_act_sub > 0.0 .or. liq_wet_factor > 0.80)

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

      !! Update rain-on-snow and snowmelt diagnostics for the downstream channel.
      !! snow_melt_m3 stores generated solid-snow melt volume. The channel-scale
      !! mm/day value should be computed later as:
      !! snow_melt = snow_melt_m3 / (upstream_area_ha * 10).
      if (ob(j)%ru_tot > 0) then
        iru = ob(j)%ru(1)                     ! lsu number
        iru = iru + sp_ob1%ru - 1
        if (ob(iru)%src_tot > 0) then
          if (ob(iru)%obtyp_out(1) == 'sdc') then
            icha = ob(iru)%obj_out(1)          ! channel object index
            icha = icha - sp_ob1%chandeg + 1  ! sequential channel index
            if (icha > 0 .and. icha <= size(sd_ch)) then
              sd_ch(icha)%snow_melt_area_ha = sd_ch(icha)%snow_melt_area_ha + ob(j)%area_ha
              sd_ch(icha)%snow_melt_m3 = sd_ch(icha)%snow_melt_m3 + snomlt * ob(j)%area_ha * 10.
              if (ros_event) then
                sd_ch(icha)%ros_area_ha = sd_ch(icha)%ros_area_ha + ob(j)%area_ha
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