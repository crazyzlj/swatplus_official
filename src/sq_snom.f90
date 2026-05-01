      subroutine sq_snom
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine predicts daily snow accumulation, snowmelt,
!!    rain-on-snow advection, snowpack cold-content compensation,
!!    liquid-water retention/release, and simple refreezing in the snowpack.

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
!!    wst(:)%weat%ts(:)  |mm H2O        |precipitation for the time step during day
!!    snofall      |mm H2O        |amount of precipitation falling as freezing rain/snow on day
!!    snomlt       |mm H2O        |amount of water in snow melt for the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Real, Sin, Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use time_module
      use hydrograph_module
      use hru_module, only : hru, ihru, precip_eff, snofall, snomlt 
      use climate_module, only:  w
      use output_landscape_module
      
      implicit none

      integer :: j = 0      !none       |HRU number
      real :: smfac = 0.    !           |
      real :: rto_sno  = 0. !none       |ratio of current day's snow water to minimum amount needed to
                            !           |cover ground completely 
      real :: snocov = 0.   !none       |fraction of HRU area covered with snow
      real :: snotmp = 0.   !deg C      |temperature of snow pack
      
      real :: precip_day = 0.  !mm         |daily effective precipitation entering the snow routine
      real :: rain_mm = 0.     !mm         |liquid fraction of daily precipitation after phase partitioning
      real :: snow_mm = 0.     !mm H2O     |solid fraction of daily precipitation after phase partitioning
      real :: rain_snow = 0.   !mm         |rainfall falling on the snow-covered fraction of the HRU
      real :: rain_bypass = 0. !mm         |rainfall falling on the snow-free fraction of the HRU
      real :: frac_rain = 0.   !none       |fraction of daily precipitation treated as rain
      real :: frac_snow = 0.   !none       |fraction of daily precipitation treated as snowfall

      real :: melt_dd = 0.     !mm H2O     |degree-day snowmelt before limiting by available snow
      real :: melt_ros = 0.    !mm H2O     |additional melt from rain-on-snow advective heat
      real :: melt_pot = 0.    !mm H2O     |potential total melt from degree-day and rain-on-snow terms
      real :: melt_act = 0.    !mm H2O     |actual melt limited by available solid snow water equivalent

      real :: liq_max_ratio = 0. !none    |maximum liquid-water holding capacity as a fraction of solid snow SWE
      real :: liq_capacity = 0.  !mm      |maximum liquid water that can be retained in the snowpack
      real :: liq_release = 0.   !mm      |liquid water released from the snowpack to runoff/infiltration
      real :: total_liquid_input = 0. !mm |rain-on-snow plus actual melt entering snowpack liquid storage
      real :: refreeze = 0.      !mm H2O  |liquid water refrozen into the solid snow/ice store

      real :: tsnow = 0.       !deg C     |temperature threshold for all-snow precipitation
      real :: train = 0.       !deg C     |temperature threshold for all-rain precipitation
      real :: train_ros = 0.   !deg C     |estimated temperature of rain falling on snow
      real :: dt_range = 0.    !deg C     |daily temperature range used for rain/snow partitioning

      real :: cpw_lf = 0.      !mm/deg C |specific heat of water / latent heat of fusion
                               !         |= 4.186 / 334 = 0.01253 mm melt equivalent per mm rain per deg C
      real :: cis_lf = 0.      !mm/deg C |specific heat of ice / latent heat of fusion
                               !         |= 2.10 / 334 = 0.00629 mm melt equivalent per mm snow per deg C
      
      real :: rain_heat_melt_eq = 0. !mm  |rain advective heat expressed as melt-equivalent water depth
      real :: cold_content = 0.      !mm  |snowpack cold content expressed as melt-equivalent water depth
      real :: heat_excess = 0.       !mm  |rain heat remaining for melt after satisfying cold content
      logical :: snow_process = .false.
      
      j = ihru

      !!reset daily variables
      snofall = 0.
      snomlt = 0.
      melt_dd = 0.
      melt_ros = 0.
      melt_pot = 0.
      melt_act = 0.
      liq_release = 0.
      refreeze = 0.
      rain_bypass = 0.
      rain_snow = 0.
      rain_mm = 0.
      snow_mm = 0.
      frac_rain = 0.
      frac_snow = 0.
      train_ros = 0.
      snow_process = .false.
      
      precip_day = precip_eff
      precip_eff = 0.
      
      cpw_lf = 0.0125
      cis_lf = 0.0063
      
      liq_max_ratio = 0.05 !todo, put this parameter in parameters.bsn
      
      !! calculate snow fall
      !use two thresholds for parcipitation phase partitioning
      tsnow = hru(j)%sno%falltmp
      train = hru(j)%sno%falltmp + 2.0
      
      if (precip_day > 0.) then
        if (w%tmax <= tsnow) then
          frac_snow = 1.0
        else if (w%tmin >= train) then
          frac_snow = 0.0
        else
          dt_range = w%tmax - w%tmin
          if (dt_range > 1.e-6) then
            !! Approximate fraction of the day below the snow threshold,
            !! assuming a linear daily temperature transition.
            frac_rain = (w%tmax - tsnow) / dt_range
            frac_rain = max(0.0, min(1.0, frac_rain))
            frac_snow = 1.0 - frac_rain
          else
            !! Fallback for zero daily temperature range.
            if (w%tave <= tsnow) then
              frac_snow = 1.0
            else if (w%tave >= train) then
              frac_snow = 0.0
            else
              frac_snow = (train - w%tave) / (train - tsnow)
              frac_snow = max(0.0, min(1.0, frac_snow))
            endif
          endif
        endif

        snow_mm = precip_day * frac_snow
        rain_mm = precip_day - snow_mm

        if (snow_mm > 0.) then
          hru(j)%sno_mm = hru(j)%sno_mm + snow_mm
          snofall = snow_mm
        endif

        !temporary value here. If a snowpack exists below, this will be overwritten
        !by rain_bypass + liquid water released from the snowpack.
        precip_eff = rain_mm
      else
        precip_eff = 0.
      endif
      
      !! process existing snowpack
      if (hru(j)%sno_mm > 0.) then
        snow_process = .true.
        !estimate internal snow pack temperature, keep it <= 0 deg C.
        !positive energy contributes to ripening/melt rather than raising snow temperature above freezing.
        snotmp = hru(j)%sno_tmp
        snotmp = snotmp * (1. - hru(j)%sno%timp) + w%tave * hru(j)%sno%timp
        snotmp = min(snotmp, 0.0)
        hru(j)%sno_tmp = snotmp
        
        !adjust for areal extent of snow cover
        if (hru(j)%sno_mm < hru(j)%sno%covmx) then
          rto_sno = hru(j)%sno_mm / hru(j)%sno%covmx
          snocov = rto_sno / (rto_sno + Exp(hru(j)%snocov1 - hru(j)%snocov2 * rto_sno))
        else
          snocov = 1.
        endif
        snocov = max(0.0, min(1.0, snocov))
        
        !split rainfall into snow-covered and snow-free portions of the HRU
        rain_snow = rain_mm * snocov
        rain_bypass = rain_mm - rain_snow
 
        !calculate temperature-index snow melt
        if (w%tmax > hru(j)%sno%melttmp) then
          !! adjust degree-day melt factor for time of year
          smfac = (hru(j)%sno%meltmx + hru(j)%sno%meltmn) / 2. + Sin((time%day - 81) / 58.09) *     &
                        (hru(j)%sno%meltmx - hru(j)%sno%meltmn) / 2.        !! 365/2pi = 58.09
          melt_dd = smfac * (((snotmp + w%tmax)/2.) - hru(j)%sno%melttmp)
          melt_dd = max(0.0, melt_dd)
          melt_dd = melt_dd * snocov
        endif
          
        !! calculate rain-on-snow advection heat
        !! warm rainfall first compensates cold content. Only excess heat melts snow.
        if (rain_snow > 0.) then
          !aovid to use w%tave, because the following situation may occurs:
          !Tmin < 0, Tmax > 0, Tave < 0, rain_mm > 0.
          train_ros = max(0.0, 0.5 * (w%tave + w%tmax))
          rain_heat_melt_eq = rain_snow * cpw_lf * train_ros

          !cold content in mm melt equivalent.
          !this uses remaining solid snow SWE and negative snow temperature.
          cold_content = max(0.0, -snotmp) * hru(j)%sno_mm * cis_lf
          heat_excess = rain_heat_melt_eq - cold_content
          
          if (heat_excess > 0.) then
            melt_ros = heat_excess
            hru(j)%sno_tmp = 0.0
          else
            melt_ros = 0.0
            !optional: warm the pack toward 0 C if rain heat partly satisfies cold content.
            if (hru(j)%sno_mm > 1.e-6) then
              hru(j)%sno_tmp = - (cold_content - rain_heat_melt_eq) /            &
                               (hru(j)%sno_mm * cis_lf)
              hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
            endif
          endif
        endif
        
        !actual melt cannot exceed available solid snow water equivalent
        melt_pot = melt_dd + melt_ros
        melt_pot = max(0.0, melt_pot)
        melt_act = min(melt_pot, hru(j)%sno_mm)
        
        hru(j)%sno_mm = hru(j)%sno_mm - melt_act
        snomlt = melt_act

        !!route rainfall and meltwater through snowpack liquid storage.
        total_liquid_input = rain_snow + melt_act
        hru(j)%sno_liq = hru(j)%sno_liq + total_liquid_input
        
        !!simple refreezing of liquid water in a still-cold snowpack
        if (hru(j)%sno_liq > 0. .and. hru(j)%sno_tmp < 0. .and. hru(j)%sno_mm > 1.e-6) then
          cold_content = max(0.0, -hru(j)%sno_tmp) * hru(j)%sno_mm * cis_lf
          refreeze = min(hru(j)%sno_liq, cold_content)

          hru(j)%sno_liq = hru(j)%sno_liq - refreeze
          hru(j)%sno_mm = hru(j)%sno_mm + refreeze

          cold_content = max(0.0, cold_content - refreeze)
          if (hru(j)%sno_mm > 1.e-6) then
            hru(j)%sno_tmp = - cold_content / (hru(j)%sno_mm * cis_lf)
            hru(j)%sno_tmp = min(hru(j)%sno_tmp, 0.0)
          else
            hru(j)%sno_tmp = 0.0
          endif
        endif

        !liquid water holding capacity after melt/refreezing
        liq_capacity = hru(j)%sno_mm * liq_max_ratio
        if (hru(j)%sno_liq > liq_capacity) then
            !!snowpack collapses, releasing excess water
            liq_release = hru(j)%sno_liq - liq_capacity
            hru(j)%sno_liq = liq_capacity
        else
            !!all liquid water retained inside the snowpack sponge
            liq_release = 0.
        endif    
        
        !!liquid water released from the snowpack becomes effective precipitation.
        precip_eff = rain_bypass + liq_release
          
        if (precip_eff < 0.) precip_eff = 0.
        
      end if !!end of sno_mm > 0 block
      
      !!clear negligible snowpack and release remaining liquid water
      if (hru(j)%sno_mm < 1.e-6) then
          if (hru(j)%sno_liq > 0.) then
            precip_eff = precip_eff + hru(j)%sno_liq
          endif
          hru(j)%sno_mm = 0.
          hru(j)%sno_tmp = 0.
          hru(j)%sno_liq = 0.
      endif  
      
      if (time%step > 1 .and. snow_process) then
        w%ts(:) = precip_eff / time%step
      end if

      return
      end subroutine sq_snom