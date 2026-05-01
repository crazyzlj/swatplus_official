subroutine sd_channel_icejam(j)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Conceptual daily minor/major ice-jam storage-release module.
!!
!!    Concept:
!!      1. Ice cover represents blockage potential, not the main water source.
!!      2. Minor jams occur more frequently but only weakly obstruct flow.
!!      3. Major jams occur under stronger breakup forcing and can store/release
!!         a larger amount of incoming snowmelt/rain-on-snow flood water.
!!      4. Water balance is preserved through:
!!            blocked  : ht1%flo -> ice_jam_stor
!!            released : ice_jam_stor -> ht1%flo
!!

      use basin_module
      use time_module
      use hydrograph_module
      use sd_channel_module
      use climate_module

      implicit none

      integer, intent(in) :: j

      !integer :: iwst = 0
      integer :: ord = 1

      real :: ch_vol_cap = 0.          !m3      |approximate bankfull channel volume
      real :: ice_cover_max = 0.       !m3      |maximum ice-cover condition state
      real :: ice_cover_min = 0.       !m3      |minimum ice cover needed for jam formation
      real :: ice_cover_major_min = 0. !m3      |minimum ice cover needed for major jam formation

      real :: q_in_rate_raw = 0.       !m3/s    |raw inflow rate before ice-jam adjustment
      real :: q_bank_rate = 0.         !m3/s    |bankfull flow rate
      real :: q_ratio = 0.             !none    |q_in_rate_raw / q_bank_rate
      real :: q_rise_rate = 0.         !none    |relative daily rise in raw inflow

      real :: tw_ice = 0.              !deg C   |water temperature proxy for ice processes
      real :: t_air = 0.               !deg C   |daily mean air temperature
      real :: tmax = 0.                !deg C   |daily maximum air temperature
      real :: t_ice_growth = 0.        !deg C   |temperature driver for ice growth
      real :: t_ice_decay = 0.         !deg C   |temperature driver for ice decay / breakup

      real :: ice_growth = 0.          !m3/day  |growth of ice-cover condition state
      real :: ice_decay = 0.           !m3/day  |thermal decay of ice-cover condition state

      real :: jam_susc = 0.            !none    |channel ice-jam susceptibility
      real :: block_frac = 0.          !none    |fraction of incoming water blocked by jam
      real :: stor_max_frac = 0.       !none    |maximum jam storage as fraction of channel volume
      real :: jam_stor_max = 0.        !m3      |maximum water storage behind ice jam
      real :: release_frac = 0.        !none    |fast release fraction during breakup
      real :: leak_frac = 0.           !none    |slow leakage fraction while jam persists

      real :: blocked = 0.             !m3/day  |water blocked into ice-jam storage
      real :: released = 0.            !m3/day  |water released from ice-jam storage

      real :: raw_flo = 0.             !m3/day  |raw inflow volume before ice-jam adjustment
      real :: adj_ratio = 1.           !none    |ratio for subdaily tsin adjustment

      logical :: minor_jam_possible = .false.
      logical :: major_jam_possible = .false.
      logical :: breakup_condition = .false.

      !! ------------------------------------------------------------------
      !! Parameters for first diagnostic implementation.
      !! Later these can be moved to basin/channel parameter files.
      !! ------------------------------------------------------------------
      !ice-cover state related parameters
      real :: ice_frz_tmp = -1.0       !deg C   |ice-cover growth threshold
      real :: ice_melt_tmp = 0.0       !deg C   |breakup temperature threshold
      real :: ice_growth_coeff = 0.03  !1/degC/day |ice-cover growth coefficient
      real :: ice_decay_coeff = 0.15   !1/degC/day |ice-cover decay coefficient
      real :: ice_max_frac = 0.20      !none       |

      !minor ice jam parameters: frequent, weak obstruction
      real :: minor_q_ratio = 0.10
      real :: minor_q_rise = 0.15
      real :: minor_block_frac = 0.07
      real :: minor_stor_max_frac = 0.40
      real :: minor_release_frac = 0.30
      real :: minor_leak_frac = 0.04

      !major jam parameters: less frequent, stronger flood-wave storage/release
      real :: major_q_ratio = 0.20
      real :: major_q_rise = 0.30
      real :: major_block_frac = 0.25
      real :: major_stor_max_frac = 2.50
      real :: major_release_frac = 0.75
      real :: major_leak_frac = 0.02

      !maximum blocking limits
      real :: minor_block_frac_max = 0.30
      real :: major_block_frac_max = 0.80

      
      ich = j
      iwst = ob(icmd)%wst
      t_air = wst(iwst)%weat%tave
      tmax = wst(iwst)%weat%tmax

      sd_ch(ich)%icejam_block = 0.
      sd_ch(ich)%icejam_release = 0.
      sd_ch(ich)%icejam_qraw = 0.
      sd_ch(ich)%icejam_qadj = 0.
      sd_ch(ich)%icejam_qratio = 0.
      sd_ch(ich)%icejam_qrise = 0.
      sd_ch(ich)%icejam_susc = 0.

      !! ------------------------------------------------------------
      !! Channel-order-based ice-jam susceptibility
      !! ------------------------------------------------------------
      ord = sd_ch(ich)%order

      !if (max_order > 1) then
      !  jam_susc = jam_susc_min + (1.0 - jam_susc_min) * &
      !             real(ord - 1) / real(max_order - 1)
      !else
      !  jam_susc = 1.0
      !endif
      select case (ord)
      case (1)
        jam_susc = 0.10
      case (2)
        jam_susc = 0.25
      case (3)
        jam_susc = 0.45
      case (4)
        jam_susc = 0.70
      case default
        jam_susc = 1.00
      end select
      if (sd_ch(ich)%chs > 1.e-9 .and. sd_ch(ich)%chs < 0.001) then
        jam_susc = min(1.0, jam_susc + 0.15)
      endif
      sd_ch(ich)%icejam_susc = jam_susc


      !raw inflow before ice-jam modification.
      q_in_rate_raw = max(0., ht1%flo) / 86400.
      sd_ch(ich)%icejam_qraw = q_in_rate_raw

      !approximate channel capacity and bankfull flow
      ch_vol_cap = sd_ch(ich)%chl * 1000. * sd_ch(ich)%chw * sd_ch(ich)%chd
      ch_vol_cap = max(ch_vol_cap, 1.)

      ice_cover_max = ice_max_frac * ch_vol_cap
      ice_cover_min = 0.05 * ice_cover_max
      ice_cover_major_min = 0.20 * ice_cover_max

      !bankfull flow rate
      if (sd_ch(ich)%bankfull_flo > 1.e-6) then
        q_bank_rate = sd_ch(ich)%bankfull_flo * ch_rcurv(ich)%elev(2)%flo_rate
      else
        q_bank_rate = ch_rcurv(ich)%elev(2)%flo_rate
      endif
      q_bank_rate = max(q_bank_rate, 1.e-6)
      q_ratio = q_in_rate_raw / q_bank_rate
      sd_ch(ich)%icejam_qratio = q_ratio

      !relative flow rise based on raw inflow
      q_rise_rate = 0.
      if (sd_ch(ich)%q_prev > 1.e-6) then
        q_rise_rate = (q_in_rate_raw - sd_ch(ich)%q_prev) / sd_ch(ich)%q_prev
      endif
      q_rise_rate = max(0., q_rise_rate)
      sd_ch(ich)%icejam_qrise = q_rise_rate

      !water temperature for ice processes
      tw_ice = sd_ch(ich)%tmp_prx

      if (tw_ice < -20. .or. tw_ice > 40.) then
        tw_ice = t_air
      endif

      !ice growth should respond to cold water/air conditions.
      t_ice_growth = min(tw_ice, t_air)

      !ice decay/breakup should not be driven by daily Tmax alone.
      !use a combined water-air thermal proxy.
      t_ice_decay = 0.5 * tw_ice + 0.5 * t_air

      !! ------------------------------------------------------------------
      !! 1. Develop ice-cover condition during cold periods.
      !! Do not directly freeze daily inflow ht1%flo.
      !! ------------------------------------------------------------------
      if (t_ice_growth < ice_frz_tmp) then
        ice_growth = ice_growth_coeff * (ice_frz_tmp - t_ice_growth) * ice_cover_max
        ice_growth = max(0., ice_growth)
        ice_growth = min(ice_growth, ice_cover_max - sd_ch(ich)%ice)
        sd_ch(ich)%ice = sd_ch(ich)%ice + ice_growth
      endif

      !! ------------------------------------------------------------------
      !! 2. Thermal deterioration of ice cover.
      !! Ice decay reduces blockage potential. 
      !! ------------------------------------------------------------------
      if (t_ice_decay > ice_melt_tmp .and. sd_ch(ich)%ice > 0.) then
        ice_decay = ice_decay_coeff * (t_ice_decay - ice_melt_tmp) * sd_ch(ich)%ice
        ice_decay = max(0., ice_decay)
        ice_decay = min(ice_decay, sd_ch(ich)%ice)
        sd_ch(ich)%ice = sd_ch(ich)%ice - ice_decay
      endif

      !! ------------------------------------------------------------------
      !! 3. Release existing jam-stored water during breakup.
      !! This is done before forming a new jam to avoid same-day block-release
      !! cancellation.
      !! ------------------------------------------------------------------
      breakup_condition = .false.

      if (sd_ch(ich)%ice_jam_stor > 0.) then
        !major break condition: warming plus moderate/high flow or rapid rise
        if (t_ice_decay > 0.) then
          if (q_ratio > major_q_ratio .or. q_rise_rate > major_q_rise) then
            breakup_condition = .true.
          endif
        end if
        select case (sd_ch(ich)%ice_jam_flag)
        case (2)
          release_frac = major_release_frac
          leak_frac = major_leak_frac
        case (1)
          release_frac = minor_release_frac
          leak_frac = minor_leak_frac
        case default
          release_frac = minor_release_frac
          leak_frac = minor_leak_frac
        end select

        if (breakup_condition) then
          released = release_frac * sd_ch(ich)%ice_jam_stor
        else
          released = leak_frac * sd_ch(ich)%ice_jam_stor
        endif

        released = max(0., released)
        released = min(released, sd_ch(ich)%ice_jam_stor)

        sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor - released
        ht1%flo = ht1%flo + released

        sd_ch(ich)%icejam_release = released

        if (sd_ch(ich)%ice_jam_stor < 1.e-6) then
          sd_ch(ich)%ice_jam_stor = 0.
          sd_ch(ich)%ice_jam_flag = 0.
        endif

      endif

      !! ------------------------------------------------------------------
      !! 4. Determine whether a new ice jam can form.
      !! Ice jam formation requires existing ice cover, warming condition,
      !! and either moderate/high flow or rapid daily flow rise.
      !! ------------------------------------------------------------------
      major_jam_possible = .false.
      minor_jam_possible = .false.

      if (jam_susc >= 0.50) then
        if (sd_ch(ich)%ice > ice_cover_major_min) then
          if (t_ice_decay > 0.) then
            if (q_ratio > major_q_ratio .or. q_rise_rate > major_q_rise) then
              major_jam_possible = .true.
            endif
          endif
        endif
      endif

      if (.not. major_jam_possible) then
        if (sd_ch(ich)%ice > ice_cover_min) then
          if (t_ice_decay > 0.) then
            if (q_ratio > minor_q_ratio .or. q_rise_rate > minor_q_rise) then
              minor_jam_possible = .true.
            endif
          endif
        endif
      endif

      !! ------------------------------------------------------------------
      !! 5.  Set jam parameters for today's new blocking process
      !! ------------------------------------------------------------------
      block_frac = 0.
      stor_max_frac = 0.
      jam_stor_max = 0.

      if (major_jam_possible) then

        block_frac = major_block_frac * jam_susc
        block_frac = max(0., min(major_block_frac_max, block_frac))
        stor_max_frac = major_stor_max_frac * jam_susc
        sd_ch(ich)%ice_jam_flag = 2.

      else if (minor_jam_possible) then

        block_frac = minor_block_frac * jam_susc
        block_frac = max(0., min(minor_block_frac_max, block_frac))
        stor_max_frac = minor_stor_max_frac * jam_susc
        if (sd_ch(ich)%ice_jam_flag < 1.) sd_ch(ich)%ice_jam_flag = 1.

      endif

      jam_stor_max = stor_max_frac * ch_vol_cap
      jam_stor_max = max(0., jam_stor_max)
      
      !! ------------------------------------------------------------------
      !! 6. Block part of incoming flood water behind the ice jam.
      !! This is the main ice-jam water storage process.
      !! ------------------------------------------------------------------
      if (block_frac > 1.e-6 .and. jam_stor_max > 1.e-6) then

        blocked = block_frac * ht1%flo
        blocked = min(blocked, jam_stor_max - sd_ch(ich)%ice_jam_stor)
        blocked = max(0., blocked)

        ht1%flo = ht1%flo - blocked
        sd_ch(ich)%ice_jam_stor = sd_ch(ich)%ice_jam_stor + blocked

        sd_ch(ich)%icejam_block = blocked

      endif

      !! ------------------------------------------------------------------
      !! 7. Final synchronization with the hydrograph used by ch_rtmusk.
      !! ch_rtmusk uses ob(icmd)%tsin(irtstep), so this must be updated
      !! after all release/blocking operations.
      !! ------------------------------------------------------------------
      ht1%flo = max(0., ht1%flo)

      if (time%step == 1) then
        ob(icmd)%tsin(1) = ht1%flo
      else
        raw_flo = max(1.e-6, q_in_rate_raw * 86400.)
        adj_ratio = ht1%flo / raw_flo
        adj_ratio = max(0., adj_ratio)
        ob(icmd)%tsin(:) = ob(icmd)%tsin(:) * adj_ratio
      endif

      sd_ch(ich)%icejam_qadj = ht1%flo / 86400.

      !! Important:
      !! q_prev must store raw inflow before ice-jam adjustment.
      sd_ch(ich)%q_prev = q_in_rate_raw

      return 

end subroutine sd_channel_icejam
