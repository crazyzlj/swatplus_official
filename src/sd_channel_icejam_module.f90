module sd_channel_icejam_module

    !! Centralized parameter and reach-scaling support for the conceptual
    !! SWAT+ river-ice / ice-jam module.
    !!
    !! Design goals:
    !!   1. Keep calibration parameters in one derived type instead of scattering
    !!      hard-coded constants inside process subroutines.
    !!   2. Replace channel-volume-normalized ice_ratio with reach-mean ice
    !!      thickness diagnostics:
    !!          sim_ice_thick = ice_volume / channel planform area
    !!          ice_maturity  = sim_ice_thick / ice_maturity_ref_thick
    !!      Seasonal phase transitions use absolute reach-mean ice-thickness
    !!      thresholds rather than a normalized maximum-ice fraction.
    !!   3. Replace stream-order hard coding of jam susceptibility with a continuous
    !!      normalized exponential transfer function of slope, sinuosity, and length.
    !!   4. Replace one-day relative q-rise triggering with a smoothed, noise-filtered
    !!      hydrologic trigger.

    implicit none

    private

    public :: icejam_param_type
    public :: icejam_reach_scale_type
    public :: icejam_default_params
    public :: icejam_validate_params
    public :: icejam_compute_reach_scale
    public :: icejam_sigmoid
    public :: icejam_clamp


    type :: icejam_param_type
        ! Ice thermodynamics and characteristic reach-mean ice-thickness scale.
        ! ice_maturity_ref_thick is a basin-scale reach-mean reference thickness
        ! used to normalize continuous process intensity. Seasonal phase
        ! transitions use absolute reach-mean ice-thickness thresholds.
        real :: ice_maturity_ref_thick = 0.30
        real :: ice_frz_tmp = -1.0
        real :: ice_melt_tmp = 0.0
        real :: ice_growth_coeff = 0.05
        real :: max_daily_ice_growth_thick = 0.015
        real :: ice_decay_coeff = 0.02
        real :: ice_freeze_inflow_frac = 0.05

        ! Freeze/thaw memory and seasonal state machine.
        real :: freeze_memory = 0.80
        real :: thaw_memory = 0.80
        real :: freezeup_freeze_dd = 3.0
        integer :: freezeup_min_days = 14
        integer :: freezeup_max_days = 45
        real :: freezeup_strong_index = 0.10
        integer :: deepwinter_min_days_before_breakup = 60
        integer :: deepwinter_max_days_before_breakup = 120
        integer :: warm_min_days_before_freezeup = 90
        integer :: breakup_min_days_before_warm = 21
        integer :: breakup_max_days_before_warm = 120
        real :: warm_season_weakening_index = 0.80
        real :: warm_storage_exit_ratio = 0.05
        real :: storage_cleanup_thaw_dd = 2.0
        real :: flush_thaw_dd = 35.0
        real :: breakup_onset_weakening_index = 0.25
        real :: jam_release_weakening_index = 0.30

        ! Absolute reach-mean ice-thickness thresholds for phase logic.
        real :: warm_ice_thick = 0.01
        real :: freezeup_ice_thick = 0.03
        real :: mobile_ice_thick = 0.03
        real :: retention_ice_thick = 0.08
        real :: deepwinter_ice_thick = 0.08
        real :: jam_material_ice_thick = 0.08
        real :: ros_ice_thick_min = 0.08

        ! Thaw and rain-on-ice forcing.
        real :: thaw_tmax_base = 4.0
        real :: thaw_tmax_base_ros = 1.5
        real :: thaw_tave_base = 0.0
        real :: ros_min_melt_mm = 0.5
        real :: warm_flush_tmp = 5.0

        ! Under-ice stable cover retention.
        real :: underice_cap_min = 0.05
        real :: underice_cap_open = 0.40
        real :: underice_cap_exp = 2.0
        real :: ice_cover_ret_frac_max = 0.80
        real :: ice_cover_ret_cap_coeff = 3.00
        real :: ice_cover_ret_stor_frac = 3.0
        real :: ice_cover_ret_q_damp = 0.30
        real :: freezeup_ret_mult = 0.50
        real :: deepwinter_ret_mult = 2.00
        real :: breakup_ret_mult = 0.70
        real :: deepwinter_ret_frac_min = 0.25
        real :: deepwinter_ret_frac_max = 0.90
        real :: deepwinter_cover_leak_frac = 0.005

        ! Flow reference. Raw flow is a water-supply term, not a BREAKUP trigger.
        real :: jam_ref_frac = 0.50

        ! BREAKUP jam-cycle memory and unified release parameters.
        integer :: release_auto_hold_days = 3
        real :: jam_leak_frac = 0.01
        real :: jam_release_frac_max = 0.15
        real :: jam_release_weak_exp = 2.0
        real :: release_ramp_day1_frac = 0.40
        real :: release_ramp_day2_frac = 0.70

        ! Unified jam-formation parameters.
        real :: jam_form_block_cap_coeff = 2.00
        real :: jam_form_stor_max_frac = 8.00
        real :: jam_form_block_frac_max = 0.60
        real :: jam_mobile_trigger_ratio = 0.02
        real :: cover_breakup_frac = 0.25
        real :: breakup_onset_cover_to_jam_frac = 0.50

        ! Mobile ice generation and routing.
        real :: mobile_q_min = 0.20
        real :: mobile_thaw_dd = 1.0
        real :: drift_mobilization_frac = 0.01
        real :: dynamic_mobilization_frac = 0.05
        real :: mobile_max_daily_frac = 0.20
        real :: mobile_breakup_drift_multiplier = 3.0
        real :: mobile_freezeup_drift_multiplier = 0.5
        real :: mobile_deepwinter_dynamic_weight = 0.3
        real :: ice_support_frac = 1.00

        ! Mobile ice capture in downstream reaches.
        real :: mobile_capture_base = 0.10
        real :: mobile_capture_susc_weight = 0.50
        real :: mobile_capture_ice_weight = 0.25
        real :: mobile_capture_depth_weight = 0.15
        real :: mobile_capture_min = 0.0
        real :: mobile_capture_max = 0.95
        real :: mobile_capture_capacity_mult = 2.0
        real :: warm_capture_max = 0.20
        real :: freezeup_capture_min = 0.70
        real :: deepwinter_capture_min = 0.80
        real :: breakup_capture_min = 0.40

        ! Continuous jam-susceptibility transfer function.
        real :: slope_ref = 0.01
        real :: slope_pow = 0.75
        real :: slope_min = 0.0
        real :: sinu_alpha = 1.50
        real :: length_ref = 1000.0
        real :: w_slope = 0.45
        real :: w_sinu = 0.30
        real :: w_inter = 0.15
        real :: w_len = 0.10
        real :: width_min = 0.10
        real :: depth_min = 0.05
        real :: jam_susc_min = 0.10
        real :: jam_susc_max = 0.95

        ! Convert susceptibility into process-specific modifiers.
        real :: jam_storage_base = 0.50
        real :: jam_storage_weight = 0.50
        real :: jam_block_base = 0.30
        real :: jam_block_weight = 0.70
        real :: jam_capture_base = 0.20
        real :: jam_capture_weight = 0.80
        real :: jam_mech_min = 0.20

    end type icejam_param_type

    type :: icejam_reach_scale_type
        real :: hyd_storage_scale = 1.0       !m3, channel hydraulic storage scale
        real :: ice_area = 1.0                !m2, planform water surface area
        real :: ice_cap_vol = 1.0             !m3, ice_maturity_ref_thick * ice_area
        real :: q_jam_ref_rate = 0.05         !m3/s
        real :: jam_susc = 0.5
        real :: jam_form_modifier = 0.5
        real :: jam_storage_modifier = 0.75
        real :: jam_block_modifier = 0.65
        real :: ice_capture_modifier = 0.60
        real :: mechanical_weakness_modifier = 1.0
        real :: mobile_order_mult = 1.0       !kept for backward-compatible mobile-ice scaling
        real :: slope_eff = 0.0
        real :: sinuosity_eff = 1.0
        real :: length_m = 1.0
        real :: susc_i_slope = 0.0
        real :: susc_i_sinu = 0.0
        real :: susc_i_inter = 0.0
        real :: susc_i_len = 0.0
        real :: susc_i_raw = 0.0
    end type icejam_reach_scale_type

contains

    real function icejam_clamp(x, xmin, xmax)
        real, intent(in) :: x, xmin, xmax
        icejam_clamp = max(xmin, min(xmax, x))
    end function icejam_clamp

    real function icejam_sigmoid(x)
        real, intent(in) :: x
        if (x > 50.) then
            icejam_sigmoid = 1.0
        else if (x < -50.) then
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

        if (p%ice_maturity_ref_thick <= 0.) stop "icejam parameter error: ice_maturity_ref_thick must be positive"
        if (p%slope_ref <= 0.) stop "icejam parameter error: slope_ref must be positive"
        if (p%slope_pow <= 0.) stop "icejam parameter error: slope_pow must be positive"
        if (p%length_ref <= 0.) stop "icejam parameter error: length_ref must be positive"
        if (p%w_slope < 0. .or. p%w_sinu < 0. .or. p%w_inter < 0. .or. p%w_len < 0.) &
                stop "icejam parameter error: jam susceptibility weights must be non-negative"
        if (p%w_slope + p%w_sinu + p%w_inter + p%w_len <= 0.) &
                stop "icejam parameter error: at least one jam susceptibility weight must be positive"
        if (p%jam_susc_min < 0. .or. p%jam_susc_max > 1. .or. p%jam_susc_min >= p%jam_susc_max) &
                stop "icejam parameter error: invalid jam susceptibility bounds"
        if (.not. (p%warm_ice_thick <= p%freezeup_ice_thick .and. &
                p%freezeup_ice_thick <= p%retention_ice_thick .and. &
                p%retention_ice_thick <= p%deepwinter_ice_thick .and. &
                p%jam_material_ice_thick <= p%ice_maturity_ref_thick)) then
            stop "icejam parameter error: absolute ice-thickness thresholds are inconsistent"
        endif
        if (p%jam_form_block_cap_coeff < 0. .or. p%jam_form_stor_max_frac < 0. .or. &
                p%jam_form_block_frac_max < 0. .or. p%jam_form_block_frac_max > 1.) &
                stop "icejam parameter error: invalid jam formation parameters"
        if (p%release_auto_hold_days < 1) &
                stop "icejam parameter error: release_auto_hold_days must be >= 1"
        if (p%ros_min_melt_mm < 0.) &
                stop "icejam parameter error: ros_min_melt_mm must be >= 0"
        if (p%breakup_onset_cover_to_jam_frac < 0. .or. p%breakup_onset_cover_to_jam_frac > 1.) &
                stop "icejam parameter error: breakup_onset_cover_to_jam_frac must be between 0 and 1"
        if (p%cover_breakup_frac < 0. .or. p%cover_breakup_frac > 1.) &
                stop "icejam parameter error: cover_breakup_frac must be between 0 and 1"
        if (p%jam_leak_frac < 0. .or. p%jam_leak_frac > 1.) &
                stop "icejam parameter error: jam_leak_frac must be between 0 and 1"
        if (p%jam_release_frac_max < 0. .or. p%jam_release_frac_max > 1.) &
                stop "icejam parameter error: jam_release_frac_max must be between 0 and 1"
        if (p%jam_release_weak_exp <= 0.) &
                stop "icejam parameter error: jam_release_weak_exp must be > 0"
        if (p%release_ramp_day1_frac < 0. .or. p%release_ramp_day1_frac > 1. .or. &
            p%release_ramp_day2_frac < 0. .or. p%release_ramp_day2_frac > 1.) &
                stop "icejam parameter error: release ramp fractions must be between 0 and 1"
        if (p%jam_leak_frac > p%jam_release_frac_max) &
                stop "icejam parameter error: jam_leak_frac must not exceed jam_release_frac_max"
        if (p%breakup_onset_weakening_index < 0. .or. p%breakup_onset_weakening_index > 1.) &
                stop "icejam parameter error: breakup_onset_weakening_index must be between 0 and 1"
        if (p%jam_release_weakening_index < 0. .or. p%jam_release_weakening_index > 1.) &
                stop "icejam parameter error: jam_release_weakening_index must be between 0 and 1"
        if (p%freezeup_strong_index < 0. .or. p%freezeup_strong_index > 1.) &
                stop "icejam parameter error: freezeup_strong_index must be between 0 and 1"
        if (p%freezeup_min_days < 0 .or. p%freezeup_max_days < p%freezeup_min_days) &
                stop "icejam parameter error: freezeup max days must be >= min days"
        if (p%deepwinter_min_days_before_breakup < 0 .or. &
            p%deepwinter_max_days_before_breakup < p%deepwinter_min_days_before_breakup) &
                stop "icejam parameter error: deepwinter max days must be >= min days"
        if (p%warm_min_days_before_freezeup < 0) &
                stop "icejam parameter error: warm_min_days_before_freezeup must be >= 0"
        if (p%breakup_min_days_before_warm < 0 .or. &
            p%breakup_max_days_before_warm < p%breakup_min_days_before_warm) &
                stop "icejam parameter error: breakup max days must be >= min days"
        if (p%warm_season_weakening_index < 0. .or. p%warm_season_weakening_index > 1.) &
                stop "icejam parameter error: warm_season_weakening_index must be between 0 and 1"
        if (p%warm_storage_exit_ratio < 0. .or. p%warm_storage_exit_ratio > 1.) &
                stop "icejam parameter error: warm_storage_exit_ratio must be between 0 and 1"
    end subroutine icejam_validate_params

    subroutine icejam_compute_reach_scale(p, ch_width, ch_length_km, ch_depth, ch_slope, &
            sinuosity, q_rate_low, q_rate_high, scale)
        type(icejam_param_type), intent(in) :: p
        real, intent(in) :: ch_width, ch_length_km, ch_depth, ch_slope
        real, intent(in) :: sinuosity, q_rate_low, q_rate_high
        type(icejam_reach_scale_type), intent(out) :: scale

        real :: width_eff, depth_eff, length_eff, slope_eff, sinu_eff
        real :: i_slope, i_sinu, i_len, i_inter, i_raw, i_den

        width_eff = max(p%width_min, ch_width)
        depth_eff = max(p%depth_min, ch_depth)
        length_eff = max(1.0, ch_length_km * 1000.0)
        slope_eff = max(0.0, ch_slope)
        sinu_eff = max(1.0, sinuosity)

        scale%hyd_storage_scale = max(1.0, length_eff * width_eff * depth_eff)
        scale%ice_area = max(1.0, length_eff * width_eff)
        scale%ice_cap_vol = max(1.0e-6, p%ice_maturity_ref_thick * scale%ice_area)

        ! Scheme B: normalized exponential weighted model.  slope=0 is valid
        ! and gives I_slope=1.0, so no divide-by-zero or log singularity occurs.
        i_slope = exp(- (slope_eff / max(p%slope_ref, 1.0e-9)) ** p%slope_pow)
        i_sinu = 1.0 - exp(- p%sinu_alpha * max(sinu_eff - 1.0, 0.0))
        i_len = 1.0 - exp(- length_eff / max(p%length_ref, 1.0e-6))
        i_inter = i_slope * i_sinu

        i_raw = p%w_slope * i_slope + p%w_sinu * i_sinu + &
                p%w_inter * i_inter + p%w_len * i_len
        i_den = max(p%w_slope + p%w_sinu + p%w_inter + p%w_len, 1.0e-6)

        scale%jam_susc = p%jam_susc_min + (p%jam_susc_max - p%jam_susc_min) * i_raw / i_den
        scale%jam_susc = icejam_clamp(scale%jam_susc, p%jam_susc_min, p%jam_susc_max)

        scale%jam_form_modifier = scale%jam_susc
        scale%jam_storage_modifier = icejam_clamp(p%jam_storage_base + &
                p%jam_storage_weight * scale%jam_susc, 0.0, 2.0)
        scale%jam_block_modifier = icejam_clamp(p%jam_block_base + &
                p%jam_block_weight * scale%jam_susc, 0.0, 2.0)
        scale%ice_capture_modifier = icejam_clamp(p%jam_capture_base + &
                p%jam_capture_weight * scale%jam_susc, 0.0, 2.0)
        scale%mechanical_weakness_modifier = max(p%jam_mech_min, scale%jam_susc)

        scale%q_jam_ref_rate = p%jam_ref_frac * sqrt(max(q_rate_low, 1.0e-6) * max(q_rate_high, 1.0e-6))
        scale%q_jam_ref_rate = max(scale%q_jam_ref_rate, 0.05)

        ! Keep the original intuition that larger/wider channels can mobilize
        ! and pass more broken ice, but compute it continuously from width.
        scale%mobile_order_mult = icejam_clamp(0.70 + 0.15 * log(max(width_eff, 1.0)), 0.70, 1.40)

        ! Diagnostics for debug output.
        scale%slope_eff = slope_eff
        scale%sinuosity_eff = sinu_eff
        scale%length_m = length_eff
        scale%susc_i_slope = i_slope
        scale%susc_i_sinu = i_sinu
        scale%susc_i_inter = i_inter
        scale%susc_i_len = i_len
        scale%susc_i_raw = i_raw / i_den
    end subroutine icejam_compute_reach_scale

end module sd_channel_icejam_module

