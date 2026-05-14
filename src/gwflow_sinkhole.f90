      subroutine gwflow_sinkhole !ljzhu 05/11/2026

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine determines the volume of groundwater that is added to the aquifer via sinkholes
!!    (recharge volumes are used in gwflow_simulate, in groundwater balance equations)

      use gwflow_module
      use hydrograph_module, only : ob,sp_ob,sp_ob1
      use hru_module, only : gwholeq

      implicit none

      integer :: i = 0                !           |counter
      integer :: j = 0                !           |counter for number of HRUs within an LSU
      integer :: k = 0                !           |counter
      integer :: n = 0                !           |counter
      integer :: s = 0                !           |solute counter
      integer :: hru_id = 0           !           |id of the HRU
      integer :: ob_num = 0           !           |object number of the HRU
      integer :: cell_id = 0          !           |id of the gwflow cell
      real :: hole_volume = 0.        !m3         |summation of sinkhole recharge from multiple HRUs
      real :: cell_hole_volume = 0.   !m3         |volume of sinkhole recharge to the cell
      real :: cell_weight = 0.

      !only proceed if conduit is active
      if (gw_sinkhole_flag == 1) then
      !use hru sinkhole recharge to calculate recharge (m3) cell values
      if (lsu_cells_link == 1) then !LSU-cell connection
        !loop through the landscape units
        !currently not implemented.
      else !proceed with HRU-cell connection
        !map hole volume from the HRUs to the grid cells
        do hru_id=1,sp_ob%hru
          ob_num = sp_ob1%hru + hru_id - 1
          if (gw_sinkhole_hruflag(hru_id) == 1 .and. gwholeq(hru_id) > 0.) then
              hole_volume = (gwholeq(hru_id)/1000.) * (ob(ob_num)%area_ha * 10000.) !m * m2 = m3
              if (gw_sinkhole_hruarea(hru_id) > 1.e-6) then
              do i=1,hru_num_cells(hru_id)
                cell_id = hru_cells(hru_id,i)
                if(gw_state(cell_id)%stat == 2) then !if boundary cell, give recharge to nearest active cell
                    cell_id = gw_bound_near(cell_id)
                endif
                if (gw_state(cell_id)%hole > 0) then
                  cell_weight = hru_cells_fract(hru_id,i) / gw_sinkhole_hruarea(hru_id)
                  cell_weight = max(0., cell_weight)  
                  cell_hole_volume = hole_volume * cell_weight
                  gw_hyd_ss(cell_id)%hole = gw_hyd_ss(cell_id)%hole + cell_hole_volume
                  gw_hyd_ss_yr(cell_id)%hole = gw_hyd_ss_yr(cell_id)%hole + cell_hole_volume !store for annual water
                  gw_hyd_ss_mo(cell_id)%hole = gw_hyd_ss_mo(cell_id)%hole + cell_hole_volume !store for monthly water

                  !if (gw_solute_flag == 1) then
                  !  do s=1,gw_nsolute !loop through the solutes
                  !    !currently not implemented.
                  !  enddo
                  !endif
                endif !gw_state(cell_id)%hole > 0
              enddo !loop hru_num_cells(hru_id)
              endif !gw_sinkhole_hruarea(hru_id) > 1.e-6
           endif ! gw_sinkhole_hruflag(hru_id) == 1 .and. gwholeq(hru_id) > 0.
        enddo !loop sp_ob%hru
        
        !write (9003,*) "after sinkhole, hole_volume:", gw_hyd_ss(3606)%hole

      endif !check for LSU-cell connection
      endif
      return
      end subroutine gwflow_sinkhole
