! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none

      integer :: time0, time1, clock_rate

      !ALP Parameters
      real(dp) :: axion_g10
      real(dp) :: axion_g
      real(dp) :: m_axion !in eV

      !Global arrays used for interpolation
      real(dp) :: mary(1001), kary(61), finit(1001*61), params(2)
      real(dp), pointer :: f1(:)
      real(dp), target :: farray(4*1001*61)
      integer :: MD


      ! these routines are called by the standard run_star check_model
      contains


      ! *************************************************************************************
      !                            MJD, FJH, RRV
      ! *************************************************************************************


      ! This function interpolates between the pre-processed grid of values stored in the work directory using MESA's interp2d function.
      ! It takes as input the value of log10(ma/kb*T) and log10(k_s/kb*T) specified in the other_neu subroutine below
      subroutine faxion_get(m_a, k_a, result, ierr, xary, yary, ND, f1)
          use interp_2d_lib_db, only: interp_rgbi3p_db, interp_mkbicub_db, interp_evbicub_db
          real(dp), intent(in) :: m_a, k_a, xary(1001), yary(61) ! Input: m_a (log10(ma/kb*T)), (k_a log10(k_s/kb*T))
                                                                 ! xary, yary - arrays of x and y values in interpolation range
          integer, intent(in) :: ND !=1 for new data, = 2 for old data
          real(dp), intent(out) :: result(6) ! output array - see interp_evbicub_db documentation for more information
          integer, intent(out) :: ierr
          real(dp), intent(inout), pointer :: f1(:) !Internal workspace for interp_evbicub_db
          integer :: i, j
          real(dp) :: xval(1), yval(1), zval(1), xmax(61), xmin(61), ymax(1001), ymin(1001)
          integer :: xsize, ysize, ilinx, iliny, ict(6)

          ict = [1,0,0,0,0,0] ! This array specifies what we would like interp_evbicub_db to calculate. See its documentation for more
                            ! information

          ! Boundary values, used for extrapolation if necessary
          ymin(:) = 0.0
          ymax(:) = zary(:, 61)
          xmin(:) = zary(1, :)
          xmax(:) = 0.0

          ! The point for which we want to know the z value. Has to be in an array for interp_evbicub_db to work.
          xval(1) = m_a
          yval(1) = k_a

          ! The size of each array
          xsize = SIZE(xary)
          ysize = SIZE(yary)

          ! Tells interp_mkbicub_db that xary and yary are linearly spaced.
          ilinx = 1
          iliny = 1

          if (ND == 1) then ! If ND==1, it implies that this is the first time this data is being used. interp_mkbicub_db needs to be run
                            ! to set up the internal workspace.
              call interp_mkbicub_db(xary, 1001, yary, 61, f1, 1001, 0, xmin, 0, xmax, 0, ymin, 0, ymax, ilinx, iliny, ierr)
          end if

          ! Call the interp_evbicub_db routine to determine the value of our function at m_a, k_a
          call interp_evbicub_db(m_a, k_a, xary, 1001, yary, 61, ilinx, iliny, f1, 1001, ict, result, ierr)
      end subroutine faxion_get



      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% other_neu => other_neu

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
      end subroutine extras_controls



      ! *************************************************************************************
      !                     MJD, FJH, RRV - based on code from 1210.1271
      ! *************************************************************************************

      subroutine other_neu(  &
                  id, k, T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, loss, sources, ierr)
          use neu_lib, only: neu_get
          use neu_def
          integer, intent(in) :: id ! id for star
          integer, intent(in) :: k ! cell number or 0 if not for a particular cell
          real(dp), intent(in) :: T ! temperature
          real(dp), intent(in) :: log10_T ! log10 of temperature
          real(dp), intent(in) :: Rho ! density
          real(dp), intent(in) :: log10_Rho ! log10 of density
          real(dp), intent(in) :: abar ! mean atomic weight
          real(dp), intent(in) :: zbar ! mean charge
          real(dp), intent(in) :: z2bar ! mean charge squared
          real(dp), intent(in) :: log10_Tlim
          logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
          real(dp), intent(out) :: loss(num_neu_rvs) ! total from all sources
          real(dp), intent(out) :: sources(num_neu_types, num_neu_rvs)
          integer, intent(out) :: ierr

          ! Parameters important for calculating energy-loss to ALP production
          real(dp) :: ye, axionz2ye, axioncsi, faxioncsi, d_faxioncsi_dT, &
             sprimakoff, d_sprimakoff_dT, d_sprimakoff_dRho, m_T, result, k_a, res(6), log10mT, log10ka
          integer :: i,j
          type (star_info), pointer :: s

          include 'formats'

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          ! Call the standard neu_get routine. We add the contribution of ALPs to this.
          call neu_get(  &
             T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, &
             loss, sources, ierr)
          if (ierr /= 0) return


          ! Step One: importing data and writing it to arrays
          ! Perform this once the first time the routine is called
          if (MD/=2) then
              ! Define arrays of log10(mu) and log10(xi) values as defined in the paper
              do i=1, 61
                  kary(i) = -4 + 0.1*(i-1)
              end do

              do i=1, 1001
                  mary(i) = -5 + 0.01*(i-1)
              end do

              ! Open the data file - this contains the array of (F(mu, xi)+G(mu)) values from the paper
              ! NB file path should be absolute if run on a cluster
              open(3, file = "ALP_Data.txt", &
                          status = 'old', access = 'sequential', form = 'formatted', action = 'read')
              read(3,*, end = 102) finit !Write the contents into the array finit
102	          close(3)

              ! Rewrite the contents of finit into the required form for the internal workspace of interp_mkbicub_db
              do i=1, 1001*61
                 farray(4*i-3) = finit(i)
                 farray(4*i-2) = 0.0
    		     farray(4*i-1) = 0.0
    		     farray(4*i) = 0.0
    	      end do
    	      f1 => farray

              ! The file bsm.txt contains the relevant ALP parameters (m_a and g10) - this was useful when running on a cluster
              ! This can be incorporated into an inlist using extras_rpar
              open(4, file = "bsm.txt", status = 'old', access = 'sequential', &
                             form = 'formatted', action = 'read')
    	      read(4,*) params
    	      close(4)
    		    axion_g10=params(1)
    		    m_axion=params(2)

    	      MD = 1
    	  end if
          ! ****** End Writing to Arrays *********


          ! Step Two: Calculate the value of epsilon_a

          ! Calculate log10(xi) as per our paper
          ! This section is modified from 1210.1271
          ye = s% ye(k)
          axionz2ye=z2bar/abar+ye
          axioncsi=  1.658d20*axionz2ye*Rho/(T*T*T)
          !.. coefficient = 4*pi*alpha/(4*(1 Kelvin)^3) * N_A * cm^(-3)
          !..    pi*(1/137.0.35)/(1/11604.5)^3 * (6.022*10^23) * (197.326*10^6*10^(-13))^3


          m_T = m_axion/(8.61733d-5*T)
          log10mT = log10_cr(m_T)


          k_a = axioncsi**(0.5)
          if (log10_cr(k_a) < 2) then
              log10ka = log10_cr(k_a)
          else
              log10ka = 1.99999 !Specify upper value for log10ka (fusion is dominant for large mass)
          end if

          ! Call routine if within the appropriate range - (don't want NANs in areas that aren't physically interesting)
          if (log10mT .LE. 5 .AND. log10mT .GE. -5) then
              call faxion_get(log10mT, log10ka, res, ierr, mary, kary, MD, f1)
              MD = 2 ! Arrays have been loaded and mkbicub run.
                     ! Change the value of global integer MD so this section will not be called again.
              faxioncsi = res(1)

              sprimakoff = 8.87736d-51*T**7*axion_g10**2*faxioncsi/(Rho)

              loss(ineu) = loss(ineu) + sprimakoff
              !loss(idneu_dT) = loss(idneu_dT) + d_sprimakoff_dT
              !loss(idneu_dRho) = loss(idneu_dRho) + d_sprimakoff_dRho
          end if
       end subroutine other_neu



      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         call system_clock(time0,clock_rate)
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         logical :: okay
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         okay = .true.
         write(*,1) 's% star_age', s% star_age
         if (s% star_age > 6e5) then
            write(*,1) 'star_age too large'
            okay = .false.
         end if
         write(*,1) 's% power_neutrinos', s% power_neutrinos
         if (s% power_neutrinos < 1e3) then
            write(*,1) 'power_neutrinos too small'
            okay = .false.
         end if
         write(*,1) 's% non_nuc_neu(s% nz)', s% non_nuc_neu(s% nz)
         if (s% non_nuc_neu(s% nz) < 2d3) then
            write(*,1) 'center non_nuc_neu too small'
            okay = .false.
         end if
         call system_clock(time1,clock_rate)
         dt = dble(time1 - time0) / clock_rate / 60
         write(*,'(/,a50,f12.2,99i10/)') 'runtime (minutes), retries, backups, steps', &
            dt, s% num_retries, s% num_backups, s% model_number
         ierr = 0
         if (okay) write(*,*) 'all tests okay'
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)
      end function extras_finish_step


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
