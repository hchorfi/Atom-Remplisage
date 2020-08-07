!***********************************************************************
!
!    Copyright (c) 2012, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Nicolas Schunck, schunck1@llnl.gov
!
!    LLNL-CODE-573953 All rights reserved.
!
!    Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                    N. Michel, J. Sarich, S. Wild
!    Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************
!==================================================================================================================================
!#START MAINPROGRAM
!==================================================================================================================================
  !================================================================================================================================
  ! PARALLEL RUN ex script uses these two lines. Do not modify them!
  Program hfbthoprog; Use HFBTHO; Use HFBTHO_utilities; Call Main_Program; End Program hfbthoprog
  !================================================================================================================================
  !
  Subroutine Main_Program
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: i,it,icount,l,noForce
    !
    !-------------------------------------------------------------------
    ! Read input and namelist data for the requested nucleus
    !-------------------------------------------------------------------
    Call initialize_HFBTHO_NAMELIST
    !
    Call read_HFBTHO_NAMELIST
    !memo: Namelist /HFBTHO_GENERAL/ number_of_shells,oscillator_length,&
    !                                proton_number,neutron_number,type_of_calculation
    !      Namelist /HFBTHO_ITERATIONS/ number_iterations, accuracy
    !      Namelist /HFBTHO_FUNCTIONAL/ functional, add_initial_pairing, type_of_coulomb
    !      Namelist /HFBTHO_CONSTRAINTS/ lambda_values, lambda_active, expectation_values
    !      Namelist /HFBTHO_BLOCKING/ proton_blocking(1:5), neutron_blocking(1:5)
    !      Namelist /HFBTHO_BLOCKING/ switch_to_THO, projection_is_on, gauge_points, delta_N, delta_P
    !      Namelist /HFBTHO_TEMPERATURE/ set_temperature, temperature
    !      Namelist /HFBTHO_DEBUG/ number_Gauss, number_Laguerre, number_Legendre, &
    !                              force_parity, print_time
    !
    n00_INI                   = number_of_shells    ! number of shells
    b0_INI                    = oscillator_length   ! oscillator length
    q_INI                     = basis_deformation   ! deformation beta_2 of the basis
    npr_INI(1)                = neutron_number      ! N
    npr_INI(2)                = proton_number       ! Z
    kindhfb_INI               = type_of_calculation ! 1: HFB, -1: HFB+LN
    !
    MAX_ITER_INI              = number_iterations   ! max number of iterations
    epsi_INI                  = accuracy            ! convergence of iterations
    inin_INI                  = restart_file        ! restart from file
    !
    skyrme_INI                = TRIM(functional)    ! functional
    Add_Pairing_INI           = add_initial_pairing ! add pairing starting from file
    icou_INI                  = type_of_coulomb     ! coul: no-(0), dir.only-(1), plus exchange-(2)
    !
    set_pairing               = user_pairing        ! pairing is defined by user if .True.
    V0n_INI                   = vpair_n             ! pairing strength for neutrons
    V0p_INI                   = vpair_p             ! pairing strength for protons
    pwi_INI                   = pairing_cutoff      ! pairing q.p. cutoff
    cpv1_INI                  = pairing_feature     ! Type of pairing: volume, surface, mixed
    !
    nkblocase(1,:)            = neutron_blocking    ! config. of neutron blocked state
    nkblocase(2,:)            = proton_blocking     ! config. of proton blocked state
    !
    iLST_INI                  = switch_to_THO       ! 0:HO, -1:HO->THO, 1:THO
    keypj_INI                 = gauge_points        ! PNP: number of gauge points
    iproj_INI                 = projection_is_on    ! projecting on different nucleus
    npr1pj_INI                = delta_N             ! its neutron number
    npr2pj_INI                = delta_Z             ! its proton number
    !
    switch_on_temperature     = set_temperature     ! switches on temperature mode
    temper                    = temperature         ! value of the temperature
    !
    ngh_INI                   = number_Gauss        ! number of Gauss-Hermite points for z-direction
    ngl_INI                   = number_Laguerre     ! number of Gauss-Laguerre points for rho-direction
    nleg_INI                  = number_Legendre     ! number of Gauss-Legendre points for Coulomb
    basis_HFODD_INI           = compatibility_HFODD ! flag to enforce same basis as HFODD
    nstate_INI                = number_states       ! total number of states in basis
    Parity_INI                = force_parity        ! reflection symmetry
    IDEBUG_INI                = print_time          ! debug
    DO_FITT_INI               = .False.             ! calculates quantities for reg.optimization
    Print_HFBTHO_Namelist_INI = .False.             ! Print Namelist
    !
    If(SUM(lambda_active).Gt.0) Then
       icount=0
       Do l=1,lambdaMax,2
          If(lambda_active(l).Gt.0) icount=icount+1
       End Do
       If(icount.Gt.0) Parity_INI=.False.
    End If
    !
    Call read_UNEDF_NAMELIST(skyrme_INI,noForce)
    ! If functional is used, projection automaticaly switched off
    If(noForce.Eq.0) iproj_INI=0
    !---------------------------------------------------------------------------
    ! GROUND STATE BLOCKING WALKER: blocking candidates are predefined
    ! by the parent nucleus and we block them one by one
    !---------------------------------------------------------------------------
    iblocase=0; bloqpdif=zero !blomax=0; blomax will be charged from the previous solution
    Do it=1,2
       If(nkblocase(it,1).Ne.0.And.nkblocase(it,2).Eq.0) Then
          If(it.Eq.1) Then
             iblocase(1)=iblocase(1)+1
             If(iblocase(1).Gt.blomax(1)) iblocase(1)=1
          Else
             If(iblocase(1).Le.1) iblocase(2)=iblocase(2)+1
          End If
          nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
          nkblo_INI(it,2)=0
       Else
          ! case of external blocking
          nkblo_INI(it,:)=nkblocase(it,:)
       Endif
    End Do
    !---------------------------------------------------------------------------
    ! MANUAL BLOCKING: manualBlocking=1 in the module (ocasionaly used)
    ! One types which level to be blocked referencing the parent nucleus
    !---------------------------------------------------------------------------
    If(manualBlocking.Ne.0) Then
       Write(*,'(a,5(1pg12.4))') 'Please print the number of the neutron level to block, num='
       Read(*,*) nkblo_INI(1,1)
       Write(*,'(a,5(1pg12.4))') 'Neutron blocked level num=',nkblo_INI(1,1)
       nkblo_INI(1,2)=0
       Write(*,'(a,5(1pg12.4))') 'Please print the number of the proton level to block, num='
       Read(*,*) nkblo_INI(2,1)
       Write(*,'(a,5(1pg12.4))') 'Proton blocked level num=',nkblo_INI(2,1)
       nkblo_INI(2,2)=0
    End If
    !--------------------------------------------------------------------
    ! Calculations for 'FITS' functional (Modifies some values if needed)
    !--------------------------------------------------------------------
    If (Trim(skyrme_INI).Eq.'FITS') Then
       !
       DMEORDER=-1; DMELDA=0; use_cm_cor=.False.
       !
       HBZERO =  20.7355300000000D0; E2CHARG =   1.4399784000000d0
       CRHO(0)=-731.2227858295098d0; CDRHO(0)= 855.6900515849785d0; CTAU(0) = -0.5439888609059821d0
       CRHO(1)= 263.7103055246761d0; CDRHO(1)=-176.8641956040411d0; CTAU(1) =-33.3618818665213400d0
       CRDR(0)= -43.2900897553154d0; CJ(0)   =   0.0000000000000d0; CRDJ(0) =-75.2608700482894100d0
       CRDR(1)=-164.1379857135440d0; CJ(1)   =   0.0000000000000d0; CRDJ(1) =-22.6528199648713000d0
       CPV0(0)=-186.1922465962490d0; CPV1(0) =   0.5000000000000d0; SIGMA   =  0.2987839827782357d0
       CPV0(1)=-206.7464168983860d0; CPV1(1) =   0.5000000000000d0; CEXPAR  =  0.6391295237623640d0
       !
       RHO_NM=0.15732716296394680d0;
       E_NM=-15.80000048487058000d0; K_NM    =225.94339488338960d0; SMASS_NM=  0.9958725808228520d0
       ASS_NM=28.3483385569865900d0; LASS_NM =40.001962979089330d0; VMASS_NM=  1.2489999532699580d0
       !
    End If
    !--------------------------------------------------------------------
    ! Run the solver in all cases EVEN/ODDS, FITS/NO-FITS
    !--------------------------------------------------------------------
    Call HFBTHO_SOLVER
    !--------------------------------------------------------------------
    ! Display error messages in case of problems
    !--------------------------------------------------------------------
    If (ierror_flag.Ne.0) Then
        Write(*,*)
        Write(*,'(a)') ' ERRORS IN HFBTHO_SOLVER'
        Do i=1,ierror_flag
           Write(*,'(a,i2,2x,a)') ' error_flag=',i,ierror_info(i)
        End Do
        Write(*,*)
        Else
        Write(*,*)
        Write(*,'(a)') ' HFBTHO_SOLVER ended without errors'
        Write(*,*)
    End If
    !
    If (lout.Lt.lfile) Close(lfile) ! close the output
    !
  End Subroutine Main_Program
!==================================================================================================================================
!#END MAINPROGRAM
!==================================================================================================================================
