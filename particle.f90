!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.3 $ $Date: 2007/02/21 13:23:06 $
!

!=======================================================================
      subroutine write_xyfile(t1,t2,fname)
!=======================================================================
!  
!     This routine combines 2 SAC traces into 1 xy file with t1 and x
!     and t2 as y. Apart from LEVEN and ITYPE, the header for t1 is
!     written. ** THE TRACES ARE ASSUMED TO BE THE SAME LENGTH **
!
!-----------------------------------------------------------------------
      use f90sac  ! use the f90sac module
!-----------------------------------------------------------------------
      implicit none
      type (SACTrace) :: t1,t2 ! the traces

      integer :: lu
      
      integer, parameter :: int4 = selected_int_kind(9) ;
      integer, parameter :: real4 = selected_real_kind(6,37) ;
      
      
      character (len=4) :: tmp_char1,tmp_char2,tmp_char3,tmp_char4
      character (len=50):: fname
      integer :: i, istatus
      type (SACtrace) :: out
      logical UNITOK, UNITOP

      real(real4) :: sacrh(70) ! SAC floating point header
      integer(int4) :: sacih(40) ! SAC floating point header
      
!     * Open the file for writing *
!  ** Open the file for writing
      open (unit=f90sac_iounit, file=fname,form='unformatted', &
            access='direct',recl=f90sac_32bit_record_length, &
            status='unknown')

      lu = f90sac_iounit

      t1 % xmaximum = maxval(t1%trace(1:t1 % npts))
      t1 % xminimum = minval(t1%trace(1:t1 % npts))

      t1 % ymaximum = maxval(t2%trace(1:t1 % npts))
      t1 % yminimum = minval(t2%trace(1:t1 % npts))

!
!     * WRITE OUT THE SAC HEADER
!
      sacrh(001)  = t1%delta              
      sacrh(002)  = t1%depmin
      sacrh(003)  = t1%depmax
      sacrh(004)  = t1%scale
      sacrh(005)  = t1%odelta
      sacrh(006)  = t1%b
      sacrh(007)  = t1%e
      sacrh(008)  = t1%o
      sacrh(009)  = t1%a
      sacrh(010)  = t1%internal0
      sacrh(011)  = t1%t0
      sacrh(012)  = t1%t1
      sacrh(013)  = t1%t2
      sacrh(014)  = t1%t3
      sacrh(015)  = t1%t4
      sacrh(016)  = t1%t5
      sacrh(017)  = t1%t6
      sacrh(018)  = t1%t7
      sacrh(019)  = t1%t8
      sacrh(020)  = t1%t9
      sacrh(021)  = t1%f
      sacrh(022)  = t1%resp0
      sacrh(023)  = t1%resp1
      sacrh(024)  = t1%resp2
      sacrh(025)  = t1%resp3
      sacrh(026)  = t1%resp4
      sacrh(027)  = t1%resp5
      sacrh(028)  = t1%resp6
      sacrh(029)  = t1%resp7
      sacrh(030)  = t1%resp8
      sacrh(031)  = t1%resp9
      sacrh(032)  = t1%stla
      sacrh(033)  = t1%stlo
      sacrh(034)  = t1%stel
      sacrh(035)  = t1%stdp
      sacrh(036)  = t1%evla
      sacrh(037)  = t1%evlo
      sacrh(038)  = t1%evel
      sacrh(039)  = t1%evdp
      sacrh(040)  = t1%mag
      sacrh(041)  = t1%user0
      sacrh(042)  = t1%user1
      sacrh(043)  = t1%user2
      sacrh(044)  = t1%user3
      sacrh(045)  = t1%user4
      sacrh(046)  = t1%user5
      sacrh(047)  = t1%user6
      sacrh(048)  = t1%user7
      sacrh(049)  = t1%user8
      sacrh(050)  = t1%user9
      sacrh(051)  = t1%dist
      sacrh(052)  = t1%az
      sacrh(053)  = t1%baz
      sacrh(054)  = t1%gcarc
      sacrh(055)  = t1%internal1
      sacrh(056)  = t1%internal2
      sacrh(057)  = t1%depmen
      sacrh(058)  = t1%cmpaz
      sacrh(059)  = t1%cmpinc
      sacrh(060)  = t1%xminimum
      sacrh(061)  = t1%xmaximum
      sacrh(062)  = t1%yminimum
      sacrh(063)  = t1%ymaximum
      sacrh(064)  = t1%unused1
      sacrh(065)  = t1%unused2
      sacrh(066)  = t1%unused3
      sacrh(067)  = t1%unused4
      sacrh(068)  = t1%unused5
      sacrh(069)  = t1%unused6
      sacrh(070)  = t1%unused7

      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(sacrh,70)
      endif

      sacih(001) =t1%nzyear
      sacih(002) =t1%nzjday
      sacih(003) =t1%nzhour
      sacih(004) =t1%nzmin
      sacih(005) =t1%nzsec
      sacih(006) =t1%nzmsec
      sacih(007) =t1%nvhdr
      sacih(008) =t1%norid
      sacih(009) =t1%nevid
      sacih(010) =t1%npts
      sacih(011) =t1%internal3
      sacih(012) =t1%nwfid
      sacih(013) =t1%nxsize
      sacih(014) =t1%nysize
      sacih(015) =t1%unused8
      sacih(016) =4
      sacih(017) =t1%idep
      sacih(018) =t1%iztype
      sacih(019) =t1%unused9
      sacih(020) =t1%iinst
      sacih(021) =t1%istreg
      sacih(022) =t1%ievreg
      sacih(023) =t1%ievtyp
      sacih(024) =t1%iqual
      sacih(025) =t1%isynth
      sacih(026) =t1%imagtyp
      sacih(027) =t1%imagsrc
      sacih(028) =t1%unused10
      sacih(029) =t1%unused11
      sacih(030) =t1%unused12
      sacih(031) =t1%unused13
      sacih(032) =t1%unused14
      sacih(033) =t1%unused15
      sacih(034) =t1%unused16
      sacih(035) =t1%unused17
      sacih(036) =0!! BECAUSE A XY FILE NOW
      sacih(037) =t1%lpspol
      sacih(038) =t1%lovrok
      sacih(039) =t1%lcalda
      sacih(040) =t1%unused18

      if (f90sac_force_byteswap) then
         call f90sac_int32_byteswap(sacih,40)
      endif

      do i=1,40
         write(lu,rec=i+70) sacih(i)
      enddo

      
      write(lu,rec=111) t1%kstnm(1:4)
      write(lu,rec=112) t1%kstnm(5:8)
      write(lu,rec=113) t1%kevnm(1:4)
      write(lu,rec=114) t1%kevnm(5:8)
      write(lu,rec=115) t1%kevnm(9:12)
      write(lu,rec=116) t1%kevnm(13:16)
      write(lu,rec=117) t1%khole(1:4)
      write(lu,rec=118) t1%khole(5:8)
      write(lu,rec=119) t1%ko(1:4)
      write(lu,rec=120) t1%ko(5:8)
      write(lu,rec=121) t1%ka(1:4)
      write(lu,rec=122) t1%ka(5:8)
      write(lu,rec=123) t1%kt0(1:4)
      write(lu,rec=124) t1%kt0(5:8)
      write(lu,rec=125) t1%kt1(1:4)
      write(lu,rec=126) t1%kt1(5:8)
      write(lu,rec=127) t1%kt2(1:4)
      write(lu,rec=128) t1%kt2(5:8)
      write(lu,rec=129) t1%kt3(1:4)
      write(lu,rec=130) t1%kt3(5:8)
      write(lu,rec=131) t1%kt4(1:4)
      write(lu,rec=132) t1%kt4(5:8)
      write(lu,rec=133) t1%kt5(1:4)
      write(lu,rec=134) t1%kt5(5:8)
      write(lu,rec=135) t1%kt6(1:4)
      write(lu,rec=136) t1%kt6(5:8)
      write(lu,rec=137) t1%kt7(1:4)
      write(lu,rec=138) t1%kt7(5:8)
      write(lu,rec=139) t1%kt8(1:4)
      write(lu,rec=140) t1%kt8(5:8)
      write(lu,rec=141) t1%kt9(1:4)
      write(lu,rec=142) t1%kt9(5:8)
      write(lu,rec=143) t1%kf(1:4)
      write(lu,rec=144) t1%kf(5:8)
      write(lu,rec=145) t1%kuser0(1:4)
      write(lu,rec=146) t1%kuser0(5:8)
      write(lu,rec=147) t1%kuser1(1:4)
      write(lu,rec=148) t1%kuser1(5:8)
      write(lu,rec=149) t1%kuser2(1:4)
      write(lu,rec=150) t1%kuser2(5:8)
      write(lu,rec=151) t1%kcmpnm(1:4)
      write(lu,rec=152) t1%kcmpnm(5:8)
      write(lu,rec=153) t1%knetwk(1:4)
      write(lu,rec=154) t1%knetwk(5:8)
      write(lu,rec=155) t1%kdatrd(1:4)
      write(lu,rec=156) t1%kdatrd(5:8)
      write(lu,rec=157) t1%kinst(1:4)
      write(lu,rec=158) t1%kinst(5:8)


!  ** if required, byteswap the first trace
      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(t1%trace,t1%npts)
      endif

!  ** Output the first trace   
       do i=1,t1 % npts
            write(lu,rec=158+i) t1%trace(i)
      enddo

!  ** if required, byteswap the first trace
      if (f90sac_force_byteswap) then
         call f90sac_real32_byteswap(t2%trace,t2%npts)
      endif

!  ** Output the first trace   
       do i=1,t2 % npts
            write(lu,rec=158+i) t2%trace(i)
      enddo


      close(lu)
      

      return           
      end subroutine write_xyfile
!=======================================================================
