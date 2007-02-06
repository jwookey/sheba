!===============================================================================
!  S H E B A - Shear-wave Birefringence Analysis
!===============================================================================
!  This software is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!===============================================================================
!
!  James Wookey, School of Earth Sciences, University of Bristol
!  CVS: $Revision: 1.2 $ $Date: 2007/02/06 14:03:49 $
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
      character (len=4) :: tmp_char1,tmp_char2,tmp_char3,tmp_char4
      character (len=50):: fname
      integer :: i, istatus
      type (SACtrace) :: out
      logical UNITOK, UNITOP
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
!     * READ IN THE SAC HEADER
!
      write(lu,rec=001) t1%delta              
      write(lu,rec=002) t1%depmin
      write(lu,rec=003) t1%depmax
      write(lu,rec=004) t1%scale
      write(lu,rec=005) t1%odelta
      write(lu,rec=006) t1%b
      write(lu,rec=007) t1%e
      write(lu,rec=008) t1%o
      write(lu,rec=009) t1%a
      write(lu,rec=010) t1%internal0
      write(lu,rec=011) t1%t0
      write(lu,rec=012) t1%t1
      write(lu,rec=013) t1%t2
      write(lu,rec=014) t1%t3
      write(lu,rec=015) t1%t4
      write(lu,rec=016) t1%t5
      write(lu,rec=017) t1%t6
      write(lu,rec=018) t1%t7
      write(lu,rec=019) t1%t8
      write(lu,rec=020) t1%t9
      write(lu,rec=021) t1%f
      write(lu,rec=022) t1%resp0
      write(lu,rec=023) t1%resp1
      write(lu,rec=024) t1%resp2
      write(lu,rec=025) t1%resp3
      write(lu,rec=026) t1%resp4
      write(lu,rec=027) t1%resp5
      write(lu,rec=028) t1%resp6
      write(lu,rec=029) t1%resp7
      write(lu,rec=030) t1%resp8
      write(lu,rec=031) t1%resp9
      write(lu,rec=032) t1%stla
      write(lu,rec=033) t1%stlo
      write(lu,rec=034) t1%stel
      write(lu,rec=035) t1%stdp
      write(lu,rec=036) t1%evla
      write(lu,rec=037) t1%evlo
      write(lu,rec=038) t1%evel
      write(lu,rec=039) t1%evdp
      write(lu,rec=040) t1%mag
      write(lu,rec=041) t1%user0
      write(lu,rec=042) t1%user1
      write(lu,rec=043) t1%user2
      write(lu,rec=044) t1%user3
      write(lu,rec=045) t1%user4
      write(lu,rec=046) t1%user5
      write(lu,rec=047) t1%user6
      write(lu,rec=048) t1%user7
      write(lu,rec=049) t1%user8
      write(lu,rec=050) t1%user9
      write(lu,rec=051) t1%dist
      write(lu,rec=052) t1%az
      write(lu,rec=053) t1%baz
      write(lu,rec=054) t1%gcarc
      write(lu,rec=055) t1%internal1
      write(lu,rec=056) t1%internal2
      write(lu,rec=057) t1%depmen
      write(lu,rec=058) t1%cmpaz
      write(lu,rec=059) t1%cmpinc
      write(lu,rec=060) t1%xminimum
      write(lu,rec=061) t1%xmaximum
      write(lu,rec=062) t1%yminimum
      write(lu,rec=063) t1%ymaximum
      write(lu,rec=064) t1%unused1
      write(lu,rec=065) t1%unused2
      write(lu,rec=066) t1%unused3
      write(lu,rec=067) t1%unused4
      write(lu,rec=068) t1%unused5
      write(lu,rec=069) t1%unused6
      write(lu,rec=070) t1%unused7
      write(lu,rec=071) t1%nzyear
      write(lu,rec=072) t1%nzjday
      write(lu,rec=073) t1%nzhour
      write(lu,rec=074) t1%nzmin
      write(lu,rec=075) t1%nzsec
      write(lu,rec=076) t1%nzmsec
      write(lu,rec=077) t1%nvhdr
      write(lu,rec=078) t1%norid
      write(lu,rec=079) t1%nevid
      write(lu,rec=080) t1%npts
      write(lu,rec=081) t1%internal3
      write(lu,rec=082) t1%nwfid
      write(lu,rec=083) t1%nxsize
      write(lu,rec=084) t1%nysize
      write(lu,rec=085) t1%unused8
      write(lu,rec=086) 4
      write(lu,rec=087) t1%idep
      write(lu,rec=088) t1%iztype
      write(lu,rec=089) t1%unused9
      write(lu,rec=090) t1%iinst
      write(lu,rec=091) t1%istreg
      write(lu,rec=092) t1%ievreg
      write(lu,rec=093) t1%ievtyp
      write(lu,rec=094) t1%iqual
      write(lu,rec=095) t1%isynth
      write(lu,rec=096) t1%imagtyp
      write(lu,rec=097) t1%imagsrc
      write(lu,rec=098) t1%unused10
      write(lu,rec=099) t1%unused11
      write(lu,rec=100) t1%unused12
      write(lu,rec=101) t1%unused13
      write(lu,rec=102) t1%unused14
      write(lu,rec=103) t1%unused15
      write(lu,rec=104) t1%unused16
      write(lu,rec=105) t1%unused17
      write(lu,rec=106) 0!! BECAUSE A XY FILE NOW
      write(lu,rec=107) t1%lpspol
      write(lu,rec=108) t1%lovrok
      write(lu,rec=109) t1%lcalda
      write(lu,rec=110) t1%unused18
      
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

!     * Output the first trace *   
       do i=1,t1 % npts
            write(lu,rec=158+i) t1%trace(i)
      enddo
      
!      print*,maxval(t1%trace(1:t1 % npts))
      
!     * Output the first trace *   
       do i=1,t2 % npts
            write(lu,rec=158+i+t2 % npts) t2%trace(i)
      enddo

 
      close(lu)
      

      return           
      end subroutine write_xyfile
!=======================================================================
