CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C  Rainfall-Runoff-Router-Rootzone Model                                   C
C  Copy Right (C) 1996 - 2017 Wildermuth Environmental, Inc.               C
C                                                                          C
C  This program is free software: you can redistribute it and/or modify    C
C  it under the terms of the GNU General Public License as published by    C
C  the Free Software Foundation, either version 3 of the License, or       C
C  (at your option) any later version.                                     C
C                                                                          C
C  This program is distributed in the hope that it will be useful,         C
C  but WITHOUT ANY WARRANTY; without even the implied warranty of          C
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           C
C  GNU General Public License for more details.                            C
C                                                                          C
C  You should have received a copy of the GNU General Public License       C
C  along with this program.  If not, see <http://www.gnu.org/licenses/>.   C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c 12/14/15 rz_v3S1.for modification
c   use multiple ET data (net) and ET weights
c 10/20/15 rz_v3s.for --- cleaned version of RZ_V3R4B
c 

c open files
c 1   input  *.con response filec
c 2   input        land use data file
c                  rootzone data file
c                  soil data data file
c            *.par HA/LU/SOIL data from runoff model
c 3   input        daily evapotranspiration data
c 4   input  *.pfl HA/LU/SOIL infiltration data from runoff model
c 6   output *.out main output file
c 11  output *.inf monthly infiltration       (af/m) by pervious HA for water budget
c 12  output *.et  monthly evapotranspiration by HA for water budget
c 13  output *.pc  monthly deep percolation   by HA for water budget, input to perc_gm model
c 14  output *.maw monthly M&I applied water  by HA for water budget
c 15  output *.aaw monthly ag applied water   by HA for water budget
c 16  output *.mdw monthly M&I delivered water by HA for water budget
c 17  output *.adw monthly ag delivered water  by HA for water budget
c 18               reserved
c 19               reserved
c 20  output *.sto storage at end of the month
c 21  output *.awl monthly applied water by land use type
c 22  output *.dwl monthly delivered water by land use types


      include 'RF_DIM_3.MAX'

      character*8  haname(MXHA),ch8,fmtin*20
      dimension    haarea(MXHA),           ! area of HA on GW modeling area
     &             xminf(MXHA),            ! monthly total infiltration (af)
     &             xmet(MXHA),             ! monthly evapotranspiration (af)
     &             xmaw(MXHA,2),           ! monthly applied water (af)
     &             xmdw(MXHA,2),           ! monthly delivered water (af)
     &             xmpc(MXHA),             ! monthly percolation (af)
     &             xmsto(MXHA),
     &             xmpcl(MXLU),
     &             xmawl(MXLU),            ! monthly applied water by land use
     &             xmdwl(MXLU),            ! monthly delivered water by land use
     &             xmpcr(MXHA),
     &             xmpca(MXHA),
     &             xmpclr(MXLU),
     &             xmpcla(MXLU),
     &             gwbr(MXHA),              ! v2r3a fraction inside groundwater basin
     &             etwt(MXHA,MXCT),
     &             dget(MXCT)               ! daily gage ET data
      real      rval(12)
      logical*1 apwrq(MXLU),debug          ! applied water required for the landuse
      logical*1 irmn(MXCT,MXLU,12)            ! irrigation month
      dimension apwef(MXCT,MXLU,2),           ! applied water efficiency
     &          rzd(MXCT,MXLU),               !  ( ,1) for root zone efficiency
     &          cropcoef(MXCT,MXLU,12),       !  ( ,2) for delivery efficiency
     &          isrc(MXLU),                ! source of applied water
     &          apad(MXCT,MXLU,3),            ! average percentage of allowable deficit
     &          csm(MXCT,MXST,MXLU,3)         ! critical soil moisture

      dimension ihaiz(MXHA)

      dimension fc(MXST), wp(MXST)   !!!!!, SMC(MXST), AWC(MXST)

      dimension ihls(MXPG,3), dac(MXPG,2),
     *          rinf(MXPG,4), rsm(MXPG,2),
     &          rapw(MXPG,2),rdp(MXPG,2),reta(MXPG,2),
     &          rfc(MXPG,2), rwp(MXPG,2),
     &          rdpr(MXPG,2), rdpa(MXPG,2),
     &          rcsm(MXPG,3,2)

      character*50 file0,ifile,ofile,bl50,clu(21)
      character*1  ch1
      character*20 fmtet

      dimension  ndinm(12)

      data ndinm/31,28,31,30,31,30,31,31,30,31,30,31/




      sqm2ac = 640./1609.**2          ! conversion factor from sq m to acres
      af2cfs = 43560./86400.          !                   from acre-ft/day to cfs
      bl50 = '                                                  '


c***** read input data

      call getcon(file0,n1)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')

      read(1,'(i5)') nsy,ismonth,npor,nha,nlu,nst,ndebug      !!!!!!!,nz
                   ! nsy     starting year for simulaiton (starts Oct. 1)
                   ! npor    number of years to simulate
      read(1,'(2i5)') nz,net
      if(net.le.1) net = 1
      debug = .false.
      if(ndebug.eq.1) debug = .true.
      if(ndebug.eq.2) then
        open(44,file='RZDebug.OUT',status='UNKNOWN')
        write(44,'(a,a)')   'Debug out for ',file0
        write(44,'(a,a,a)')
     &    '                        ',
     &    'Irrigated Pervious Area                   ',
     &    'Non-Irrigated Pervious Area        '
        write(44,'(a,a,a,/)')
     &    'YYYYMMDD IHA ILU IST   K',
     &    ' Infilt Apl-W.   E-ET    D-P     SM  DP-AI',
     &    ' Infilt   E-ET    D-P     SM  DP-AI'
      end if

      ! main output file
      call openofile(1,6,ofile,'RZ Model main output file',file0)
      write(6,'(i5)') nsy,ismonth,npor,nha,nlu,nst,ndebug,nz,net

      ! read zone number for HSA
      call openifile(1,2,ifile)
      read(2,'(a)') fmtin
      do iha=1,nha
        read(2,fmtin) iii,ihaiz(iha)
        if (iii.ne.iha) then
          write(*,'(//,a,a)') ' Problem in ',ifile
          stop
        end if
CCC      write(6,'(2i5)') iha,ihaiz(iha)
      end do
      close(2)
      

      etwt = 1.
      if(net.gt.1)  then
        call openifile(1,2,ifile)
        write(6,'(/,a)') ' ET weights for each HSA'
        write(6,'(a)') '   ID  HANAME   IZ  ET WEIGHTS' 
        read(2,'(a)') fmtet
        do iha=1,nha
          read(2,fmtet) j,ch8,(etwt(iha,i),i=1,net)
          write(6,'(i5,2x,a8,i5,10f8.2)') j,ch8, ihaiz(iha),
     &                (etwt(iha,i),i=1,net)
        end do
        close(2)
      end if        

      ! WELU land use property        !v2
      do iz=1,nz
      call openifile(1,2,ifile)
      write(6,'(/,a,i3)') ' Root zonevariables for zone',iz
      do ilu=1,nlu
c       read(2,'(i6,36x,7x,7x,5f7.0,i7)') j,
        read(2,'(i6,36x,28x,5f7.0,i7)') j,
     &               apad(iz,ilu,1),              ! allowable deficit for SMc
     &               apad(iz,ilu,2),              ! allowable deficit for irrigation to start
     &               apad(iz,ilu,3),              ! allowable deficit for irrigation to end
     &               apwef(iz,ilu,1),             ! irrigation efficiency
     &               apwef(iz,ilu,2),             ! delivery efficiency
     &               isrc(ilu)                    ! source 1 for M&I, 2 for Ag
        write(6,'(i6,5f7.0,i7)') j,
     &               apad(iz,ilu,1),              ! allowable deficit for SMc
     &               apad(iz,ilu,2),              ! allowable deficit for irrigation to start
     &               apad(iz,ilu,3),              ! allowable deficit for irrigation to end
     &               apwef(iz,ilu,1),             ! irrigation efficiency
     &               apwef(iz,ilu,2),             ! delivery efficiency
     &               isrc(ilu)                    ! source 1 for M&I, 2 for Ag

      end do
      read(2,'(a1)') ch1
      read(2,'(i6)') klunv      ! LU type for native vegetation
      read(2,'(i6)') klubar     ! LU type for barren
      read(2,'(i6)') klunrf     ! LU type for no-runoff
      read(2,'(i6)') kluimp     ! LU type for impervious
      write(6,'(a,5i5)') ' klunv, klubar, klunrf, kluimp = ',
     &                    klunv, klubar, klunrf, kluimp 
      close(2)
      end do !iz

      ! WELU rootzone property
      do iz=1,nz
      call openifile(1,2,ifile)
      write(6,'(/,a,i3)') 'Root Zone Parameters for zone',iz

      do ilu=1,nlu

        if(apwef(iz,ilu,1).gt.0.) then
          apwef(iz,ilu,1) = 100./apwef(iz,ilu,1)        ! convert % to ratio
          apwef(iz,ilu,2) = 100./apwef(iz,ilu,2)
          apwrq(ilu) = .true.      ! applied water required
          if(isrc(ilu).ne.1 .and. isrc(ilu).ne.2) then
            write(*,'(//,a,a,a)') ' Error in ',ifile,' -- Stop'
            stop
          end if
        else
          apwrq(ilu) = .false.
          isrc(ilu) = 0
        end if

        read(2,'(i8,a50,20f8.0)') j,clu(ilu),
     &                            rzd(iz,ilu),
     &                            (cropcoef(iz,ilu,icm),icm=1,12)
        if(j.ne.ilu) stop
        do icm=1,12
          irmn(iz,ilu,icm) = .false.
          if(cropcoef(iz,ilu,icm).lt.0.) then
            cropcoef(iz,ilu,icm) = -cropcoef(iz,ilu,icm)
          else
            if(apwrq(ilu)) irmn(iz,ilu,icm) = .true.
          end if
        end do

        write(6,'(i6,2x,a50,f8.2,16f8.2)') j,clu(ilu),
     &                             rzd(iz,ilu),
     &                             (cropcoef(iz,ilu,icm),icm=1,12)

        rzd(iz,ilu) = rzd(iz,ilu)*12.            ! convert (ft) to (inches)
      end do
      close(2)
      end do !iz

      ! soil property
      call openifile(1,2,ifile)

      write(6,'(/,a)') 'Soil Parameters'
      do ist=1,nst
        read(2,'(i8,3f8.0)') j,fc(ist),wp(ist)
        if(j.ne.ist) stop
        write(6,'(i8,3f8.3)') j,fc(ist),wp(ist)
        do iz=1,nz
        do ilu=1,nlu
          csm(iz,ist,ilu,1) = wp(ist)
     &                     + (fc(ist) - wp(ist))*apad(iz,ilu,1)/100.
          csm(iz,ist,ilu,2) = wp(ist)
     &                     + (fc(ist) - wp(ist))*apad(iz,ilu,2)/100.
          csm(iz,ist,ilu,3) = wp(ist)
     &                     + (fc(ist) - wp(ist))*apad(iz,ilu,3)/100.
        end do
        end do !iz
      end do
      close(2)

      do iz=1,nz
      write(6,'(a,i2)') 'Zone = ',iz
      write(6,'(a)') 'SOIL  LU   RZD    WP    FC   SMC   SMI   SMC'
      do ilu=1,nlu
        do ist=1,nst
          write(6,'(2i4,f6.1,10f6.2)')
     &              ist,ilu,rzd(iz,ilu),wp(ist),fc(ist),
     &              (csm(iz,ist,ilu,i),i=1,3)
        end do
      end do
      end do !iz
      ! csm(,,1) --- critical soil moisture
      ! csm(,,2) --- soil moisture to start irrigation
      ! csm(,,3) --- soil moisture to end irrigation

      ! HA-LU-ST data
      call openifile(1,2,ifile)
      write(6,'(/,a,a)') 'File = ',ifile
      read(2,'(i8)') iha
      if(iha.ne.nha) stop
      read(2,'(25a8)') (haname(iha),iha=1,nha)

      read(2,'(i8)') kpg
      if(kpg.gt.MXPG) then
        write(*,'(//,a//)') ' Increase dimension parameter MXPG'
        stop
      end if
      
      read(2,'(50i4)') (ihls(kp,1),kp=1,kpg)      ! iha
      read(2,'(50i4)') (ihls(kp,2),kp=1,kpg)      ! ilu
      read(2,'(50i4)') (ihls(kp,3),kp=1,kpg)      ! ist
c     read(2,'(25f10.0)') (dac(kp),kp=1,kpg)      ! dac(iha,ilu,ist)
      read(2,'(25f10.0)') (dac(kp,1),dac(kp,2),kp=1,kpg)
      close(2)

      write(*,'(a,a)') ' done with ',ifile !v2r3a delete

      gwbr = 1.

      write(6,'(a,i5)') 'NHA = ',nha
      write(6,'(a,i5)') 'KPG = ',kpg

      do iha=1,nha
        haarea(iha) = 0.
        xmpc(iha)   = 0.      ! to hold pervious area temporarily
        xmaw(iha,1) = 0.      ! to hold pervious M&I area temporarily
        xmaw(iha,2) = 0.      ! to hold pervious AG area temporarily
      end do

      write(6,'(/,a)') 'Polygons with area < 0.1 acres will be'
      write(6,'(  a)') 'removed to avoid numerical problem'

      do ilu=1,nlu
         xmpcl(ilu) = 0.
      end do

      do kp=1,kpg
        iha = ihls(kp,1)
        ilu = ihls(kp,2)
        ist = ihls(kp,3)
        iz = ihaiz(iha)
        write(6,'(i6,2x,a8,3i4,2f10.2)') kp,haname(iha),
     &                        iha,ilu,ist,dac(kp,1),dac(kp,2)
        if(dac(kp,1).lt.0.1) then
          if(dac(kp,1).gt.0.) write(6,'(16x,a)') 'Removed 1'
          dac(kp,1) = 0.
        end if
        if(dac(kp,2).lt.0.1) then
          if(dac(kp,2).gt.0.) write(6,'(16x,a)') 'Removed 2'
          dac(kp,2) = 0.
        end if

        if(ist.eq.0) then
          rfc(kp,1) = 0.
          rwp(kp,1) = 0.
          rcsm(kp,1,1) = 0.
          rcsm(kp,2,1) = 0.
          rcsm(kp,3,1) = 0.
          rfc(kp,2) = 0.
          rwp(kp,2) = 0.
          rcsm(kp,1,2) = 0.
          rcsm(kp,2,2) = 0.
          rcsm(kp,3,2) = 0.
        else
          ! irrigated pervious area
          rfc(kp,1) = rzd(iz,ilu) * fc(ist)
          rwp(kp,1) = rzd(iz,ilu) * wp(ist)
          rcsm(kp,1,1) = rzd(iz,ilu)*csm(iz,ist,ilu,1)             ! allowable critical soil moisture content
          rcsm(kp,2,1) = rzd(iz,ilu)*csm(iz,ist,ilu,2)             ! soil moisture content to start irrigation
          rcsm(kp,3,1) = rzd(iz,ilu)*csm(iz,ist,ilu,3)             ! soil moisture content to end irrigation
          ! non-irrigated pervious area
          rfc(kp,2) = rzd(iz,klubar) * fc(ist)
          rwp(kp,2) = rzd(iz,klubar) * wp(ist)
          rcsm(kp,1,2) = rzd(iz,klubar)*csm(iz,ist,klubar,1)
          rcsm(kp,2,2) = rzd(iz,klubar)*csm(iz,ist,klubar,2) ! this is useless, since no irrigation area
          rcsm(kp,3,2) = rzd(iz,klubar)*csm(iz,ist,klubar,3) ! this is useless, since no irrigation area
        end if

        haarea(iha) = haarea(iha) + dac(kp,1) + dac(kp,2)
        xmpcl(ilu) = xmpcl(ilu) + (dac(kp,1) + dac(kp,2))*gwbr(iha)
        if(ist.gt.0) then        ! don't need this if
          xmpc(iha) = xmpc(iha) + dac(kp,1) + dac(kp,2)
          if(apwrq(ilu)) then
            rsm(kp,1) = rfc(kp,1)                     ! initial soil moisture
            xmaw(iha,isrc(ilu)) = xmaw(iha,isrc(ilu)) + dac(kp,1)
          else
            rsm(kp,1) = (rwp(kp,1) + rfc(kp,1))/2.
          end if
          rsm(kp,2) = (rwp(kp,2) + rfc(kp,2))/2.    ! initial soil moisture
        end if                   ! don't need this if
      end do ! kp

      write(6,'(/,a)') 'Area summarized by HA'
      write(6,'(a)') '    ID  HANAME     TOTAL   XMAW_1    XMAW_2'
      do iha=1,nha
c        if(haarea(iha).gt.0.) write(6,'(i6,2x,a8,f10.2)') iha,
c     &                        haname(iha),haarea(iha)
        if(haarea(iha).gt.0.) write(6,'(i6,2x,a8,10f10.2)') iha,
     &             haname(iha),haarea(iha),xmaw(iha,1),xmaw(iha,2)
      end do

      do ilu=1,nlu
        write(6,*) ilu,isrc(ilu),apwrq(ilu)
      end do

      ! HA-LU-ST daily infiltration data
      call openifile(1,4,ifile)
      read(4,'(i4,2i2))') m11,m21,m31
      if(m11.ne.nsy .or. m21.ne.ismonth .or. m31.ne.1) then
        write(*,'(a,a50,i4,2i2.2)') ' Start date of file, ',ifile,
     &                              m11,m21,m31
        write(*,'(a,i4,2i2.2)') ' Simulation starts ',nsy,ismonth,1
        stop
      end if
      backspace(4)


      ! Daily evapotranspiration data
      call openifile(1,3,ifile)
      read(3,'(a20)') fmtet

      ! set the file at the simulation starting date
      do while (.true.)
        read(3,fmtet,end=990) m11,m21,m31
        if(m11.eq.nsy .and. m21.eq.ismonth .and. m31.eq.1) then
          backspace(3)
          goto 40
        end if
      end do
      write(*,'(//,a,/)') ' Check Reference ETo data file.'
      stop
 40   continue

      read(1,'(a)') ch1

      call openofile(1,11,ofile,'Monthly infiltration to RZ',file0)
      call openofile(1,12,ofile,'Monthly evapotranspiration',file0)
      call openofile(1,13,ofile,'Monthly deep percolation  ',file0)

      call openofile(1,14,ofile,'Monthly M&I applied water to HA',
     &                                   file0)
      call openofile(1,15,ofile,'Monthly ag applied water to HA',
     &                                   file0)
      call openofile(1,16,ofile,'Monthly M&I delivered water to HA',
     &                                   file0)
      call openofile(1,17,ofile,'Monthly ag delivered water to HA',
     &                                   file0)
      read(1,'(a1)') ch1
      read(1,'(a1)') ch1
      call openofile(1,20,ofile,'Storage in Rootzone at the end',file0)
      call openofile(1,21,ofile,'Monthly applied by landuse',file0)
      call openofile(1,22,ofile,'Monthly delivered water by landuse',
     &                                   file0)

      do iu=11,20
        if(iu.ne.18 .and. iu.ne.19)
     &  write(iu,'(a6,25a8,/,(6x,25a8))') 'yyyymm',
     &                (haname(iha),iha=1,nha)
      end do
      write(13,'(a6,25f8.1,/,(6x,25f8.1))') 'Area T',
     &                (haarea(iha),iha=1,nha)
      write(13,'(a6,25f8.1,/,(6x,25f8.1))') 'Area P',
     &                (xmpc(iha),iha=1,nha)
      write(14,'(a6,25f8.1,/,(6x,25f8.1))') 'Area P',
     &                (xmaw(iha,1),iha=1,nha)
      write(15,'(a6,25f8.1,/,(6x,25f8.1))') 'Area P',
     &                (xmaw(iha,2),iha=1,nha)

      do iu=21,22
        write(iu,'(a)') 'Groundwater basin area only'
        write(iu,'(a)') 'WE Land Use Types'
        write(iu,'(a6,25i8)') 'Date  ',(ilu,ilu=1,nlu)
        write(iu,'(a,25f8.0)') 'Area  ',(xmpcl(ilu),ilu=1,nlu)
      end do

      if(debug) then
        call openofile(1,91,ofile,'Monthly Deep Percolation',
     &                                   file0)
        call openofile(1,92,ofile,'Monthly Deep Perc. of Rainfall',
     &                                   file0)
        call openofile(1,93,ofile,'Monthly Deep Perc. of App. Water',
     &                                   file0)
        do iu=91,93
          write(iu,'(a)') 'Groundwater basin area only'
          write(iu,'(a,25i8)')   'YearMn',(ilu,ilu=1,nlu)
          write(iu,'(a,25f8.0)') 'Area  ',(xmpcl(ilu),ilu=1,nlu)
        end do

        read(1,'(a8,2i5)') ch8,ist,ilu
        do i=1,nha
          if(ch8.eq.haname(i)) then
            iha = i
            goto 50
          end if
        end do
        write(*,'(//,a,a,a)') ' HA ID ',ch8,' is not available, stop'
        stop
 50     jpg = 0
        do k=1,kpg
          if(iha.eq.ihls(k,1) .and. ilu.eq.ihls(k,2) .and.
     &       ist.eq.ihls(k,3)) then
            jpg = k
            goto 55
          end if
        end do
        write(*,'(//,a,a8,2i7,a)') ' debugging polygon ',ch8,ist,ilu,
     &                             ' is not available, stop'
        stop
 55     read(1,'(i5)') idsy
        read(1,'(i5)') idey
        call openofile(1,94,ofile,'Debug out',file0)
        write(94,'(a,a8)') 'DA Name ',ch8
        write(94,'(a,i5)') 'DA ID   ',iha
        write(94,'(a,i5)') 'Landuse ',ilu
        write(94,'(a,i5)') 'Soil    ',ist
        write(94,'(a,i5)') 'DA #    ',jpg
        write(94,'(a8,a)') 'YearMnDy',
     &           '  Infilt      AW      ET      DP Storage'

      end if  !debug


      ! initial rootzone soil moisture storage
      do iha=1,nha
        xmsto(iha) = 0.
      end do
      do k=1,kpg
        iha = ihls(k,1)
        ilu = ihls(k,2)
        ist = ihls(k,3)
        if(ist.gt.0) then  ! ignore impervious area
          xmsto(iha) = xmsto(iha) +
     &                 (rsm(k,1)*dac(k,1) + rsm(k,2)*dac(k,2)) / 12.
c v3      xmsto(iha) = xmsto(iha) + rsm(k) * dac(k) / 12.
        end if
      end do
      icy = nsy
      icm = ismonth-1
      write(20,8081) icy,icm,(xmsto(iha),iha=1,nha)


      do 500 ipy=1,npor !===== YEAR LOOP ===============================

      icy = nsy+ipy-1     ! starts at January 1
      write(*,'(/,i5)') icy


      do 400 iwm=1,12   !===== MONTH LOOP ==============================

C      write(*,'(5x,i2)') iwm
      icm = iwm + ismonth - 1
      if(icm.gt.12) then
        icm = icm - 12
        if(icm.eq.1) icy = icy+1
      end if

      mnday = ndinm(icm)
      if(icm.eq.2) then
        if(mod(icy,4).eq.0.) mnday = 29
      end if
      do iha=1,nha
        xminf(iha) = 0.
        xmet(iha) = 0.
        xmaw(iha,1) = 0.         ! 1 M&I
        xmaw(iha,2) = 0.         ! 2 AG
        xmdw(iha,1) = 0.
        xmdw(iha,2) = 0.
        xmpc(iha) = 0.
c       xmimp(iha) = 0.
        xmsto(iha) = 0.
        xmpcr(iha)   = 0.
        xmpca(iha)   = 0.
      end do
      do ilu=1,nlu
        xmpcl(ilu) = 0.
        xmawl(ilu) = 0.
        xmdwl(ilu) = 0.
        xmpclr(ilu) = 0.
        xmpcla(ilu) = 0.
      end do

      do 300 idy=1,mnday    !===== DAY LOOP ============================


      read(3,fmtet,end=990) m10,m20,m30,(dget(i),i=1,net)
        
      read(4,'(i4,2i2,/,(25f8.0))',end=990) m11,m21,m31,
     &                    ((rinf(k,i),i=1,4),k=1,kpg)  ! (inches/day)
                ! rinf(k,1) effective rainfall to irrigated pervious area
                ! rinf(k,2)      "       "     to nonirrigated pervious area (barren area)
                ! rinf(k,3) infiltration       in irrigated pervious area
                ! rinf(k,4)      "             in nonirrigated pervious area (barren area)

      do 200 k=1,kpg   !===== HA/LU/ST LOOP ============================
        iha = ihls(k,1)
        ilu = ihls(k,2)
        ist = ihls(k,3)
        iz  = ihaiz(iha)

        RefETo = 0.
        do i=1,net
          RefETo = RefETo + etwt(iha,i) * dget(i)
        end do

        if(dac(k,1).eq.0.) goto 100
        etp = cropcoef(iz,ilu,icm) * RefETo   ! maximum ET potential
        rsmn = rsm(k,1) + rinf(k,3)
        rapw(k,1) = 0.
        rapw(k,2) = 0.
        rdp(k,1) = 0.

        rdpr(k,1) = 0.
        rdpa(k,1) = 0.
        rapw(k,1) = 0.
        rapw(k,2) = 0.
        deficit = 0.

        if(rsmn.gt.rfc(k,1)) then  ! land is oversaturated with rain
          rdpr(k,1) = rsmn - rfc(k,1)
          rsmn = rfc(k,1)
          etd = etp
        else
          if(apwrq(ilu).and.irmn(iz,ilu,icm)) then
            if(rsmn.lt.rcsm(k,2,1)) then    ! soil is dry below irrigation requirement soil moisture
              deficit = rcsm(k,3,1) - rsmn
              rapw(k,1) = deficit * apwef(iz,ilu,1)
              rapw(k,2) = rapw(k,1) * apwef(iz,ilu,2)
              rsmn = rsmn + rapw(k,1)
              rdpa(k,1) = rsmn - rcsm(k,3,1)
              rsmn = rcsm(k,3,1)
            end if
          end if
        end if

        if(rsmn.gt.rcsm(k,1,1)) then
          etd = etp
        elseif(rsmn.gt.rwp(k,1)) then
          etd = etp * (rsmn-rwp(k,1)) / (rcsm(k,1,1) - rwp(k,1))
        else
          etd = 0.
        end if

        rsmn = rsmn - etd             ! evapotranspiration

        rsm(k,1) = rsmn
        reta(k,1) = etd               ! actual et
        rdp(k,1) = rdpr(k,1) + rdpa(k,1)    ! deep perc out of root zone

        ! monthly summary
                                                                ! monthly
        xminf(iha) = xminf(iha) + rinf(k,3) * dac(k,1) / 12.
        v1 = rdp(k,1) * dac(k,1) / 12.
        xmpc(iha) = xmpc(iha) + v1     
        xmpcl(ilu) = xmpcl(ilu) + v1*gwbr(iha) 

        v1 = rdpr(k,1) * dac(k,1) / 12.
        xmpcr(iha) = xmpcr(iha) + v1
        xmpclr(ilu) = xmpclr(ilu) + v1*gwbr(iha)

        v1 = rdpa(k,1) * dac(k,1) / 12.
        xmpca(iha) = xmpca(iha) + v1
        xmpcla(ilu) = xmpcla(ilu) + v1*gwbr(iha)

        v1 = rapw(k,1) * dac(k,1) / 12.
        v2 = rapw(k,2) * dac(k,1) / 12.

        if(v1.gt.0.) then
          xmaw(iha,isrc(ilu)) = xmaw(iha,isrc(ilu)) + v1  
          xmdw(iha,isrc(ilu)) = xmdw(iha,isrc(ilu)) + v2 
          xmawl(ilu) = xmawl(ilu) + v1*gwbr(iha) 
          xmdwl(ilu) = xmdwl(ilu) + v2*gwbr(iha)  
        end if

  90    xmet(iha) = xmet(iha) + reta(k,1) * dac(k,1) / 12.          ! ET by HA

 100    continue


        if(dac(k,2).eq.0.) goto 195

        etp = cropcoef(iz,klubar,icm) * RefETo   ! maximum ET potential
        rsmn = rsm(k,2) + rinf(k,4)                                   ! 0329 rsm(k) need to be adjusted every year
        rdp(k,2) = 0.

        rdpr(k,2) = 0.           ! deep perc of rain water
        rdpa(k,2) = 0.           ! deep perc of applied water
        deficit = 0.

        if(rsmn.gt.rfc(k,2)) then  ! land is oversaturated with rain
          rdpr(k,2) = rsmn - rfc(k,2)
          rsmn = rfc(k,2)
          etd = etp
        end if

        if(rsmn.gt.rcsm(k,1,2)) then
          etd = etp
        elseif(rsmn.gt.rwp(k,2)) then
          etd = etp * (rsmn-rwp(k,2)) / (rcsm(k,1,2) - rwp(k,2))
        else
          etd = 0.
        end if

        rsmn = rsmn - etd               ! evapotranspiration

        rsm(k,2) = rsmn
        reta(k,2) = etd                 ! actual et

        rdp(k,2) = rdpr(k,2)

        ! monthly summary
        xminf(iha) = xminf(iha) + rinf(k,4) * dac(k,2) / 12.    ! infiltration by HA
        v1 = rdp(k,2) * dac(k,2) / 12.
        xmpc(iha) = xmpc(iha) + v1                              ! deep percolation by HA
        xmpcl(ilu) = xmpcl(ilu) + v1*gwbr(iha)    !v2r3a

        v1 = rdpr(k,2) * dac(k,2) / 12.
        xmpcr(iha) = xmpcr(iha) + v1
        xmpclr(ilu) = xmpclr(ilu) + v1*gwbr(iha)  !v2r3a

        v1 = rdpa(k,2) * dac(k,2) / 12.
        xmpca(iha) = xmpca(iha) + v1
        xmpcla(ilu) = xmpcla(ilu) + v1*gwbr(iha)  !v2r3a


 190    xmet(iha) = xmet(iha) + reta(k,2) * dac(k,2) / 12.       ! ET by HA

 195    continue
 
        if(ndebug.eq.2) then
          rv1 = rdp(k,1) * dac(k,1)
          rv2 = rdp(k,2) * dac(k,2)
          rv3 = rdpr(k,1)+rdpa(k,1)
          write(44,'(i4,2i2.2,4i4,5f7.4,f7.2,4f7.4,f7.2)')
     &      icy,icm,idy,iha,ilu,ist,k,
     &      rinf(k,3),rapw(k,1),reta(k,1),rv3,rsm(k,1),rv1,
     &      rinf(k,4),reta(k,2),rdpr(k,2),rsm(k,2),rv2
        end if


 200  continue ! k - HA/LU/ST loop
        if(debug .and.
     &     icy.ge.idsy .and. icy.le.idey)
     &       write(94,'(i4,2i2,10f8.2)') icy,icm,idy,
     &                rinf(jpg,3),
     &           rapw(jpg,1),reta(jpg,1),rdp(jpg,1),rsm(jpg,1),
     &                rinf(jpg,4),
     &           rapw(jpg,2),reta(jpg,2),rdp(jpg,2),rsm(jpg,2)

 300  continue ! day loop


      ! rootzone soil moisture storage
      do k=1,kpg
        iha = ihls(k,1)
        ilu = ihls(k,2)
        ist = ihls(k,3)
        if(ist.gt.0) then  ! ignore impervious area
          xmsto(iha) = xmsto(iha) +
     &          (rsm(k,1) * dac(k,1) + rsm(k,2) * dac(k,2)) / 12.
        end if
      end do

      ! print monthly variables
      write(11,8081) icy,icm,(xminf(iha),iha=1,nha)
      write(12,8081) icy,icm,(xmet(iha),iha=1,nha)
      write(13,8081) icy,icm,(xmpc(iha),iha=1,nha)
      write(14,8081) icy,icm,(xmaw(iha,1),iha=1,nha)
      write(15,8081) icy,icm,(xmaw(iha,2),iha=1,nha)
      write(16,8081) icy,icm,(xmdw(iha,1),iha=1,nha)
      write(17,8081) icy,icm,(xmdw(iha,2),iha=1,nha)


      write(20,8081) icy,icm,(xmsto(iha),iha=1,nha)
      write(21,8081) icy,icm,(xmawl(ilu),ilu=1,nlu)
      write(22,8081) icy,icm,(xmdwl(ilu),ilu=1,nlu)

      if(debug) then
        write(91,8081) icy,icm,(xmpcl(ilu),ilu=1,nlu)
        write(92,8081) icy,icm,(xmpclr(ilu),ilu=1,nlu)
        write(93,8081) icy,icm,(xmpcla(ilu),ilu=1,nlu)
      end if
 400  continue ! month loop

 500  continue ! year loop

 990  close(1)
      close(3)
      close(4)
      close(6)
      do iu=11,22
        close(iu)
      end do

 8081 format(i4,i2.2,25f8.2,/,(6x,25f8.2))

      stop
      end

