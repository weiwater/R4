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

! rf_v7a.for   1/21/2016
!    wspi2 and wspn2 in infiltration is returned to rf_v6.for
! 12/21/2015
! rf_v7 is modified from rf_v6s, to extract initial abstract on pervious area infiltration
!
! 10/20/15 rf_v6s.for  --- cleaned version of rf_v6r, with larger dimensions
!

!
! Iunit  Input/Out *.EXT  Description
!  1     Input            control file *** not closed - modify
!  2     Input            HA definition and precipitation station weights
!  2     Input            for each HA, soil-kN table ID number, CN multiplier and rain multiplier
!  2     Input            Soil Curve Number Tables
!  2     Input            HA-LU-Soil area data
!  3     Input            daily rainfall data
!  4     Input            daily evaporation data
!  5     Input            daily evapotranspiration data
!  2     Input            if qopt, runoff-TDS and runoff-TIN tables
!  6     Output           main output file

!  11    Output    *.PAR  HA ID, LU ID, ST ID, and area for deep perc model
!  12    Output    *.PFL  daily perc data to root zone from HA-LU-ST area

!  21    Output           daily runoff data from HA for router model
!  22    Output           if qopt, runoff TDS conc, input to route model
!  23    Output           if qopt, runoff TIN conc, input to route model

!  31    Output    *.MPR  HA area, xmpr, in inches
!  34    Output    *.MRO  HA area, xmro, in inches
!  35    Output    *.MIF  HA area, xmpci, in inches
!  37    Output    *.MEV  HA area, xmev, in inches
!  38    Output    *.MET  HA area, xmet, in inches



      include 'rf_dim_3.max'           !  maximum dimensions

      character*50 file0,ifile,ofile,bl50
      character*40 fmtrn, fmtev, fmtet,fmtin
      character    ch1*1, ch8*8

      ! property for area (iha,ilu,ist)
      common /havar/haname(MXHA),haarea(MXHA)
      character*8  haname

      common /big_ar/ cn1(MXHA,MXLU,4),      ! CN for AMC I
     &          cn2(MXHA,MXLU,4),            ! CN for AMC II
     &          cn3(MXHA,MXLU,4),            ! CN for AMC III
     &          ssc(MXHA,MXLU,4),            ! soil storage capacity
     &          sscmax(MXHA,MXLU,4),
     &          sscmin(MXHA,MXLU,4),
     &          dacls(MXHA,MXLU,4),          ! DA total for LU and soil     
     &          dacid(MXHA,MXLU,4),          ! DCIA area of (iha,ilu,ist)   
     &          daciu(MXHA,MXLU,4),          ! UCIA area                    
     &          dacpi(MXHA,MXLU,4),          ! pervious irrigated area      
     &          dacpn(MXHA,MXLU,4),          ! pervious non-irrigated area  
     &          cnpi(MXHA,MXLU,4),           ! CN of pervious irrigated area
     &          cnpn(MXHA,MXLU,4),           ! CN of pervious non-irrigated 
     &          sspi(MXHA,MXLU,4),           ! S of pervious irrigated area 
     &          sspn(MXHA,MXLU,4),           ! S of pervious non-irrigated  
     &          soilcn(MXCT,MXLU,4),         ! landuse - soil curve number
     &          fimp(MXCT,MXLU),             ! fraction of impervious area for land use type ilu
     &          firg(MXCT,MXLU),             ! fraction of irrigated area of pervious area of irrigated land use  RF_V2
     &          fDCIA(MXCT,MXLU),            ! fraction of DCIA
     &          fUCIAi(MXCT,MXLU),           ! fraction of runoff from UCIA to irrigated area 
     &          CNi(MXCT),                   ! impervious area CN
     &          Si(MXCT),                    ! impervious area S
     &          rIai(MXCT)                   ! impervious area Ia

      ! property for area(iha)
      dimension haluarea(MXHA, MXLU)   ! used for summary report purpose only

      dimension kct(MXHA),             ! cn curve number table
     &          rmcn(MXHA),            ! multiplier for sensitivity analysis
     &          rmrain(MXHA,2),        ! multiplier for sensitivity analysis
     &          roha(MXHA),            ! runoff from HSA
     &          rmcz(10  ),            ! multiplier for CN in CN zone
     &          rmrz(10  ,2)           ! multiplier for rain in rain zone
                                       ! rainfall threshold in rain zone
      dimension gmarea(MXHA),
     &          gmfrac(MXHA)           ! fraction of HA area overlying GW basin

      dimension prwt(MXHA,MXPR),       ! weight for precipitation stations
     &          prst(MXPR),            ! precipitation recorded at rain gage
     &          prha(MXHA)             ! precipitation on HA

      dimension evwt(MXHA,MXEV),       ! weight for evaporation stations
     &          evst(MXEV),            ! daily pan evaporation data
     &          evha(MXHA)             ! weighted evaporation data for iha area

      dimension etwt(MXHA,MXEV),       ! ETp weight
     &          etst(MXEV),            ! daily ETp at Cimis Station
     &          etha(MXHA)             ! weighted ETp for iha area

      dimension cha(MXHA,2),           ! concentration in runoff 1 for TDS, 2 for TIN
     &          ntconc(0:MXLU,2),
     &          tconc(0:MXLU,2,MXTT,2)
                !  subscript 1 - for welu
                !            2 - 1 for TDS, 2 for TIN
                !            3 - table entry, serial
                !            4 - 1 for runoff, 2 for conc

      ! monthly summary data
      dimension xmro(MXHA), xmpr(MXHA), xmev(MXHA), xmpci(MXHA),
     &      xmet(MXHA)

      dimension iv(MXPG,3), rv(MXPG,4),rval(10)
                ! used once to record          
                ! rv(,1) - dacpi(iha,ilu,ist)  
                ! rv(,2) - dacpn(iha,ilu,ist)                 
                !
                ! used daily to record         
                ! rv(,1) - effective rainfall to dacpi(iha,ilu,ist)
                ! rv(,2) - effective rainfall to danpn 
                ! rv(,3) - infiltration in dacpi       
                ! rv(,4) - infiltration in dacpn       

      logical*1 qopt                   ! not activated for this version
                ! 1 for TDS
                ! 2 for TIN




      sqm2ac = 640./1609.**2           ! conversion factor from sq m to acres
      af2cfs = 43560./86400.           !                   from acre-ft/day to cfs
      bl50 = '                                                  '
      dacmin = 0.1                     ! ignore polygon noise

      ! initialize variables
      do iha=1,MXHA
        haname(iha) = '        '
        haarea(iha) = 0.
        rmcn(iha) = 1.
        rmrain(iha,1) = 1.
        rmrain(iha,2) = 5.
        kct(iha) = 1

        gmarea(iha) = 0.
        gmfrac(iha) = 0.

        do ilu=1,MXLU
          haluarea(iha,ilu) = 0.
          do ist = 1,4
            ! default values
            cn1(iha,ilu,ist) = 60.
            cn2(iha,ilu,ist) = 75.
            cn3(iha,ilu,ist) = 90.

            ssc(iha,ilu,ist) = 0.
            sscmax(iha,ilu,ist) = 0.
            sscmin(iha,ilu,ist) = 0.

            dacls(iha,ilu,ist) = 0.   ! DA total for LU and soil   
            dacid(iha,ilu,ist) = 0.   ! DCIA area of (iha,ilu,ist) 
            daciu(iha,ilu,ist) = 0.   ! UCIA area
            dacpi(iha,ilu,ist) = 0.   ! pervious irrigated area    
            dacpn(iha,ilu,ist) = 0.   ! pervious non-irrigated area

            cnpi(iha,ilu,ist) = 0.    ! CN of pervious irrigated area
            cnpn(iha,ilu,ist) = 0.    ! CN of pervious non-irrigated 
            sspi(iha,ilu,ist) = 0.    ! S of pervious irrigated area 
            sspn(iha,ilu,ist) = 0.    ! S of pervious non-irrigated  

          end do
        end do
      end do

      call getcon(file0,n2)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')

      read(1,'(i5)') nha,nlu,nsy,ismonth,npor,
     &               nprwts,nevwts,netwts,nct
                   ! nha        number of hydraulic area to compute runoff
                   ! nlu        number of land use types
                   ! nsy        starting year of runoff histories
                   ! ismonth    starting month
                   ! npor       number of years in simulation
                   ! nprwts     number of precipitation stations
                   ! nevwts     numer of evaporation stations
                   ! netwts     numer of ETp stations
                   ! nct        numer of soil-curve number tables
      qopt = .false.
c     read(1,'(f5.0)') rval(1)
      read(1,'(2f5.0)') rval(1),rval(2)
      if(rval(1).ne.0.) qopt = .true.
      idebug = rval(2)

      ! Check dimensions
      if(nprwts.lt.nha .and. nprwts.gt.MXPR) then
        write (*,'(//,a,//)') ' Check MXPR'
        stop
      end if

      if(nha.gt.MXHA .or.
     &   nlu.gt.MXLU .or.
     &   nevwts.gt.MXEV .or.
     &   netwts.gt.MXEV .or.
     &   npor.gt.MXYR) then
        write (*,'(//,a,//)') ' Check maximum dimensions.'
        stop
      end if


      call openofile(1,6,ofile,'Runoff main output file',file0)
      write(6,'(i5)') nha,nlu,nsy,ismonth,npor,nprwts,nevwts,netwts,nct
      write(6,'(f5.0)') rval(1)

      ! read HA name and precipitation station weighting factors
      ! note that the order of precipitation stations in this file
      ! should be same as in daily precipitation data file
      if(nprwts.lt.nha) then
        write(6,'(a)') 'Precipitation gaging station weights'
        write(6,'(/,5x,a8,48i5)') 'HA Name ',(j,j=1,nprwts), 
     &                                     (j,j=1,nevwts),   
     &                                     (j,j=1,netwts)    
      else 
        write(6,'(/,5x,a8,48i5)') 'HA Name ',(j,j=1,nevwts), 
     &                                       (j,j=1,netwts)  
      end if                                                 
      call openifile(1,2,ifile)
      read(2,'(a40)') fmtin
      do iha=1,nha
        if(nprwts.lt.nha) then
          read(2,fmtin,err=992,end=992) k,
     &             haname(iha),
     &             (prwt(iha,j),j=1,nprwts),
     &             (evwt(iha,j),j=1,nevwts),
     &             (etwt(iha,j),j=1,netwts)
          call ctrim8r(haname(iha))
          write(6,'(i4,1x,a8,48f5.2)')
     &             iha, haname(iha),
     &             (prwt(iha,j),j=1,nprwts),
     &             (evwt(iha,j),j=1,nevwts),
     &             (etwt(iha,j),j=1,netwts)
        else
          read(2,fmtin,err=992,end=992) k,
     &             haname(iha),
     &             (evwt(iha,j),j=1,nevwts),
     &             (etwt(iha,j),j=1,netwts)
          call ctrim8r(haname(iha))
          write(6,'(i4,1x,a8,48f5.2)')
     &             iha, haname(iha),
     &             (evwt(iha,j),j=1,nevwts),
     &             (etwt(iha,j),j=1,netwts)
        end if
      end do
      close(2)

      ! For each HA, read area, soil-CN table ID number, CN multiplier and
      ! and rain multiplier.
      ! Multiple CN tables may be needed if the soil surveys for
      ! different counties were done by different geologists and
      ! classifications are slightly different.
      ! Multipliers for CN and rainfall are for sensitivity analysis.

      call openifile(1,2,ifile)
      write(6,'(/,a,a)') 'Reading ',ifile
      read(2,'(a40)') fmtin
      do j=1,nha
        read(2,fmtin,err=992,end=992) ii,ch8,
     &             areaha, ii1,ii2,rval(3)
        ! find HAID
        call ctrim8r(ch8)
        do i=1,nha
          iha=i
          if(ch8.eq.haname(i)) goto 5
        end do
        write(*,'(//,a,a,a,//)') ' HA ID ',ch8, ' is not available'
        stop
 5      kct(iha) = ii1
        if(ii1.gt.nct) then
          write(*,'(//,a)') ' Kct > Nct'
          stop
        end if
        haarea(iha) = areaha * sqm2ac    !v1
        iv(iha,1) = ii1                    ! CN zone number
        iv(iha,2) = ii2                    ! rain zone number
        gmfrac(iha) = rval(3)
        write(6,'(i4,1x,a8,f10.0,i8,4f8.3)')
     &        iha,haname(iha),haarea(iha), kct(iha),
     &        rmcn(iha),rmrain(iha,1),rmrain(iha,2),gmfrac(iha)
      end do

      close(2)


      call openifile(1,2,ifile)
      write(6,'(/,a,a)') 'Reading ',ifile 
      read(2,'(i8)') ncnz        ! number of CN zones 
      if(ncnz.ne.nct) then
        write(*,'(a)') ' Con file and Parameter files do not match'
        write(*,'(a)') ' for number fo CN zones'
        stop
      end if
      if(ncnz.gt.MXPR) then
        write(*,'(a)') ' Dimension problem for ncnz'
        stop
      end if
      do i=1,ncnz 
        read(2,'(i8,f8.0)') j,rmcz(i) 
      end do
      call skipline(2,'*')
      read(2,'(i8)') nrainz        ! number of rain zones
      do i=1,nrainz
        read(2,'(i8,2f8.0)') j,rmrz(i,1),rmrz(i,2)
      end do 
      close(2)

      do iha=1,nha
        rmcn(iha) = rmcz(iv(iha,1))              ! multiplier
        rmrain(iha,1) = rmrz(iv(iha,2),1)        ! multiplier
        rmrain(iha,2) = rmrz(iv(iha,2),2)        ! threshold value
      end do


      write(6,'(/,a)') 'Soil Curve Number Tables'
      do ict=1,nct        ! read % impervious and CNs for WELU
        call openifile(1,2,ifile)
        do ilu=1,nlu
          read(2,'(i6,36x,4f7.0)',err=992,end=992)
     &             j, (soilcn(ict,ilu,ist),ist=1,4)
          if(ilu.ne.j) then
            write(*,'(//,a,a,//)') ' Error in ',ifile
            stop
          end if
        end do
        ! read CN for impervous area                  
        read(2,'(42x,f7.0)', err=992,end=992) CNi(ict)
        ! assume impervious area CN is fixed          
        Si(ict) = 1000./CNi(ict) - 10.                
        rIai(ict) = 0.2*Si(ict)                       


        close(2)
        call openifile(1,2,ifile)           ! land use properties
        do ilu=1,nlu
          read(2,'(i6,36x,6f7.0,i7)',err=992,end=992)
     &             j,fimp(ict,ilu),
     &             fDCIA(ict,ilu), fUCIAi(ict,ilu),
     &             firg(ict,ilu)

        end do
        read(2,'(a1)') ch1      ! blank line
        read(2,'(i6)') klunv    ! native vegetation coverage type
        read(2,'(i6)') klubar   ! barren area land use type
        read(2,'(i6)') klunrf   ! land use type for no runoff/no percolation
        close(2)

        write(6,'(/,a,i6,6x,a50)') 'ICT =   ',ict,ifile
        write(6,'(a)') '  WELU                            Curve Numbers'
        write(6,'(a,a)') ' Class %Impv %DCIA %Dist %Irld',
     &                   '     A     B     C     D'

        do ilu=1,nlu
          write(6,'(i6,4f6.1,4f6.0)')
     &             ilu,fimp(ict,ilu),
     &             fDCIA(ict,ilu), fUCIAi(ict,ilu),
     &             firg(ict,ilu),
     &             (soilcn(ict,ilu,ist),ist=1,4)
          fDCIA(ict,ilu)  = fDCIA(ict,ilu) /100.
          fUCIAi(ict,ilu) = fUCIAi(ict,ilu)/100.
          fimp(ict,ilu)  = fimp(ict,ilu)/100.     ! fraction of impervious
          firg(ict,ilu)  = firg(ict,ilu)/100.     ! fraction of pervious area that is irrigated
        end do
        write(6,'(1h )')
c       write(6,'(i6,a)') kluimp,' impervious area land use type  '
        write(6,'(i6,a)') klunv ,' native vegetation coverage type'
        write(6,'(i6,a)') klubar,' barren area land use type      '
        write(6,'(i6,a)') klunrf,' no runoff land use type      '
      end do ! ict

      ! open drainage area - LU - soil data file

      call openifile(1,2,ifile)              ! *.rain
      write(6,'(/,a,a)') ifile,'is ready'


      do while (.true.)
        read(2,'(5x,a8,i5,4f12.0)',err=992,end=20)
     &            ch8,ilu,(rval(ist),ist=1,4)
        if(ch8.eq.'00000000' .or. ch8.eq.'        ') goto 20
        call ctrim8r(ch8)

        if(ilu.eq.0) ilu = klunv             ! default

        ! find HAID
        do i=1,nha
          iha=i
          if(ch8.eq.haname(i)) goto 10
        end do
        write(*,'(//,a,a,a,//)') ' HA ID ',ch8, ' is not avaiable'
        stop
 10     ict = kct(iha)
        do ist=1,4
          if(rval(ist).gt.0.) then
            rval(ist) = rval(ist) * sqm2ac
            areap = rval(ist)*(1.-fimp(ict,ilu))
            areai = rval(ist)*fimp(ict,ilu)

            dacls(iha,ilu,ist) = dacls(iha,ilu,ist)        ! LU-Soil type area
     &                       + rval(ist)
            dacid(iha,ilu,ist) = dacid(iha,ilu,ist)        ! DCIA area
     &                       + areai * fDCIA(ict,ilu)
            daciu(iha,ilu,ist) = daciu(iha,ilu,ist)        ! UCIA area
     &                       + areai * (1. - fDCIA(ict,ilu))
            dacpi(iha,ilu,ist) = dacpi(iha,ilu,ist)
     &                       + areap * firg(ict,ilu)
            dacpn(iha,ilu,ist) = dacpn(iha,ilu,ist)
     &                       + areap * (1. - firg(ict,ilu))
          end if
        end do
      end do
 20   continue
      close(2)

       if(idebug.eq.1) write(*,'(a)') ' pass 0b'

      ! assign and summarize HA-LU-ST area data
      kpt = 0
      do iha=1,nha
        ict = kct(iha)
        do ist=1,4
          cn2(iha,klubar,ist) = soilcn(ict,klubar,ist) * rmcn(iha)
          cn1(iha,klubar,ist) = cn2(iha,klubar,ist) /
     &                          (2.27 - 0.0125*cn2(iha,klubar,ist))
          cn3(iha,klubar,ist) = min(100.,
     &                              cn2(iha,klubar,ist) /
     &                              (0.44 + 0.0055*cn2(iha,klubar,ist)))
        end do
        do ilu=1,nlu
          if(ilu.eq.klunrf) goto 220
          do ist=1,4
            if(idebug.eq.1) write(*,'(a)') ' pass 0b1'
            cn2(iha,ilu,ist) = soilcn(ict,ilu,ist) * rmcn(iha)
            cn2(iha,ilu,ist) = max(0.,min(100.,cn2(iha,ilu,ist)))
            if(firg(ict,ilu).gt.0.) then
              cn1(iha,ilu,ist) = cn2(iha,ilu,ist)
              ! limit the irrigated land AMC to AMC2
            else
              cn1(iha,ilu,ist) = cn2(iha,ilu,ist)
     &                           / (2.27 - 0.0125*cn2(iha,ilu,ist))
            end if

            cn3(iha,ilu,ist) = min(100.,
     &                           cn2(iha,ilu,ist)
     &                           / (0.44 + 0.0055*cn2(iha,ilu,ist)))
            sscmax(iha,ilu,ist) = 1000. / cn1(iha,ilu,ist) - 10.
            if(idebug.eq.1) write(*,'(a)') ' pass 0b5'!

            sscmin(iha,ilu,ist) = 1000. / cn3(iha,ilu,ist) - 10.
            if(idebug.eq.1) write(*,'(a)') ' pass 0b6'

            ! initial soil storage

            haluarea(iha,ilu) = haluarea(iha,ilu) + dacls(iha,ilu,ist)

            if(idebug.eq.1) write(*,'(a)') ' pass 0b7'
            cnpi(iha,ilu,ist) = cn2(iha,ilu,ist)             
            cnpn(iha,ilu,ist) = cn1(iha,klubar,ist)          
            sspi(iha,ilu,ist) = 1000./cnpi(iha,ilu,ist) - 10.
            sspn(iha,ilu,ist) = 1000./cnpn(iha,ilu,ist) - 10.

            if(idebug.eq.1) write(*,'(a)') ' pass 0b8'
            if(dacls(iha,ilu,ist).gt.dacmin) then
              if(gmfrac(iha).gt.0.) then
                kpt = kpt+1
                iv(kpt,1) = iha
                iv(kpt,2) = ilu
                iv(kpt,3) = ist
                gmarea(iha) = gmarea(iha)                   
     &                      + dacls(iha,ilu,ist)*gmfrac(iha)
                rv(kpt,1) = dacpi(iha,ilu,ist)*gmfrac(iha)  
                rv(kpt,2) = dacpn(iha,ilu,ist)*gmfrac(iha)  
                ! this two line is modified on 1/17/2012
              end if
              if(idebug.eq.1) write(*,'(a)') ' pass 0b9'
            end if
          end do ! ist
 220      continue
        end do ! ilu
      end do ! iha

      if(idebug.eq.2) then
        open(88,file='RFdebug.out',status='UNKNOWN')
        write(88,'(a,a)') 'Debug out for ',file0
        write(88,'(a,a,a,a)')
     &           ' HA LU ST Fimp DCIA Firg',
     &           '  DACLS  DACID  DACIU  DACPI  DACPN',
     &           '  CN1  CN2  CN3 CNPI  SSPI SScMn SScMx',
     &           '  CN1  CN2  CN3 CNPN  SSPI SScMn SScMx'
        do iha=1,nha
          do ilu=1,nlu
            ict = kct(iha)
            do ist=1,4
              if(dacls(iha,ilu,ist).gt.0.)
     &          write(88,'(3i3,3f5.2,5f7.1, 2(4f5.1,3f6.1))')
     &              iha,ilu,ist,
     &              fimp(ict,ilu),fDCIA(ict,ilu),firg(ict,ilu),
     &               dacls(iha,ilu,ist),
     &               dacid(iha,ilu,ist),daciu(iha,ilu,ist),
     &               dacpi(iha,ilu,ist),dacpn(iha,ilu,ist),
     &              cn1(iha,ilu,ist),cn2(iha,ilu,ist),cn3(iha,ilu,ist),
     &              cnpi(iha,ilu,ist),
     &              sspi(iha,ilu,ist),
     &              sscmin(iha,ilu,ist),sscmax(iha,ilu,ist),
     &                cn1(iha,klubar,ist),cn2(iha,klubar,ist),
     &                cn3(iha,klubar,ist),cnpn(iha,ilu,ist),
     &                sspn(iha,ilu,ist),
     &                sscmin(iha,klubar,ist),sscmax(iha,klubar,ist)
            end do
          end do
        end do
        write(88,'(/,a,a,a,a,a,/)') 'YYYYMMDD HA LU ST',
     &     ' Precip    ETo roi_in Iai_in',
     &     '  Stor3   CNpi  prha3 ro3_in perci3   ETa3  Dprc3',
     &     '  Stor4   CNpn  prha4 ro4_in perci4   ETa4  Dpec4'

c    &     '  prha3   CNpi ro3_in Io3_in perci3 ro3_AI Io3_AI prc3AI',
c    &     '  prha4   CNpn ro4_in Io4_in perci4 ro4_AI Io4_AI prc4AI'
        write(88,'(1h )')
      end if

      if(idebug.eq.1) write(*,'(a)') ' pass 0c'


      if(kpt.gt.MXPG) then
        write(*,'(a,i7)') ' Make MXPG > ',kpt
        stop
      end if

      write(6,'(/,a,i5)') 'Area data check'
      write(6,'(a,25i10)') '   ID HAID    ',(i,i=1,nlu)
      do iha=1,nha
        write(6,'(i5,1x,a8,25f10.0)') iha,haname(iha),
     &                              (haluarea(iha,ilu),ilu=1,nlu)
        sum1 = 0.
        do ilu=1,nlu
          do ist=1,4
c           sum1 = sum1 + dac(iha,ilu,ist)
            sum1 = sum1 + dacls(iha,ilu,ist)
          end do
        end do
        pdiff = (sum1 - haarea(iha))/haarea(iha)
        if (pdiff.gt.0.02)
     &     write(6,'(i5,1x,a8,2f8.0,f8.3)')
     &          iha,haname(iha),haarea(iha),sum1,pdiff
      end do
      write(6,'(1h )')

      call openhdat(1,3,nsy,ismonth,ifile,fmtrn)
      call openhdat(1,4,nsy,ismonth,ifile,fmtev)
      call openhdat(1,5,nsy,ismonth,ifile,fmtet)

      call openmfile(1,21,nha,ofile,'Daily runoff (cfs)',file0)          ! *.FLW

      call openmfile(1,31,nha,ofile,'Monthly Precipitation, (inches)',   ! *.MPR
     &                           file0)
!     call openmfile(1,32,nha,ofile,'Monthly Mun App Water, (inches)',   ! *.MWM
!    &                           file0)
!     call openmfile(1,33,nha,ofile,'Monthly Agg App Water, (inches)',   ! *.MWA
!    &                           file0)
      call openmfile(1,34,nha,ofile,'Monthly Runoff, (inces)',           ! *.MRO
     &                           file0)
      call openmfile(1,35,nha,ofile,'Monthly Infiltration, (inches)',    ! *.MIF
     &                           file0)
!     call openmfile(1,36,nha,ofile,'Monthly Percolation, (inches)',     ! *.MDP
!    &                           file0)
      call openmfile(1,37,nha,ofile,'Monthly Evaporation,  (inches)',    ! *.MEV
     &                           file0)
      call openmfile(1,38,nha,ofile,'Monthly ET, (inches)',              ! *.MET
     &                           file0)

      call openofile(1,11,ofile,'DACpi & DACpn(iha,ilu,ist)',file0)
      write(11,'(i8)') nha
      write(11,'(25a8)') (haname(iha),iha=1,nha)

      write(11,'(i8)') kpt
      write(11,'(50i4)') (iv(kp,1),kp=1,kpt)                ! iha
      write(11,'(50i4)') (iv(kp,2),kp=1,kpt)                ! ilu
      write(11,'(50i4)') (iv(kp,3),kp=1,kpt)                ! ist
      write(11,'(25E10.4)') (rv(kp,1),rv(kp,2), kp=1,kpt)
         ! (dacpi(iha,ilu,ist),dacpn(iha,ilu,ist)           ! dacpn
      close(11)

      write(6,'(a,a)') ofile,'is done'

      call openofile(1,12,ofile,'Daily infil. by HA-LU-ST-IN',file0)  ! *.PFL
      write(6,'(a,a)') ofile,'is ready'


      if(qopt) then
        do ich=1,2
          call openifile(1,2,ifile)                ! runoff-TDS or TIN tables
          do ilu=1,nlu
            read (2,'(2i8)') i,j
            if(i.ne.ilu) stop
            ntconc(ilu,ich) = j
            do i=1,j
              read(2,'(2f8.0)') tconc(ilu,ich,i,1),tconc(ilu,ich,i,2)
            end do
          end do
          close(2)
        end do

        call openofile(1,22,ofile,'       1  daily values','') ! runoff-TDS output file
        write(6,'(a,a)') ofile,'is ready'
        call openofile(1,23,ofile,'       1  daily values','') ! runoff-TIN output file
        write(6,'(a,a)') ofile,'is ready'
      end if

!
!     start daily runoff computation with Ia = Iamax
!

      do 90 ipy=1,npor !===== Year Loop ================================

      icy = nsy + (ipy - 1)

      write(*,'(/,i5)') icy


      do 80 iwm=1,12   !===== Month Loop ===============================

      icm = iwm + (ismonth - 1)
      if(icm.gt.12) then
        icm = icm - 12
        if (icm.eq.1) icy = icy+1
      end if


      k4 = mndays(icy,icm)

      do iha=1,nha
        xmpr(iha)=0.
        xmro(iha)=0.
        xmev(iha)=0.
        xmet(iha)=0.
        xmpci(iha)=0.
      end do

      do 70 idy=1,k4   !===== Day Loop =================================

      if(nprwts.lt.nha) then
        read(3,fmtrn,end=990) m10,m20,m30,(prst(j),j=1,nprwts)
      else
        read(3,fmtrn,end=990) m10,m20,m30,(prha(j),j=1,nprwts)
      end if
      if(m10.eq.0.) goto 990

      read(4,fmtev,end=990) m11,m21,m31,(evst(j),j=1,nevwts)
      if(m11.eq.0.) goto 990
      read(5,fmtet,end=990) m12,m22,m32,(etst(j),j=1,netwts)
      if(m12.eq.0.) goto 990

      if(m10.ne.icy .or. m11.ne.icy .or.
     &   m20.ne.icm .or. m21.ne.icm .or.
     &   m30.ne.idy .or. m31.ne.idy ) then
        write(*,'(//,a,a)') ' Dates in precipitation and evaporation',
     &                      ' data files do not match'
        write(*,'(a,3i5)') ' icy,icm,idy   = ',icy,icm,idy
        write(*,'(a,3i5)') ' precipitation = ',m10,m20,m30
        write(*,'(a,3i5)') ' evaporation   = ',m11,m21,m31
        write(*,'(a,3i5)') ' ETo           = ',m12,m22,m32
        stop
      end if
      if(idebug.eq.1) write(*,'(a)') ' Pass day loop  '

      kp = 0
      do 60 iha=1,nha  !===== Drainage Area Loop =======================

      ict = kct(iha)
      roha(iha) = 0.
      cha(iha,1) = 0.
      cha(iha,2) = 0.

      if(nprwts.lt.nha) then
        prha(iha)=0.0
        sumwts=0.0
        do j=1,nprwts
          sumwts = sumwts+prwt(iha,j)
          prha(iha) = prha(iha)+prwt(iha,j)*prst(j)
        end do
      end if

      if(prha(iha).gt.rmrain(iha,2))
     &      prha(iha) = rmrain(iha,2)
     &                  + (prha(iha) - rmrain(iha,2)) * rmrain(iha,1)
c v5  prha(iha) = prha(iha) * rmrain(iha)
ccc   if(sumwts.gt.0.0) prha(iha) = prha(iha)/sumwts * rmrain(iha)

      xmpr(iha) = xmpr(iha) + prha(iha)    ! 9/9/08 moved out of ilu loop

       if(idebug.eq.1) write(*,'(a)') ' Pass HA loop 1'

      evha(iha)=0.
      sumwts=0.0
      do j=1, nevwts
        sumwts = sumwts+evwt(iha,j)
        evha(iha) = evha(iha)+evwt(iha,j)*evst(j)
      end do
      if(sumwts.gt.0.0) evha(iha) = evha(iha)/sumwts
      if(idebug.eq.1) write(*,'(a)') ' Pass HA loop 2'

      etha(iha)=0.
      sumwts=0.0
      do j=1, netwts
        sumwts = sumwts+etwt(iha,j)
        etha(iha) = etha(iha)+etwt(iha,j)*etst(j)
      end do
      if(sumwts.gt.0.0) etha(iha) = etha(iha)/sumwts
      if(idebug.eq.1) write(*,'(a)') ' Pass HA loop 3'

      ! calculate impervious area runoff                               
      roi_in = 0.                                                      
      if(prha(iha).gt.rIai(ict))                                       
     &  roi_in = (prha(iha) - rIai(ict))**2 / (prha(iha) + 0.8*Si(ict))
      evi_in =     prha(iha) - roi_in                                  
      if(idebug.eq.1) write(*,'(a)') ' Pass HA loop 4'

      do 50 ilu=1,nlu  !===== Land Use Class Loop ======================

      ro_af  = 0.
      pci_af = 0.
      pcd_af = 0.
      ev_af  = 0.
      et_af  = 0.

      ro1_af = 0.       ! runoff from DCIA                       
      ro3_af = 0.       ! runoff from pervious/irrigated area    
      ro4_af = 0.       ! runoff from pervious/non-irrigated area
      ev1_af = 0.       ! ev from DCIA
      ev2_af = 0.       ! ev from UCIA

      et3_af = 0.       ! et from pervious/irrigated area 
      et4_af = 0.       ! et from pervious/non-irrigated area
      pci3_af = 0.      ! infiltration in pervious/irrigated area    
      pci4_af = 0.      ! infiltration in pervious/non-irrigated area

      if(ilu.eq.klunrf) goto 50

      do 40 ist=1,4    !===== Soil Type Loop ===========================

        if(dacls(iha,ilu,ist).le.dacmin) goto 40
        if (gmfrac(iha).gt.0.) kp = kp + 1
        if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass ST loop 1'

      ! runoff from DCIA in (acre-ft) from (iha,ilu,ist) area
      ro1_ai = dacid(iha,ilu,ist) * roi_in
      ev1_ai = dacid(iha,ilu,ist) * evi_in

      ro1_af = ro1_af + ro1_ai / 12.
      ev1_af = ev1_af + ev1_ai / 12.
      if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass ST loop 2'

      ! runoff from UCIA in (acre-in)
      ro2_ai = daciu(iha,ilu,ist) * roi_in
      ev2_ai = daciu(iha,ilu,ist) * evi_in
      ev2_af= ev2_af + ev2_ai / 12.
      ! ro from UCIA goes to pervious area, so no need to calculate ro2_af

      if(idebug.eq.1) write(*,'(a)') ' Pass LU ST loop 1'

      ! rerouted runoff from UCIA to pervious irrigated area
      ro23 = ro2_ai * fUCIAi(ict,ilu)
      prha3 = prha(iha)
      if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass loop 4a'
      if(dacpi(iha,ilu,ist) .le. 0.) goto 35
        prha3 = prha3 + ro23/dacpi(iha,ilu,ist)

        CNpi(iha,ilu,ist) = CN1(iha,ilu,ist)
     &              + (CN3(iha,ilu,ist) - CN1(iha,ilu,ist))
     &              * (SSpi(iha,ilu,ist) - SScmin(iha,ilu,ist))
     &              / (SScmax(iha,ilu,ist) - SScmin(iha,ilu,ist))

        if(idebug.eq.1) write(*,'(a,3i4,2f5.1,2e12.4,f5.1)')
     &                       ' Pass 1b',iha,ilu,ist,
     &                       CN1(iha,ilu,ist),CN3(iha,ilu,ist),
     &                       SSpi(iha,ilu,ist), SScmax(iha,ilu,ist),
     &                       CNpi(iha,ilu,ist)
        WSpi    = 1000./CNpi(iha,ilu,ist) - 10.
        WSpi2   = 0.2 * WSpi
        if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass loop 4b'
        ro3_in = 0.
        if(prha3 .gt. WSpi2)
     &     ro3_in = (prha3 - WSpi2)**2 / (prha3 + 0.8 * WSpi)

         perci3 = max(0., prha3 - ro3_in)                ! rf_v7a reverse modification
!        perci3 = max(0., prha3 - ro3_in - WSpi2)        ! rf_v7 modification

        ! which may be OK in agricultural and urban irrigated area.

        ! now soil zone
        SSpio = SSpi(iha,ilu,ist)  ! remember initial value
        SSn = SSpio + perci3
        eta3 = etha(iha) * (SSn - SScmin(iha,ilu,ist))
     &                   / (SScmax(iha,ilu,ist) - SScmin(iha,ilu,ist))
        eta3 = min(eta3,etha(iha))

        dprc3 = max(0., SSn - SScmax(iha,ilu,ist) - eta3)
        SSpi(iha,ilu,ist) = SSpio + perci3 - dprc3 - eta3

        if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass loop 4c'
        cnst = dacpi(iha,ilu,ist) / 12.   
        ro3_af = ro3_af + ro3_in * cnst   
        pci3_af = pci3_af + perci3 * cnst 
c       ev3_af = ev3_af + evap3 * cnst    
        et3_af = et3_af + eta3 * cnst     
        rv(kp,1) = prha3                  
        rv(kp,3) = perci3                 
        if(idebug.eq.1)write(*,'(3i5,a)') iha,ilu,ist,' pass loop 4d'
 35   continue

        if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass loop 5'

      ! pervious non-irrigated area

      ! rerouted runoff from UCIA to pervious non-irrigated area
      ro24 = ro2_ai - ro23
      prha4 = prha(iha)
        if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass loop 6'
      if(dacpn(iha,ilu,ist) .le. 0.) goto 38
      if(idebug.eq.1) write(*,'(a,3i4,5e12.4)') ' Pass 2a',iha,ilu,ist,
     &                                           SScmax(iha,klubar,ist)
      if(idebug.eq.1) write(*,'(3i5,a)') iha,ilu,ist,' pass ST loop 6a'
        prha4 = prha4 + ro24/dacpn(iha,ilu,ist)
        if(idebug.eq.1) write(*,'(a,3i4,2f5.1,2e12.4,f5.1)') ' Pass 6a',
     &                       iha,ilu,ist,
     &                       CN1(iha,klubar,ist),CN3(iha,klubar,ist),
     &                       SSpn(iha,ilu,ist), SScmax(iha,klubar,ist),
     &                       CNpn(iha,ilu,ist)

        CNpn(iha,ilu,ist) = CN1(iha,klubar,ist)
     &            + (CN3(iha,klubar,ist) - CN1(iha,klubar,ist))
     &            * (SSpn(iha,ilu,ist) - SScmin(iha,klubar,ist))
     &            / (SScmax(iha,klubar,ist) - SScmin(iha,klubar,ist))
        WSpn  = 1000./CNpn(iha,ilu,ist) - 10.
        WSpn2 = 0.2 * WSpn
        ro4_in = 0.
        if(prha4 .gt. WSpn2)
     &     ro4_in = (prha4 - WSpn2)**2 / (prha4 + 0.8 * WSpn)

        perci4 = max(0., prha4 - ro4_in)                 ! rf_v7a reverse modification
!       perci4 = max(0., prha4 - ro4_in - WSpn2)         ! rf_v7 modification
        SSn = SSpn(iha,ilu,ist) + perci4
        eta4 = etha(iha) * (SSn - SScmin(iha,klubar,ist))
     &            / (SScmax(iha,klubar,ist) - SScmin(iha,klubar,ist))
        eta4 = min(eta4,etha(iha))

        dprc4 = max(0., SSn - SScmax(iha,klubar,ist) - eta4 )
        SSpno = SSpn(iha,ilu,ist)
        SSpn(iha,ilu,ist) = SSpn(iha,ilu,ist) + perci4 - dprc4 - eta4

        cnst = dacpn(iha,ilu,ist) / 12.
        ro4_af = ro4_af + ro4_in * cnst
        pci4_af = pci4_af + perci4 * cnst

        et4_af = et4_af + eta4 * cnst
        rv(kp,2) = prha4
        rv(kp,4) = perci4

        if(idebug.eq.1)write(*,'(3i5,a)') iha,ilu,ist,' pass ST loop 6d'
 38     continue
        if(idebug.eq.1)write(*,'(3i5,a)') iha,ilu,ist,' pass ST loop 7'


c       if(idebug.eq.2 .and. prha(iha).gt.0.0) then
        if(idebug.eq.2                       ) then
c         vv1 = ro3_in * dacpi(iha,ilu,ist)
c         vv3 = perci3 * dacpi(iha,ilu,ist)
c         vv4 = ro4_in * dacpn(iha,ilu,ist)
c         vv6 = perci4 * dacpn(iha,ilu,ist)

          write(88,'(i4,2i2.2,3i3,4f7.4,
     &          2(f7.4,f7.1,5f7.4))')
     &            icy,icm,idy, iha,ilu,ist,
     &            prha(iha),etha(iha),roi_in,evi_in,                    ! inches
     &    SSpio, CNPi(iha,ilu,ist),prha3,ro3_in,perci3, ETa3,dprc3,
     &    SSpno, CNPn(iha,ilu,ist),prha4,ro4_in,perci4, ETa4,dprc4
        end if

 40   continue  ! ist

      if(idebug.eq.1) write(*,'(a)') ' Pass LU ST loop 3a'
      ro_af = ro1_af + ro3_af + ro4_af
      pci_af = pci3_af + pci4_af
c     ev_af = ev1_af + ev2_af + ev3_af + ev4_af   
cccc  pr_af = prha(iha) * dacls(iha,ilu,ist) / 12.
      et_af = et3_af + et4_af                     
                                                  
      roha(iha) = roha(iha) + ro_af

cccc  xmpr(iha) = xmpr(iha) + prha(iha)    ! 9/9/08 move out of ilu loop
      xmro(iha) = xmro(iha) + ro_af
      xmpci(iha) = xmpci(iha) + pci_af    ! infiltration
      xmev(iha) = xmev(iha) + ev_af
      xmet(iha) = xmet(iha) + et_af

c     if(icy.eq.1994 .and. icm.eq.10 .and. idy.eq.7)
c    &  write(*,'(2i5,7f9.4)')
c    &      iha,ilu,ro1_af,ro2_af,ro3_af,ro4_af,pci3_af,pci4_af

      if(idebug.eq.1) write(*,'(a)') ' Pass LU ST loop 4'

      if(qopt .and. ro_af.gt.0.) then
        ro_inch = ro_af*12./haluarea(iha,ilu)
        conc = 0.                                                 
        if(idebug.eq.1)write(*,'(2i4,f10.1,2e14.4,a)') iha,ilu,
     &       haluarea(iha,ilu),ro_af,ro_inch, ' Pass LU ST loop 4'
        do ich=1,2                                                
          nt = ntconc(ilu,ich)                                    
            do i=1,nt                                             
              if(ro_inch .le. tconc(ilu,ich,i,1)) then            
                conc = tconc(ilu,ich,i-1,2)                       
     &               + (tconc(ilu,ich,i,2) - tconc(ilu,ich,i-1,2))
     &               / (tconc(ilu,ich,i,1) - tconc(ilu,ich,i-1,1))
     &               * (ro_inch            - tconc(ilu,ich,i-1,1))
                goto 30                                           
              end if                                              
            end do                                                
            write(*,'(a)') ' Error in conc interpolation'         
            write(*,'(3i5,e12.4)') iha,ilu,ich,ro_inch
            stop                                                  
 30       continue                                                
          cha(iha,ich) = cha(iha,ich) + ro_af*conc                
        end do !ich
      end if !qopt

      if(idebug.eq.1) write(*,'(a)') ' Pass LU ST loop 5'
 50   continue   ! ilu               !----------------------------------


      if(qopt) then                                                 
        if(roha(iha).gt.0.) then                                    
          do ich=1,2                                                
            if(cha(iha,ich).eq.0.) then                             
              write(*,'(a)') ' Error in conc calc 1'                
              write(*,'(2i5,2e12.4)') iha,ich,roha(iha),cha(iha,ich)
              stop                                                  
            else                                                    
              cha(iha,ich) = cha(iha,ich) / roha(iha)               
            end if                                                  
          end do                                                    
        else                                                        
          do ich=1,2                                                
            if(cha(iha,ich).gt.0.) then                             
              write(*,'(a)') ' Error in conc calc 2'                
              write(*,'(2i5,2e12.4)') iha,ich,roha(iha),cha(iha,ich)
              stop                                                  
            end if                                                  
          end do                                                    
        end if                                                      
      end if !qopt                                                  

      if(idebug.eq.1) write(*,'(a)') ' Pass HA loop 5'

      ! convert af/day to cfs
      roha(iha) = roha(iha) * af2cfs

 60   continue         !----- HA Loop ----------------------------------
!ccccc      write(6,'(i4,2i2.2,i7)') icy,icm,idy,kp !delete

      write(12,'(i4,2i2.2,/,(25f8.5))') icy,icm,idy,
     &                                  ((rv(kp,i),i=1,4),kp=1,kpt)
c v5 &                                  (rv(kp),kp=1,kpt)


      write(21,'(i4,2i2.2,25f8.2,/,(8x,25f8.2))') icy,icm,idy,
     &                                            (roha(iha),iha=1,nha)
      if(qopt) then 
        write(22,'(i4,2i2.2,25f8.2,/,(8x,25f8.2))') icy,icm,idy,
     &                                          (cha(iha,1),iha=1,nha)
        write(23,'(i4,2i2.2,25f8.2,/,(8x,25f8.2))') icy,icm,idy,
     &                                          (cha(iha,2),iha=1,nha)

      end if !qopt 



 70   continue         !----- Day Loop ---------------------------------


      if(idebug.eq.1) write(*,'(a)') ' Pass out of day loop'

      ! output data for water budget on monthly or larger time steps
      do iha=1,nha
        cnst = haarea(iha)/12.
        xmro(iha) = xmro(iha)/cnst
        xmpci(iha) = xmpci(iha)/cnst
        xmev(iha) = xmev(iha)/cnst
        xmet(iha) = xmet(iha)/cnst  
      end do

      write(31,8081) icy,icm,(xmpr(iha),iha=1,nha)
      write(34,8081) icy,icm,(xmro(iha),iha=1,nha)
      write(35,8081) icy,icm,(xmpci(iha),iha=1,nha)
      write(37,8081) icy,icm,(xmev(iha),iha=1,nha)
      write(38,8081) icy,icm,(xmet(iha),iha=1,nha)

!8080 format(i4,2i2.2,25f8.1,/,(8x,25f8.1))
 8081 format(i4,i2.2,25f8.3,/,(6x,25f8.3))

 80   continue         !----- Month Loop -------------------------------

 90   continue         !----- Year Loop --------------------------------

 100  continue

 990  write(*,'(//,a,/)') ' End of Execution'
 991  close(3)
      close(4)
      close(6)
      close(7)
      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(21)
      if(qopt) then !vg1
        close(22) !vg1
        close(23) !vg1
      end if !vg1
      stop

 992  write(*,'(//,a,a)') ' Data error in ',ifile
      stop
      end


