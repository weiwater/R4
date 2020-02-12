! modified from qrz_v3s7a to run to use daily ET data on HSA, on 12/2/79
! this modification is necessary to use climate change information.
! modified from qrz_v3s7 to run witout quality operion

! qrz_v3s7   modified from qrz_v3s6, RW monthly supply left-over transfer is removed

! qrz_v3s6.for is modified from qrz_v3s5.for, for
!   the situation that, 1 HSA can receive RW from RW supplier (it can be WSA), 
!   and the suppler can have multiple RW pressure zone, which can be supplied by different RPs

c Modified 10/7/2018, RW supply to zone 
c    Frwz(irwz,1)   monthly supply to urban/park area (af/month)
c    Frwz(irwz,2)   monthly supply to ag use (af/month)
c     
c  qrz_v3s3.for     WEI, JH, October 2018

c qrz_v3s1.for is modified from rz_v3s18.for to include water quality simulation.
c
c rz_v3s18.FOR is modified from rz_v3s3.for to increase MXLU to 20  11/15/2014
c
c 12/14/15 rz_v3S1.for modification
c   use multiple ET data (net) and ET weights
c 10/20/15 rz_v3s.for --- cleaned version of RZ_V3R4B
c 
c RZ_V3r3b.for to RZ_V3r4b.for   for dimension
c   modified from RZ_v4r3a.for
c   user can specified irrigation/nonirrigation month for irrigated land use
c   1/3/2011
c RZ_V3r3a.FOR
c   modified from RZ_V3r3.for
c   to allow multiple zone for irrigation efficiency and crop coefficient
c   for example, older communities have larger lots and different plants and
c   older irrigation method.
c   12/28/2011
c
c RZ_V3r3.for

c RZ_V3r2.for
c   modified from rz_v3r1.for on 7/26/2011
c     to add nK/bug = 2, for daily water budget output
c     to remove da_chrc.prn file
c   need to add an external file to read gwbr(iha)
c
c

c Rootzone soil moisture accounting program
c Modified from Q:\Chino\perc\for\RZ\RZ_V7.for
c2/22/07
c This program reads daily data from
c   runoff model, version rohhq2v_.for, output for daily percolation,
c   daily evapotranspiration data
c and calculates deep percolation out of root zone area.
c
c Rz_v2r2 is modified from rz_v2r1 August 2007
c updates are:
c   1 debug features for rootzone soil moisture accounting check
c   2 critial soil moisture and irrigation timing is separated
c
c RZ_v2r3 is modified from rz_v2r2`
c 8/28/07 for
c to allow depressed irrigation for vineyard
c to allow dairy waste irrigation on land use type 2
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
c 33               monthly TDS deep perc
c 53               monthly N deep perc
c 14  output *.mwh monthly M&I applied water  by HA for water budget
c 34               monthly average TDS of M&I
c 54               monthly average N 
c 15  output *.awh monthly ag applied water   by HA for water budget
c 35               monthly average TDS of AG supply
c 55               monthly average N 
c 16  output *.rwh monthly recycled water by HA for water budget
c 36               monthly average TDS of recycled water
c 56               monthly average N 
c 17               monthly dairy washwater application by HA
c 18               reserved
c 19               reserved
c 20  output *.sto storage at end 
c 21  output *.mws monthly applied water by WSA
cremove 22  output *.aws monthly applied water by AG water source - no need
cremove 23  output *.rws monthly recycled water by HSA
cremove 24  output *.dws monthly dairy washwater by source
c 25  output *.law monthly total AW by lu type ????????????????????????????????
c 26  output *.smc SM in HSA at month end
c 27  output *.smt TDS in soil at month end
c 28  output *.smn N in soil at month end
c 91  output *.pcl total percolation by LU type
c 92  output *.pcr total rainfall percolation by LU type
c 93  output *.pca total applied water percolation by LU type 

!s1
      include 'RF_DIM_4.MAX'

      character*1  ch1, cstar
      character*8  haname(MXHA),haid8(MXHA),ch8
      character*40 fmtin, fmtet,fmym80,fmym81,fmym82,fma6a8,fm6x8i
      character*50 file0,ifile,ofile,bl50,clu(21)

      dimension    haarea(MXHA),           ! area of HA on GW modeling area
     &             xminf(MXHA),            ! monthly total infiltration (af)
     &             xmet(MXHA),             ! monthly evapotranspiration (af)
     &             xmaw(MXHA,4),           ! monthly applied water (af) 1 MW 2 AW 3 RW 4 dairy ww+AW
     &             xmpc(MXHA),             ! monthly percolation (af)
     &             xmpcr(MXHA),            ! deep perc of rainfall 
     &             xmpca(MXHA),            ! deep perc of applied water
     &             xmsto(MXHA)

      dimension    xmmws(MXWS)             ! monthly water supply by WSA
cremove     &             xmrws(MXRW) 
     
      dimension    cTxmpc(MXHA), cNxmpc(MXHA)
      dimension    cTxmaw(MXHA,3), cNxmaw(MXHA,3)
     
      dimension    xmpcl(MXLU),
     &             xmawl(MXLU),            ! monthly applied water by land use
     &             xmpclr(MXLU),
     &             xmpcla(MXLU)

      dimension    gwbr(MXHA),             ! fraction of HSA, on top of groundwater basin
     &             etwt(MXHA,MXCT),
     &             dget(MXHA)               ! daily gage ET data, modification for qrz_v3s7b.for
c     &             dget(MXCT)              ! daily gage ET data
     
      logical*1 apwrq(MXLU),debug,qopt     ! AW required for the landuse, debug, water quality option
      logical*1 irmn(MXCT,MXLU,12)         ! irrigation month
      logical*1 ihach(MXHA)
      dimension apwef(MXCT,MXLU),          ! irrigation efficiency
     &          rzd(MXCT,MXLU),            !
     &          cropcoef(MXCT,MXLU,12),    !
     &          isrc(MXCT,MXLU),           ! source of applied water (1 through 4)
     &          apad(MXCT,MXLU,3),         ! average percentage of allowable deficit
     &          csm(MXCT,MXST,MXLU,3)      ! critical soil moisture
                                                
      dimension ihaiz(MXHA), ihapt(MXHA)
                                           ! ihaiz() Rootzone definition zone
                                           ! ihapt(ii) ii=1,nhapr, ihapt(ii) = iha number to print

      dimension fc(MXST), wp(MXST)   !!!!!, SMC(MXST), AWC(MXST)

      ! variables for HSA-LU-ST overlay polygons
      dimension ihls(MXPG,3),                   ! 1 HSA number, 2 land use type, 3 soil type for polygon kp
     &          dac(MXPG,2),                    ! 1 irrigated pervious area, 2 non-irrigated pervious area
     &          rinf(MXPG,4),                   ! 1 ER on dac(,1), 2 ER on dac(,2), ! 3 and 4 for RI
     &          rsm(MXPG,2),                    ! 1 SM of dac(,1), 2 SM of dac(,2)
     &          rfc(MXPG,2),                    ! 1 FC of dac(,1), 2 FC of dac(,2), rootzone depth * fc
     &          rwp(MXPG,2),                    ! 1 WP of dac(,1), 2 WP of dac(,2), rootzone dpeth * wp
     &          rapw(MXPG),                     !   AW with irrigation efficiency on dac(,1)
     &          rdp(MXPG,2),
     &          rdpr(MXPG,2),                   ! 1 DP of rainfall on dac(,1), 2 DP of rainfall on dac(,2)
     &          reta(MXPG,2),                   ! 1 actual ET from dac(,1), 2 actual ET from dac(,2)
     &          rdpa(MXPG,2)                    ! 1 DP of applied water to dac(,1), 2 is reserved 
      dimension rcsm(MXPG,3,2)                  ! (,1) allowable SMC of dac(,1), 
                                                ! (,2) start irrig.,  
                                                ! (,3) end irrig.
                                                
                     ! ER effective rainfall, RI infall infiltration, AW applied water, DP deep perc
                     ! FC field capacity WP wilting point, SM Soil moisture, SMC critical soil moisture
      
                            
!s2     
      real      fhawsa(MXHA,MXWS),rv(50) ! fraction of HSA iha, which belongs to WSA iws
      integer   krwhsa(MXHA),            ! RW zone number for HSA iha
     &          kaghsa(MXHA)             ! HSA ID number in ag water supplied by local GW 
     
      real      atmdep(2,2),             ! atm deposit, i=1 TDS, 2 N, j=1 dry (lb/acre/day), 2 wet (mg/L)
     &          dmanu(2,MXLU),           ! dairy manure, i=1,TDS, 2 N, l = land use (lb/acre/day)
     &          ftapp(2,MXLU),           ! fertilizer,   i=1,TDS, 2 N, l = land use lb/acre/day)
     &          ftappC(2,MXLU,12),       ! urban fertilizer application (lb/acre/day)
     &          rNuptake(MXLU,12),       ! nitrogen uptake by LU and month (lb/acre/day)
     &          tnar(2,2,2)              ! historical TDS N application by manure and fertilizer to ag area
     
      real      sNloss(2)                ! fraction of daily loss from soil surface (1) to air, (2))to RZ  
                                         ! applie to AND !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!and dairy washwater
      
      real      cTmws(MXWS),             ! MW supply TDS, update monthly
     &          cNmws(MXWS),             !    ""     N       ""
     &          cTaws(MXHA),             ! local AG water supply TDS,      ""
     &          cNaws(MXHA),             !    ""                 N,        ""
     &           Frwz(MXRW,2),           ! RW supply to RW Zone (af/month), 1 for Urban/park, 2 for ag
     &           Qrwz(MXRW,2),           ! left-over transfer to following month
     &          cTrp(MXRW),
     &          cNrp(MXRW),
     *          fharp(MXHA,MXRW)
!v6     &          cTrwz(0:MXRW),           ! RW supply TDS,      ""
!v6     &          cNrwz(0:MXRW)            !    ""     N,        ""
                                         !&          rw2ag(MXHA)              ! RW supply (in/day) to AG in HSA(0 for no RW, 999 for unlimit)
                ! RW supply to HSA, and quality, updated monthly
cdelete      integer   krwt(MXLU)                
                
      real      cTmwh(MXHA), cNmwh(MXHA), 
     &          cTrwh(MXHA), cNrwh(MXHA), 
     &          cTawh(MXHA), cNawh(MXHA)
                
cccccccccccc      real       Qmwh(MXHA),                 ! monthly MW supply to HSAs
cccccccccccc     &           Qrwh(MXHA),                 ! monthly RW supply to HSAs
cccccccccccc     &           Qawh(MXHA),                 ! monthly GW supply to HSAs (Ag water)
cccccccccccc     &           Qdwh(MXHA)                  ! monthly Dairy Wash W supply to HSAs
                
      real sTmass(MXPG,2), sNmass(MXPG,2),   ! mass on surface of polygon k
     &     rTMass(MXPG,2), rNmass(MXPG,2),   ! mass in rootzone of polygon K
     &     rTleach(MXPG,2), rNleach(MXPG,2)  ! mass in deep perc boloe root zone
     
      real TNtr(MXHA,4,2)    ! 1 trasnfer from surface, 2 applied water, manure, fertilizer,  3 uptake,   4 deep perc
                                ! 1 TDS, 2 N      

cdelete      data krwt/0,2,2,2,0, 0,1,1,0,0, 1,1,0,0,1, 0,1,1,1,2/                   !!!!!!!!!!, 1,1,1,1,0/
     
!     ------------------------------------------------------------------
      cstar = '*'
      sqm2ac = 640./1609.**2                 ! conversion factor from sq m to acres
      af2cfs = 43560./86400.                 !                   from acre-ft/day to cfs
      cfact1 = 43560.*28.32/(12.*453600.)    ! conv factor from (in/day*mg/L to lb/ac-day)
      cfact2 = 453600. / (43560. * 28.32)    ! convert lbs/acre-ft to mg/L
      bl50 = '                                                  '
      ndwf = 1
      fmym80 = '(i4,i2.2,500f8.0)'
      fmym81 = '(i4,i2.2,500f8.1)'
      fmym82 = '(i4,i2.2,500f8.2)'
      fma6a8 = '(a6,500a8)'
      fm6x8i = '(6x,500i8)'
      ! read input data

      call getcon(file0,n1)

      open(1,file=file0,status='OLD')
      call skipline(1,'*')

      read(1,'(i5)') isyr,ismn,npor,nha,nlu,nst
        
                   ! isyr     starting year for simulaiton (starts Oct. 1)
                   ! npor    number of years to simulate
      isdy = 1
      qopt = .FALSE.
      read(1,'(2i5)') ndebug, iv
      if(iv.gt.0) qopt=.TRUE.
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
     &    'yyyym DD IHA ILU IST   K',
     &    ' Infilt Apl-W.   E-ET    D-P     SM  DP-AI',
     &    ' Infilt   E-ET    D-P     SM  DP-AI'
      end if

!s3      
      ! main output file
      call OpnOF2(1,6,ofile,'RZ Model main output file',file0)
      write(6,'(i5)') isyr,ismn,npor,nha,nlu,nst,ndebug,nz,net

      ! read zone number for HSA
      ihach = .FALSE.
      call OpnIF2(1,2,ifile,cstar)
      nhapr = 0
      read(2,'(a)') fmtin
      do iha=1,nha
        read(2,fmtin) iii,haid8(iha),haarea(iha),
     &                ihaiz(iha),gwbr(iha), iv
        if (iii.ne.iha) then
          write(*,'(//,a,a)') ' Problem in ',ifile
          stop
        end if
        if(iv.gt.0) then
          nhapr = nhapr+1               ! number of HSA to print (excluding arlinton and upper Temescal)
          ihapt(nhapr) = iha
          ihach(iha) = .TRUE.
        end if
      end do
      close(2)
      
      etwt = 1.
c      if(net.gt.1)  then   ! modification for qrz_v3s7b
      if(net.lt.nha)  then      
        call OpnIF2(1,2,ifile,cstar)
        write(6,'(/,a)') ' ET weights for each HSA'
        write(6,'(a)') '   ID  HANAME   IZ  ET WEIGHTS' 
        read(2,'(a)') fmtet
        do iha=1,nha
          read(2,fmtet) j,ch8,(etwt(iha,i),i=1,net)
          write(6,'(i5,2x,a8,i5,10f8.2)') j,ch8, ihaiz(iha),
     &                (etwt(iha,i),i=1,net)
        end do
        close(2)
      else   ! net = nha
        read(1,'(a)') ch1          ! blank line     
      end if        

!s4      
      ! read HSA - WSA data, the last wsa should be 'default'
      call OpnIF2(1,2,ifile,cstar)
      read(2,'(i4)') nwsa
      call skipline(2,cstar)
      fmtin = '(i4,8x,50f8.0)'
      do iha=1,nha
        read(2,fmtin) iii,(fhawsa(iha,iwsa), iwsa=1,nwsa)
        if(iii.ne.iha) stop
        sum = 0.
        do iwsa=1,nwsa-1
          sum = sum + fhawsa(iha,iwsa)
        end do
        if(sum.gt.1.0) then
          do iwsa=1,nwsa-1
            fhawsa(iha,iwsa) = fhawsa(iha,iwsa) / sum
          end do
          fhawsa(iha,nwsa) = 0.
        else 
          fhawsa(iha,nwsa)  = 1. - sum
        end if
      end do
      close(2)

      
      ! WELU land use property
      do iz=1,nz
        call OpnIF2(1,2,ifile,cstar)
        read(2,'(a1)') ch1          ! format
        write(6,'(/,a,i3)') ' Root zone variables for zone',iz
        do ilu=1,nlu
          read(2,'(i6,36x,28x,4f7.0,i7)') j,
     &                 apad(iz,ilu,1),              ! allowable deficit for SMc
     &                 apad(iz,ilu,2),              ! allowable deficit for irrigation to start
     &                 apad(iz,ilu,3),              ! allowable deficit for irrigation to end
     &                 apwef(iz,ilu),               ! irrigation efficiency
     &                 isrc(iz,ilu)                 ! source 1 M&I, 2 Ag, 3 phreatophyte
          write(6,'(i6,4f7.0,i7)') j,
     &                 apad(iz,ilu,1),              ! allowable deficit for SMc                 
     &                 apad(iz,ilu,2),              ! allowable deficit for irrigation to start 
     &                 apad(iz,ilu,3),              ! allowable deficit for irrigation to end   
     &                 apwef(iz,ilu),               ! irrigation efficiency
     &                 isrc(iz,ilu)                 ! source 1 M&I, 2 Ag
        end do
!s5     
        read(2,'(a1)') ch1
        read(2,'(i6)') klunv      ! LU type for native vegetation
        read(2,'(i6)') klubar     ! LU type for barren
        read(2,'(i6)') klunrf     ! LU type for no-runoff
        read(2,'(i6)') kluimp     ! LU type for impervious
        write(6,'(a,5i5)') ' klunv, klubar, klunrf, kluimp = ',
     &                      klunv, klubar, klunrf, kluimp 
        close(2)
      end do !iz

      ! WELU rootzone property - crop coefficients
      do iz=1,nz
        call OpnIF2(1,2,ifile,cstar)
        write(6,'(/,a,i3)') 'Root Zone Parameters for zone',iz
        
        do ilu=1,nlu
          if(apwef(iz,ilu).gt.0.) then
            apwef(iz,ilu) = 100./apwef(iz,ilu)        ! convert % to ratio
            apwrq(ilu) = .true.      ! applied water required
          else
            apwrq(ilu) = .false.
            isrc(iz,ilu) = 0
          end if
        
          read(2,'(i8,a50,20f8.0)') j,clu(ilu),
     &                              rzd(iz,ilu),
     &                              (cropcoef(iz,ilu,icm),icm=1,12)
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
     &                               rzd(iz,ilu),
     &                               (cropcoef(iz,ilu,icm),icm=1,12)
          rzd(iz,ilu) = rzd(iz,ilu)*12.            ! convert (ft) to (inches)
        end do  ! ilu
        close(2)
      end do !iz

!s6
      ! soil property
      call OpnIF2(1,2,ifile,cstar)

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
      ! csm(,,,1) --- critical soil moisture
      ! csm(,,,2) --- soil moisture to start irrigation
      ! csm(,,,3) --- soil moisture to end irrigation

      ! HA-LU-ST data transferred from Runoff Model.
      call OpnIF2(1,2,ifile,cstar)                ! PAR file
      write(6,'(/,a,a)') 'File = ',ifile
      read(2,'(i8)') iha
      if(iha.ne.nha) stop
      read(2,'(25a8)') (haname(iha),iha=1,nha)
      do iha=1,nha
        call ctrim8r(haid8(iha))
        if(haid8(iha).ne.haname(iha)) stop
      end do

      read(2,'(i8)') kpg
      if(kpg.gt.MXPG) then
        write(*,'(//,a//)') ' Increase dimension parameter MXPG'
        stop
      end if
      
      read(2,'(50i4)') (ihls(kp,1),kp=1,kpg)      ! iha
      read(2,'(50i4)') (ihls(kp,2),kp=1,kpg)      ! ilu
      read(2,'(50i4)') (ihls(kp,3),kp=1,kpg)      ! ist
      read(2,'(25f10.0)') (dac(kp,1),dac(kp,2),kp=1,kpg)       ! in pervious area in (acres)
      close(2)
      
      write(*,'(a,a)') ' done with ',ifile !v2r3a delete

      do kp=1,kpg                                             !!!!delete
        iha = ihls(kp,1)                                      !!!!delete
        xmet(iha) = xmet(iha) + dac(kp,1) + dac(kp,2)         !!!!delete
      end do                                                  !!!!delete
      
!s7      
      write(6,'(a,i5)') 'NHA = ',nha
      write(6,'(a,i5)') 'KPG = ',kpg

      ! initialize arrays for temporary use
      xmpc = 0.      ! to hold pervious area data temporarily
      xmaw = 0.      ! to hold pervious M&I area data temporarily
                     ! to hold pervious AG area data temporarily

      write(6,'(/,a,a)') 'Polygons with area < 0.1 acres will be ',
     &                   'removed to eleminate GIS operation noise.'


      xmpcl = 0.

      do kp=1,kpg
        iha = ihls(kp,1)
        ilu = ihls(kp,2)
        ist = ihls(kp,3)
        iz  = ihaiz(iha)
c        write(6,'(i6,2x,a8,3i4,2f10.2)') kp,haname(iha),
c     &                        iha,ilu,ist,dac(kp,1),dac(kp,2)
        write(6,*) kp,haname(iha),
     &             iha,ilu,ist,dac(kp,1),dac(kp,2),ihach(iha)
        if(dac(kp,1).lt.0.1) then
          if(dac(kp,1).gt.0.) write(6,'(16x,a)') 'Removed 1'
          dac(kp,1) = 0.
        end if
        if(dac(kp,2).lt.0.1) then
          if(dac(kp,2).gt.0.) write(6,'(16x,a)') 'Removed 2'
          dac(kp,2) = 0.
        end if

        if(ihach(iha))xmpcl(ilu) = xmpcl(ilu) + (dac(kp,1)+dac(kp,2))        !*gwbr(iha) is applied in RF model

        if(ist.eq.0) then
          rfc(kp,1) = 0.
          rwp(kp,1) = 0.
          rcsm(kp,1,1) = 0.
          rcsm(kp,2,1) = 0.
          rcsm(kp,3,1) = 0.
          rfc(kp,1) = 0.
          rfc(kp,2) = 0.
          rwp(kp,1) = 0.
          rwp(kp,2) = 0.
          rcsm(kp,1,1) = 0.
          rcsm(kp,1,2) = 0.
        else
          ! irrigated pervious area
          rfc(kp,1) = rzd(iz,ilu) * fc(ist)
          rwp(kp,1) = rzd(iz,ilu) * wp(ist)
          rcsm(kp,1,1) = rzd(iz,ilu)*csm(iz,ist,ilu,1)            ! allowable critical SMC
          rcsm(kp,2,1) = rzd(iz,ilu)*csm(iz,ist,ilu,2)            ! SMC to start irrigation
          rcsm(kp,3,1) = rzd(iz,ilu)*csm(iz,ist,ilu,3)            ! SMC to end irrigation
          ! non-irrigated pervious area
          rfc(kp,2) = rzd(iz,klubar) * fc(ist)
          rwp(kp,2) = rzd(iz,klubar) * wp(ist)
          rcsm(kp,1,2) = rzd(iz,klubar)*csm(iz,ist,klubar,1)      ! SM on non-irrigated Perv area

          xmpc(iha) = xmpc(iha) + dac(kp,1) + dac(kp,2)           ! temporary holding for pervious area
          if(apwrq(ilu)) then
            rsm(kp,1) = (rcsm(kp,2,1) + rfc(kp,1))/2.              ! initial soil moisture  10/17/18
            xmaw(iha,isrc(iz,ilu)) = xmaw(iha,isrc(iz,ilu)) + dac(kp,1)
          else
            rsm(kp,1) = (rwp(kp,1) + rcsm(kp,1,1))/2.
          end if
          rsm(kp,2) = (rwp(kp,2) + rcsm(kp,1,2))/2.                 ! initial soil moisture 10/17/18
        end if  ! ist.gt.0
      end do ! kp

!s8      
      write(6,'(/,a)') 'Area summarized by HA'
      write(6,'(a)') '    ID  HANAME     TOTAL   XMAW_1    XMAW_2'
      do iha=1,nha
        if(haarea(iha).gt.0.) write(6,'(i6,2x,a8,10f10.0)') iha,
     &             haname(iha),haarea(iha) !!!!!!!!!!!!!!!!!!!!!!!!!!!!,(xmaw(iha,i),i=1,4)
      end do

      do iz=1,nz
        do ilu=1,nlu
          write(6,*) ilu,isrc(iz,ilu),apwrq(ilu)
        end do
      end do

      ! HA-LU-ST daily infiltration data transferred from Runoff model
      call OpnIF2(1,4,ifile,cstar)                        ! PFL file
      Call LocTS1(4,isyr,ismn,isdy)

      ! Daily evapotranspiration data
      call OpnIF2(1,3,ifile,cstar)
      read(3,'(a40)') fmtet
      ! set the file at the simulation starting date
      call LocTS1(3,isyr,ismn,isdy)                         ! assuming (i4,2i2) for year, month and day
      
      if(qopt) then ! * read water quality parameters       ! <><><><><><><><><><><><><><><><><><><><><><>
        call OpnIF2(1,2,ifile,' ')
        
        do ic = 1,2
          call skipline(2,cstar)
          read(2,'(a)') fmtin
          do jm=1,12
            read(2,fmtin) im,(ftappC(ic,il,jm), il=1,nlu)   ! fertilizer appliaction(lb/ac/day) 
          end do                                            ! ic=1 TDS, ic=2 N
        end do
        
        call skipline(2,cstar)
        read(2,'(a)') fmtin
        do jm=1,12
          read(2,fmtin) im,(rNuptake(il,jm), il=1,nlu)      ! nitrogen uptake by plants (lb N/ac/day)
        end do                                              ! 
        
        call skipline(2,cstar)
        read(2,'(a)') fmtin
        read(2,fmtin) cTdwinc, cNdwinc                      ! dairy WW quality increment over local GW
        ! use this data if historical data is not available.

        call skipline(2,cstar)
        read(2,'(a)') fmtin
        read(2,fmtin) sNloss(1)               ! fraction of TDS N in soil surface lost to volatilization        
        read(2,fmtin) sNloss(2)               ! fraction of TDS, N in soil surface lost to soil (root) zone
        read(2,fmtin) rNloss                  ! fraction of TDS, N in rootzone lost to volatilization
        ! all nitrogen and TDS loading are estimated after volatilization, so these rates are not important
        close(2)                                         ! applied to all pervious lands.
          
        !-------------------------------------------------------------------
        !
        ! Open monthly input data files
           ! 71  Monthly TDS concentration in MW supply, iwsa=1,nwsa                           
           ! 72    "     N
           ! 73    "     TDS concentration in AG supply, for all IHA with AG, iagh = 1,nagh
           ! 74    "     N
           ! 2   RW zone data
           ! 75  Monthly RW supply to RW, irwz=1,nrwz 
           ! 76    "     TDS concentration in RW, irwz=1,nrwz
           ! 77    "     N   
           ! 78    "     Dairy washwater spray quantity, and TDS and N, if available
           ! 79    "     AG area, manure, fertilzer TSD and N loading
           ! 80    "     atmospheric TDS and nitrogen deposition in dry and wet condition
 
        do iu=71,72                        ! MW TDS and N 
          call OpnIF2(1,iu,ifile,cstar)
          read(iu,'(i4)') iv
          if(iv.ne.nwsa) stop
          read(iu,'(a1)') ch1            !!call skipline(iu,'*')
          call LocTS3(iu,isyr,ismn)
        end do

        ! TDS and N in local GW for AG supply
        iu = 73
          call OpnIF2(1,iu,ifile,cstar)      
          read(iu,'(i4)') nagh                ! # of HSAs with ag water supply + default
          read(iu,fm6x8i) (kaghsa(iagh),iagh=1,nagh)
          do iagh=1,nagh                      ! just to read serial numbers
            if(kaghsa(iagh).ne.iagh) stop
          end do
          read(iu,fm6x8i) (kaghsa(iagh),iagh=1,nagh)  ! hsa ID nubmers to supply groundwater from GWB below HSA
          if(kaghsa(nagh).ne.MXHA) then
            write(*,'(a,i3,a)') ' File with IU',iu,
     &                          ' problem. Stop execution'
            stop
          end if
          read(iu,'(a1)') ch1               ! skip HSA ID line
          call LocTS3(iu,isyr,ismn)
          iu = 74
          call OpnIF2(1,iu,ifile,cstar)      
          read(iu,'(i4)') iv                ! # of HSAs with ag water supply + default
          if(iv.ne.nagh) stop
          do i=1,3        ! assume file structurea are same for TDS and N
            read(iu,'(a1)') ch1
          end do
          call LocTS3(iu,isyr,ismn)
                      
        ! RW suply for irrigation and AG use
        
        ! RW supply zone data
        call OpnIF2(1,2,ifile,cstar)        ! RW zone - HSA data      
        read(2,'(i4)') nrwz
        read(2,'(a)') fmtin
        krwhsa = 0
        do iha=1,nha
          read(2,fmtin) kha,(rv(i),i=1,nrwz+1)
          if(kha.ne.iha) stop
          if(rv(nrwz+1).gt.0.) then
            rvmax = 0.
            do i=1,nrwz
              if(rv(i).gt.rvmax) then
                krwhsa(iha) = i
                rvmax = rv(i)
              end if
            end do
          end if     
          write(6,'(2i5,20f8.2)') iha,krwhsa(iha),(rv(i),i=1,nrwz+1)            
        end do
        close(2)
        
        ! RP data for each HSA              ! this is added for planning period
        call OpnIF2(1,2,ifile,cstar)        ! RP - HSA
        read(2,'(i4)') nrp
        read(2,'(a)') fmtin
        do iha=1,nha
          read(2,fmtin) kha,(fharp(iha,irp),irp=1,nrp)
          if(kha.ne.iha) stop
        end do
        close(2)        
        
        iu = 75                             ! 75 for supply monthly volume
          call OpnIF2(1,iu,ifile,cstar)     ! RW zone 
          read(iu,'(i4)') iv
          if(iv.ne.nrwz) then
            write(*,'(a,i3,a)') ' Problem with File IU',iu,
     &                          ' Stop execution'
            stop
          end if
          call skipline(iu,cstar)
          call LocTS3(iu,isyr,ismn)
        
        do iu = 76,77                       ! RP quality data, 76 TDS, 77 N
          call OpnIF2(1,iu,ifile,cstar)     ! RW zone 
          read(iu,'(i4)') iv
          if(iv.ne.nrp) then
            write(*,'(a,i3,a)') ' Problem with File IU',iu,
     &                          ' Stop execution'
            stop
          end if
          call skipline(iu,cstar)
          call LocTS3(iu,isyr,ismn)
        end do      
        
      end if  !(qopt)                                   ! <><><><><><><><><><><><><><><><><><><><><><>        
        
      ndwf = 1 ! this ndwf =1 is current option, may change with available data ***
      call OpnIF2(1,78,ifile,cstar)
      call LocTS3(78,isyr,ismn)

      if(qopt) then ! * read water quality parameters       ! <><><><><><><><><><><><><><><><><><><><><><>        
        ! manure TDS and N and fertilizet TDS and N, for irrigated and non-irrigated AG area      
        call OpnIF2(1,79,ifile,cstar)      ! irrigated, non-irrigated area manure and fertilizer TDS and N application
        call LocTS3(79,isyr,ismn)        
        
        call OpnIF2(1,80,ifile,cstar)       ! atmospheric wet and dry, TDS and N
        call LocTS3(80,isyr,ismn)
      end if  !(qopt)                                   ! <><><><><><><><><><><><><><><><><><><><><><>
        
!s9
      ! Open output files   
      
      call OpnOF2(1,11,ofile,'Monthly infiltration to RZ',file0)        ! by HSA on GWM
      call OpnOF2(1,12,ofile,'Monthly evapotranspiration',file0)        ! by HSA on GWM

      call OpnOF2(1,13,ofile,'Monthly total deep percolation  ',file0)  ! by HSA on GWM
      if(qopt) call OpnOF2(1,33,ofile,
     &              'Monthly Average TDS in Deep Percolaiton',file0) 
      if(qopt) call OpnOF2(1,53,ofile,
     &              'Monthly Average TIN in Deep Percolaiton',file0)

      write(13,'(a6,500f8.0)') 'dacsum', (xmet(ihapt(i)),i=1,nhapr)     ! delete=====================
      do iha=1,nha                                                      ! delete=====================
        xmet(iha) = haarea(iha) / 4045.                                 ! delete=====================
      end do                                                            ! delete=====================
      write(13,'(a6,500f8.0)') 'T_Area', (xmet(ihapt(i)),i=1,nhapr)     ! delete===================== 
      xmet = 0.                                                         ! delete=====================
       
      call OpnOF2(1,14,ofile,'Monthly M&I applied water to HA',file0)   ! by HSA on GWM
cccccccc      if(qopt) call OpnOF2(1,34,ofile,
cccccccc     &              'Monthly Average TDS in M&I Supply',file0)
cccccccc      if(qopt) call OpnOF2(1,54,ofile,
cccccccc     &              'Monthly Average TIN in M&I Supply',file0)
     
      call OpnOF2(1,15,ofile,'Monthly AG applied water to HA',file0)    ! by HSA on GWM
cccccccc      if(qopt) call OpnOF2(1,35,ofile,
cccccccc     &              'Monthly Average TDS in AG Supply',file0)
cccccccc      if(qopt) call OpnOF2(1,55,ofile,
cccccccc     &              'Monthly Average TIN in AG Supply',file0)
      
      if(qopt) call OpnOF2(1,16,ofile,
     &         'Monthly Recylced water to HA',file0)      ! by HSA on GWM
cccccccc      if(qopt) call OpnOF2(1,36,ofile,
cccccccc     &              'Monthly Average TDS in RW Supply',file0)
cccccccc      if(qopt) call OpnOF2(1,56,ofile,
cccccccc     &              'Monthly Average TIN in RW Supply',file0)

      call OpnOF2(1,17,ofile,'Monthly Dairy WW to HA',file0)            ! by HSA on GWM
cccccccc      if(qopt) call OpnOF2(1,37,ofile,
cccccccc     &              'Monthly Average TDS in DW Supply',file0)
cccccccc      if(qopt) call OpnOF2(1,57,ofile,
cccccccc     &              'Monthly Average TIN in DW Supply',file0)
 
      call OpnOF2(1,21,ofile,'Monthly Total Supply by WSAs',file0)      ! by WSA
cremove      call OpnOF2(1,22,ofile,'Monthly Total AG Supply',file0)           ! by HSA with AG supply      
cremove      if(qopt) call OpnOF2(1,23,ofile,'Monthly Total RW Supply',file0)           ! by HSA
cremove      call OpnOF2(1,24,ofile,'Monthly Total WW disposal',file0)         ! by HSA
      call OpnOF2(1,25,ofile,'Monthly Total AW by LU type',file0)       ! by LU
      call OpnOF2(1,26,ofile,'SM in HSA at month end',file0)            ! by HSA
      call OpnOF2(1,91,ofile,'Monthly Total Perc by LU type',file0)     ! by HSA
      call OpnOF2(1,92,ofile,'Monthly Total Rain Perc by LU',file0)     ! by HSA
      call OpnOF2(1,93,ofile,'Monthly Total A.W. Perc by LU',file0)     ! by HSA
      do iu=91,93
        write(iu,'(6x,50i8)') (ilu,ilu=1,nlu)
      end do

      if(qopt) then
        call OpnOF2(1,27,ofile,'TDS in soil zone at month end',file0)     ! by HSA
        call OpnOF2(1,28,ofile,'N in soil zone at month end',file0)       ! by HSA
 
        call OpnOF2(1,81,ofile,'TDS mass in Rain',file0)
        call OpnOF2(1,82,ofile,'N   mass in Rain',file0)
        call OpnOF2(1,83,ofile,'TDS mass in AW and M & F',file0)
        call OpnOF2(1,84,ofile,'N   mass in AW and M & F',file0)
        call OpnOF2(1,85,ofile,'TDS mass in Uptake',file0)
        call OpnOF2(1,86,ofile,'N   mass in Uptake',file0)
        call OpnOF2(1,87,ofile,'TDS mass in Deep Perc',file0)
        call OpnOF2(1,88,ofile,'N   mass in Deep Perc',file0)      
      end if

      do iu=11,15
        write(iu,fm6x8i) (ihapt(i),i=1,nhapr)
        write(iu,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)
      end do
      if (qopt) then
        write(16,fm6x8i) (ihapt(i),i=1,nhapr)
        write(16,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)
      end if
      write(17,fm6x8i) (ihapt(i),i=1,nhapr)
      write(17,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)

      if (qopt) then      
        write(33,fm6x8i) (ihapt(i),i=1,nhapr)
        write(33,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)
        write(53,fm6x8i) (ihapt(i),i=1,nhapr)
        write(53,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)
      end if
      write(21,fm6x8i) (iwsa, iwsa=1,nwsa)
cremove      write(22,fm6x8i) (kaghsa(iagh),iagh=1,nagh)    
cremove      if(qopt) write(23,fm6x8i) (irwz,irwz=1,nrwz)
cremove      write(24,fm6x8i) (idwz,idwz=1,1)  !!!!!!!!!!!!!!

      write(25,fm6x8i) (ilu, ilu=1,nlu)       
      write(26,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)      

      if(qopt) then
        do iu=27,28
          write(iu,fm6x8i) (ihapt(i),i=1,nhapr)        
          write(iu,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)      
        end do
      end if

      if(qopt) then
        do iu=81,88
          write(iu,fm6x8i) (ihapt(i),i=1,nhapr)        
          write(iu,fma6a8) 'yyyym ',(haid8(ihapt(i)),i=1,nhapr)      
        end do
      end if
      
!s11
      ! initial rootzone soil moisture storage
      do iha=1,nha
        xmsto(iha) = 0.
      end do
      do kp=1,kpg
        iha = ihls(kp,1)
        ilu = ihls(kp,2)
        ist = ihls(kp,3)
        if(ist.gt.0) then  ! ignore impervious area 
          xmsto(iha) = xmsto(iha) + 
     &                (rsm(kp,1)*dac(kp,1) + 
     &                 rsm(kp,2)*dac(kp,2)) / 12.
        end if
      end do
      icy = isyr
      icm = ismn-1
	  if(icm.eq.0) then
	    icy = icy-1
		icm = 12
      end if
      write(*,*) icy,icm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     write(26,8081) icy,icm,(xmsto(ihapt(i)),i=1,nhapr)
      write(26,fmym82) icy,icm,(xmsto(ihapt(i)),i=1,nhapr)

      do 500 ipy=1,npor !===== YEAR LOOP ===============================

      icy = isyr+ipy-1     
      write(*,'(/,i5)') icy

      do 400 ifm=1,12   !===== MONTH LOOP ==================================================

      icm = ifm + ismn - 1
      if(icm.gt.12) then
        icm = icm - 12
        if(icm.eq.1) icy = icy+1
      end if
  
      if(debug) write(*,*) '  pass m1a  '      
      Qrwz = Frwz    ! save any left-over RW supply
      if(qopt) then       !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! read monthly water quality data
      read(71,fmym80) iyr,imn,(cTmws(iwsa), iwsa=1,nwsa)      ! Municipal WS TDS (mg/L) 
      if(iyr.ne.icy .or. imn.ne.icm) stop
      if(debug) write(*,*) '  pass m1a1  '   
      read(72,fmym80) iyr,imn,(cNmws(iwsa), iwsa=1,nwsa)      ! Municipal WS N (mg/L)
      if(iyr.ne.icy .or. imn.ne.icm) stop
      if(debug) write(*,*) '  pass m1a2  '   
      read(73,fmym80) iyr,imn,(cTaws(iagh), iagh=1,nagh)     ! GW AG supply TDS (mg/L) by HSA
      if(iyr.ne.icy .or. imn.ne.icm) stop
      if(debug) write(*,*) '  pass m1a3  '   
      read(74,fmym80) iyr,imn,(cNaws(iagh), iagh=1,nagh)     ! GW AG supply N (mg/L) by HSA
      if(iyr.ne.icy .or. imn.ne.icm) stop
      if(debug) write(*,*) '  pass m1a4  '   
      read(75,fmym80) iyr,imn,((Frwz(irwz,irwt),irwz=1,nrwz),! RW supply (af/month) by RW Zone
     &                                irwt=1,2)              ! irwt 1 urban/park, 2 ag 
      if(iyr.ne.icy .or. imn.ne.icm) stop
      if(debug) write(*,*) '  pass m1b  '            
      read(76,fmym80) iyr,imn,(cTrp(irp), irp=1,nrp)     ! RW supply TDS (mg/L) by RW Zone
      if(iyr.ne.icy .or. imn.ne.icm) stop
      read(77,fmym80) iyr,imn,(cNrp(irp), irp=1,nrp)     ! RW supply N (mg/L) by RW Zone
      if(iyr.ne.icy .or. imn.ne.icm) stop
      end if ! qopt  !<><><><><><><><><><><><>

      read(78,fmym80) iyr,imn,dwflw, dwNinc, dwTinc          ! dairy ww spray (acre-in/(acre-day), N and TDS increment (mg/L)
      if(qopt. and. dwTinc.eq.0.) dwTinc = cTdwinc
      if(qopt. and. dwNinc.eq.0.) dwNinc = cNdwinc
      if(iyr.ne.icy .or. imn.ne.icm) stop

      if(qopt) then !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                 ! constant for a month for all dairy wash water appplication area
      read(79,fmym80) iyr,imn,(((tnar(itn,imf,irg),        
     &                           itn=1,2),imf=1,2),irg=1,2)
            ! TDS and N application in Ag area by manure and fertilizer application
            ! itn = 1 TDS       2 N
            ! imf = 1 manure    2 fertilizer
            ! irg = 1 irrigated 2 Non-irrigated
      if(iyr.ne.icy .or. imn.ne.icm) stop          
      if(debug) write(*,*) '  pass m1c  '             

      read(80,fmym80) iyr,imn,((atmdep(i,j), i=1,2),j=1,2)
      if(iyr.ne.icy .or. imn.ne.icm) stop            
      if(debug) write(*,*) '  pass y1a  '
      
      ! set manure and fertilzier application rates for the month
      ! ilu = 1,5                         51,non-irrig ag
      ! ilu = 2,3,4,20                    52,irrigated ag
      ! ilu = 7,10,11,12,15,16,17,18,19   53,urban
      ! ilu = 6,8,9,13,14,                54,no manure and fertilizer application
      dmanu = 0.
      ftapp = 0.
      do ilu=1,20
        goto (51,52,52,52,51,
     &        54,53,54,54,53,
     &        53,53,54,54,53,
     &        53,53,53,53,52) ilu
 51     dmanu(1,ilu) = tnar(1,1,2)    ! non-irrigated ag
        dmanu(2,ilu) = tnar(2,1,2)
        ftapp(1,ilu) = tnar(1,2,2)
        ftapp(2,ilu) = tnar(2,2,2)        
        goto 54
 52     dmanu(1,ilu) = tnar(1,1,1)     ! irrigated ag
        dmanu(2,ilu) = tnar(2,1,1)
        ftapp(1,ilu) = tnar(1,2,1)
        ftapp(2,ilu) = tnar(2,2,1)
        goto 54
 53     ftapp(1,ilu) = ftappC(1,ilu,icm) ! urban fertilizer only, no manure
        ftapp(2,ilu) = ftappC(2,ilu,icm)
 54     continue        
      end do
      if(debug) write(*,*) '  pass y1b  '
      end if ! qopt <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      ! reset variables
      xminf = 0.    ! (iha)     ! inflow
      xmet = 0.     ! (iha)     ! evapotranspiration
      xmaw = 0.     ! (iha,1)   !  1 M&I
                    ! (iha,2)   !  2 AG
                    ! (iha,3)   !  3 RW
                    ! (iha,4)   !  4 Dairy WW
      xmpc = 0.     ! (iha)     !    monthly applied water (af) 1 MW 2 AW 3 RW 4 dairy ww+AW
      xmsto = 0.    ! (iha)     !   
      xmpcr = 0.    ! (iha) 
      xmpca = 0.    ! (iha) 

      xmpcl  = 0.   !(ilu)
      xmawl  = 0.   !(ilu)
      xmpclr = 0.   !(ilu)
      xmpcla = 0.   !(ilu)

      xmmws = 0.    ! (MXWS),
c      xmgws         ! since GW supply its own HSA, don't need this      
c      xmrws = 0.    ! (MXRW)
c      xmdws = 0.    ! one dairy washwater to collect all for water budget

      if(qopt) then       !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      cTxmpc = 0.   ! (MXHA), 
      cNxmpc = 0.   ! (MXHA),
      cTxmaw = 0.   ! (MXHA,3),
      cNxmaw = 0.   ! (MXHA,3)      

ccccccccccccccccccc      Qmwh = 0.     ! monthly MW supply to HSAs
ccccccccccccccccccc      Qrwh = 0.     ! monthly RW supply to HSAs
ccccccccccccccccccc      Qawh = 0.     ! monthly GW supply to HSAs (Ag water)
ccccccccccccccccccc      Qdwh = 0.     ! monthly Dairy Wash W supply to HSAs


      ! reset
      cTmwh = 0.
      cNmwh = 0.
      cTrwh = 0. 
      cNrwh = 0. 
      cTawh = 0.
      cNawh = 0.
      
      TNtr = 0.
      
      if(debug) write(*,*) '  pass y1c  '      
      ! quality of supply water to the HSAs
      ! Municipal Water quality
      do iha=1,nha
        ! average water supply quality to the HSAs
        do iwsa = 1,nwsa      
          cTmwh(iha) = cTmwh(iha) + fhawsa(iha,iwsa) * cTmws(iwsa)
          cNmwh(iha) = cNmwh(iha) + fhawsa(iha,iwsa) * cNmws(iwsa)
        end do
        ! demand will be calculated
      end do 
      ! Local GW for ag use
      cTawh = cTaws(nagh)
      cNawh = cNaws(nagh)
      do iagh = 1,nagh-1
        cTawh(kaghsa(iagh)) = cTaws(iagh)   ! assign input values
        cNawh(kaghsa(iagh)) = cNaws(iagh)
      end do
      ! ag demand will be calculated
      ! RW
      ! RW supply to rwz = frwz(irwz,2)     1 for urban, 2 for af
!v6      do iha=1,nha
!v6        krw = krwhsa(iha)
!v6        if(krw.gt.0) cTrwh(iha) = cTrwz(krw)    ! assign input values
!v6        if(krw.gt.0) cNrwh(iha) = cNrwz(krw)
!v6      end do  
      do iha=1,nha
        do irp=1,nrp
          cTrwh(iha) = cTrwh(iha) + cTrp(irp) * fharp(iha,irp)
          cNrwh(iha) = cNrwh(iha) + cNrp(irp) * fharp(iha,irp)
        end do
      end do
      write(6,'(i4,i2,500f8.0)') icy,icm,
     &        (ctrwh(kaghsa(iagh)),iagh=1,nagh)
      
      end if ! qopt <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      if(debug) write(*,*) '  pass y1d  '          

      mnday = mndays(icy,icm)

	   
      do 300 icd = 1,mnday    !===== DAY LOOP ============================

      read(3,fmtet,end=990) m10,m20,m30,(dget(i),i=1,net)             ! ok for daily ET on HSAs
      if(m10.ne.icy .or. m20.ne.icm .or. m30.ne.icd) stop  
      read(4,'(i4,2i2,/,(25f8.0))',end=990) m11,m21,m31,
     &                    ((rinf(kp,i),i=1,4),kp=1,kpg)  ! (inches/day)
      if(m11.ne.icy .or. m21.ne.icm .or. m31.ne.icd) stop  
                ! rinf(kp,1) effective rainfall to irrigated pervious area
                ! rinf(kp,2)      "       "     to nonirrigated pervious area (barren area)
                ! rinf(kp,3) infiltration       in irrigated pervious area
                ! rinf(kp,4)      "             in nonirrigated pervious area (barren area)

      if(debug) write(*,*) '  pass y2a  '  
      rTleach = 0.
      rNleach = 0.
                
      do 200 kp=1,kpg   !===== HA/LU/ST LOOP ============================
        iha = ihls(kp,1)
        ilu = ihls(kp,2)
        ist = ihls(kp,3)
        iz  = ihaiz(iha)
        if(debug) write(*,*) '  pass y2b  '  
        RefETo = 0.
        if(net.lt.nha) then
          do i=1,net
            RefETo = RefETo + etwt(iha,i) * dget(i)
          end do
        else
          RefETo = dget(iha)
        end if

        etp = cropcoef(iz,ilu,icm) * RefETo   ! maximum ET potential

        if(.not.ihach(iha)) goto 195    ! if the HSA is not on the RZ simulaiton list, skip
        if(dac(kp,1).eq.0.) goto 100
        
        ! irrigated pervious area *********************************************************
        dcg = dac(kp,1)                                          ! * gwbr(iha) is applied in RF model
        dcg12 = dcg / 12.
        
        rapw(kp) = 0.               ! applied water (in)
        rdp(kp,1) = 0.              ! deep perc of irrigated polygon k
        rdpr(kp,1) = 0.             ! deep perc of rainfall  (in)
        rdpa(kp,1) = 0.             ! deep perc of applied water  (in)
        deficit = 0.                ! soil moisture deficit  (in)

        qmw = 0.
        qgw = 0.
        qrw = 0.        
        qdw = 0.                                     ! (af)
        qdwin = 0.                                   ! (in)
        
        if(debug) write(*,*) '  pass y3a  '  

        ! update soil moisture
        rsmn = rsm(kp,1) + rinf(kp,3)                ! (in) add rainfall infiltration 

        if(ilu.eq.20) then ! apply dairy ww
          qdwin = dwflw                              ! (in) 
          rsmn = rsmn + qdwin
          qdw = qdwin * dcg12                        ! (af)
          xmaw(iha,4) = xmaw(iha,4) + qdw            ! (af)
c          xmdws = xmdws + qdw
          xmawl(ilu) = xmawl(ilu) + qdw
        end if
        if(debug) write(*,*) '  pass y3b  '  
       
        if(rsmn.gt.rfc(kp,1)) then                    ! land is oversaturated with rain
          rdp(kp,1) = rsmn - rfc(kp,1)                     ! rain and wash water deep percolation *** (in)
cccc          rval2 = rinf(kp,3)+qdwin
          if(debug) write(*,*) '  pass y3b1  ',iha,ilu,kp,      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     &          rsm(kp,1),rsmn, rfc(kp,1),
     &          rinf(kp,3),qdwin  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!,rval2 

          rdpr(kp,1) = rdp(kp,1) * rinf(kp,3)/(rinf(kp,3)+qdwin)
          rdpa(kp,1) = rdp(kp,1) - rdpr(kp,1)
          if(debug) write(*,*) '  pass y3b1a  '

          rsmn = rfc(kp,1)
          etd = etp
          
        else
          if(apwrq(ilu).and.irmn(iz,ilu,icm)) then    ! (in) irrigated area, irrigated month
            if(rsmn.lt.rcsm(kp,2,1)) then             ! (in) soil is dry below irrigation requirement
                                                      !      this rsmn is already contains rinf(kp,3) and dwflw.
              if(debug) write(*,*) '  pass y3b2  '  
              deficit = rcsm(kp,3,1) - rsmn           ! (in)
              rapw(kp) = deficit * apwef(iz,ilu)      ! (in) irrigation efficiency (converted), applied water quantity.
              rapwaf = rapw(kp) * dcg12               ! (af) irrigation water volume (acre-ft)
              xmawl(ilu) = xmawl(ilu) + rapwaf
              if(debug) write(*,*) '  pass y3b2a  '  
  
              ! what is this applied water
              
              irwz = krwhsa(iha)                       ! rw supply zone
              irwt = isrc(iz,ilu)                      ! krwt(ilu)  1 for urban use, 2 for ag use 3 phreatophyte
              if(irwz.gt.0) then                       ! use up RW
                qrw = 0.
c               if(irwt.gt.0) qrw = min(Frwz(irwz,irwt), rapwaf)
                if(irwt.eq.1 .or. irwt.eq.2) grw = 
     &                        min(Frwz(irwz,irwt), rapwaf)
                if(qrw.le.0.) qrw = 0.
                Frwz(irwz,irwt) = Frwz(irwz,irwt) - qrw        ! remaining RW for the month for the RW supply zone
                rapwaf = rapwaf - qrw
 
                xmaw(iha,3) = xmaw(iha,3) + qrw
cremove                xmrws(irwz) = xmrws(irwz) + qrw
              end if
              if(rapwaf.gt.0.) then 
                if (isrc(iz,ilu).eq.1) then             ! municipal WSA
                  qmw = rapwaf
                  xmaw(iha,1) = xmaw(iha,1) + qmw
                  do iwsa=1,nwsa
                    xmmws(iwsa) = xmmws(iwsa) + qmw * fhawsa(iha,iwsa)
                  end do
                elseif(isrc(iz,ilu).eq.2) then          ! local GW
                  qgw = rapwaf
                  xmaw(iha,2) = xmaw(iha,2) + qgw
                  ! since local gw supplies only HSA, don't need to calculate                  
                end if                
                rapwaf = 0.              
              end if !rapwaf  
              if(debug) write(*,*) '  pass y3b5  '  
                
              rsmn = rsmn + rapw(kp)                    ! (in)
              rdp(kp,1) = rsmn - rfc(kp,1)              ! (in) applied water percolation
              rval = rinf(kp,3) + qdwin + rapw(kp)
              rdpr(kp,1) = rdp(kp,1) * rinf(kp,3) / rval
              rdpa(kp,1) = rdp(kp,1) - rdpr(kp,1)
              ! update rsmn
              rsmn = rfc(kp,1)
              
              if(debug) write(*,*) '  pass y3b6  ' 
              if(debug) write(*,*)  iha,ilu,kp,rsm(kp,1),
     &               rfc(kp,1),rcsm(kp,2,1),rcsm(kp,3,1),
     &               rapw(kp),rdp(kp,1),rsmn                   !!!!!!!!!!!! 
            end if  ! rsmn.lt.rcsm(kp,2,1)
          end if  ! apwrq(ilu)
        end if  ! rsmn
        if(debug) write(*,*) '  pass y3c  '  
        
        
c       ! total inflow
c            rinf(kp,3) + qdwin + rapw(kp)  (in)
        ! total applied water
        tapwin = qdwin + rapw(kp)                   ! (in)
        tapwaf = qmw + qgw + qrw + qdw              ! (af/day)
        qdwout = 0.
        qmwout = 0.
        qrwout = 0.
        qgwout = 0.
        if(tapwaf.gt.0.) then
          qmwout = rdpa(kp,1) * qmw / tapwaf            ! (in) ignore volume in storage
          qgwout = rdpa(kp,1) * qgw / tapwaf
          qrwout = rdpa(kp,1) * qrw / tapwaf
          qdwout = rdpa(kp,1) * qdw / tapwaf
        end if
        if(debug) write(*,*) '  pass y3c1  '                 
!dl1
        ! evapotranspiration
        if(rsmn.gt.rcsm(kp,1,1)) then                  ! this rsmn is storage + inflow - deep perc
          etd = etp
        elseif(rsmn.gt.rwp(kp,1)) then
          etd = etp * (rsmn-rwp(kp,1)) / (rcsm(kp,1,1) - rwp(kp,1))
        else
          etd = 0.  ! it can happen in non-irrigated area
        end if

        rval = max(0., rsmn - rwp(kp,1))
        etd = min(etd,rval)
        rsmn = rsmn - etd

        ! total water mass for concentration calculation
        rsma = rsm(kp,1) + rinf(kp,3) + qdwin + rapw(kp)          ! (in)
        
        rsm(kp,1) = rsmn
        reta(kp,1) = etd                                         ! (in) actual et
        rdp(kp,1) = rdpr(kp,1) + rdpa(kp,1)                        ! (in) total deep perc of rainfall and applied water

        
        if(debug) write(*,*) '  pass y3d  '  
        

        if(qopt) then       !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        ! update TDS and nitrogen +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! soil surface
        sTmass(kp,1) = sTmass(kp,1) + atmdep(1,1)                 ! dry atmospheric deposit (lb/ac))
        sNmass(kp,1) = sNmass(kp,1) + atmdep(2,1)
        ! loss to atmosphere
        qsNloss = sNmass(kp,1) * sNloss(1)
        qsTloss = qsNloss * 62./14.
        sTmass(kp,1) = sTmass(kp,1) - qsTloss
        sNmass(kp,1) = sNmass(kp,1) - qsNloss
        ! transfer of surface material to the rootzone
        sN2rN = sNmass(kp,1) * sNloss(2)
c        sT2rT = sN2rN * 62./14.
        sT2rT = sTmass(kp,1) * sNloss(2)
        sTmass(kp,1) = sTmass(kp,1) - sT2rT                        ! (lb/ac)
        sNmass(kp,1) = sNmass(kp,1) - sN2rN          

        !
        ! applied water quality

        ! rootzone =======================================
        ! transfer from surface        
        rTmass(kp,1) = rTmass(kp,1) + sT2rT                   ! (lb/ac), this is dry atmdep only
        rNmass(kp,1) = rNmass(kp,1) + sN2rN

        TNtr(iha,1,1)= TNtr(iha,1,1) + sT2rT * dcg
        TNtr(iha,1,2)= TNtr(iha,1,2) + sN2rN * dcg
        
        if(rinf(kp,3).gt.0.) then   ! precipitation infiltration
          Ctran = rinf(kp,3) * atmdep(1,2) * cfact1    ! (lb/ac)
          rTmass(kp,1) = rTmass(kp,1) + Ctran
          TNtr(iha,1,1)= TNtr(iha,1,1) + Ctran * dcg
                   
          Ctran = rinf(kp,3) * atmdep(2,2) * cfact1
          rNmass(kp,1) = rNmass(kp,1) + Ctran
          TNtr(iha,1,2)= TNtr(iha,1,2) + Ctran * dcg
        end if
        
        ! applied water
		cTdwh = 0.
		cNdwh = 0.
        if(qdw.gt.0.) then		
          cTdwh = cTawh(iha) + dwTinc
		  cNdwh = cNawh(iha) + dwNinc
		end if

        cTawt = 0.
        cNawt = 0.
          if(debug) write(*,*) '  pass y3e' 
          if(debug) write(*,'(i5,i2,i3,2i4,i5,5(e10.3,f8.1))')    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     &      icy,icm,icd,iha,ilu,kp,                  
     &      qmw,ctmwh(iha),qgw,ctawh(iha),qrw,ctrwh(iha),qdw,ctdwh,
     &      tapwaf,dac(kp,1)     
        if(tapwaf.gt.0.) then ! applied water conc

          cTawt = (qdw * cTdwh + qmw * cTmwh(iha)
     &           + qrw * cTrwh(iha) + qgw * cTawh(iha) ) / tapwaf
          cNawt = (qdw * cNdwh + qmw * cNmwh(iha)
     &           + qrw * cNrwh(iha) + qgw * cNawh(iha) ) / tapwaf
          if(debug) write(*,'(2f8.1)') cTawt, cNawt
        end if
 
        !applied water        
        Ctran = rapw(kp) * cTawt * cfact1     ! (lb/acre/day)
        rTmass(kp,1) = rTmass(kp,1) + Ctran
        TNtr(iha,2,1) = TNtr(iha,2,1) + Ctran * dcg
        Ctran = rapw(kp) * cNawt * cfact1
        rNmass(kp,1) = rNmass(kp,1) + Ctran
        TNtr(iha,2,2) = TNtr(iha,2,2) + Ctran * dcg
        
        ! dairy manure application        
        rTmass(kp,1) = rTmass(kp,1) + dmanu(1,ilu)                  ! daily application (lb/ac/day)
        rNmass(kp,1) = rNmass(kp,1) + dmanu(2,ilu)
        
        ! fertilizer application
        rTmass(kp,1) = rTmass(kp,1) + ftapp(1,ilu)
        rNmass(kp,1) = rNmass(kp,1) + ftapp(2,ilu)
               
        TNtr(iha,2,1) = TNtr(iha,2,1) + 
     &                  (dmanu(1,ilu) + ftapp(1,ilu)) * dcg
        TNtr(iha,2,2) = TNtr(iha,2,2) + 
     &                  (dmanu(2,ilu) + ftapp(2,ilu)) * dcg
                
        ! nitrogen uptake
        cNuptake = min(rNmass(kp,1), rNuptake(ilu,icm))
        cTuptake = min(rTmass(kp,1), cNuptake*62./14.)  ! assume nitrogen salt loss
        
        rTmass(kp,1) = rTmass(kp,1) - cTuptake 
        rNmass(kp,1) = rNmass(kp,1) - cNuptake
        
        TNtr(iha,3,1) = TNtr(iha,3,1) + cTuptake * dcg
        TNtr(iha,3,2) = TNtr(iha,3,2) + cNuptake * dcg
        
        ! deep perc
        ! the last rTmass and rNmass is the total mass subject to leach;
        rTleach(kp,1) = rTmass(kp,1) * rdp(kp,1) / rsma      ! (lb/ac)
        rNleach(kp,1) = rNmass(kp,1) * rdp(kp,1) / rsma
        rTmass(kp,1) = rTmass(kp,1) - rTleach(kp,1)
        rNmass(kp,1) = rNmass(kp,1) - rNleach(kp,1)
        
        TNtr(iha,4,1) = TNtr(iha,4,1) + rTleach(kp,1) * dcg
        TNtr(iha,4,2) = TNtr(iha,4,2) + rNleach(kp,1) * dcg
        
c        if(debug) write(*,*) '  pass y3g  '          
      end if ! qopt <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        ! monthly summary ==========================================================================


        xminf(iha) = xminf(iha) + rinf(kp,3) * dcg12             ! infiltration of rainfall, in (acre-ft)
        xmet(iha) = xmet(iha) + reta(kp,1) * dcg12               ! ET by HA
        xmpc(iha) = xmpc(iha) + rdp(kp,1) * dcg12                ! deep percolation by HA in (acre-ft)   
        if(debug) write(*,'(i5,2i3,2i5,a5,5e10.3)') icy,icm,icd,iha,
     *      kp,'  irr',xmpc(iha),rdp(kp,1),rdpr(kp,1),rdpa(kp,1),dcg12
        if(qopt) then  
          cTxmpc(iha) = cTxmpc(iha) + rTleach(kp,1) * dcg        ! in (lbs)
          cNxmpc(iha) = cNxmpc(iha) + rNleach(kp,1) * dcg        ! 
        end if

        xmpcr(iha) = xmpcr(iha) + rdpr(kp,1) * dcg12             ! deep perc of rainfall 
        xmpca(iha) = xmpca(iha) + rdpa(kp,1) * dcg12             ! deep perc of applied water
        
        ! 7/1/19
        xmpcl(ilu)  =  xmpcl(ilu)  + rdp(kp,1) * dcg12 
        xmpclr(ilu) =  xmpclr(ilu) + rdpr(kp,1) * dcg12
        xmpcla(ilu) =  xmpcla(ilu) + rdpa(kp,1) * dcg12

		
 100    continue
        ! non-irrigated area
        if(debug) write(*,*) '  pass y3h  ' 
                      
        if(dac(kp,2).eq.0.) goto 195 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! non-irrigaed area
        ! non-irrigated pervious area
        dcg = dac(kp,2)                                          ! * gwbr(iha) is applied in RF model
        dcg12 = dcg / 12. 

        rsmn = rsm(kp,2) + rinf(kp,4)                 ! rainfall infiltration

        qdw = 0.
        if(ilu.eq.20) then ! apply dairy ww
          qdwin = dwflw
          rsmn = rsmn + qdwin                        ! dairy wash water
          qdw = qdwin * dcg12
          xmaw(iha,4) = xmaw(iha,4) + qdw            ! (af)
c          xmdws = xmdws + qdw          
          xmawl(ilu) = xmawl(ilu) + qdw
        end if

        rdp(kp,2) = 0.
        rdpr(kp,2) = 0.                              ! deep perc of rain water
        deficit = 0.

        rsma = rsmn                                  ! volume to calculate concentration
        if(rsmn.gt.rfc(kp,2)) then                   ! land is oversaturated with rain
          rdp(kp,2) = rsmn - rfc(kp,2)
          rsmn = rfc(kp,2)
          etd = etp
          rdpr(kp,2) = rdp(kp,2) * rinf(kp,4) / (rinf(kp,4) + qdwin)
          rdpa(kp,2) = rdp(kp,2) - rdpr(kp,2)
        elseif(rsmn.gt.rcsm(kp,1,2)) then
          etd = etp
        elseif(rsmn.gt.rwp(kp,2)) then
          etd = etp * (rsmn-rwp(kp,2)) / (rcsm(kp,1,2) - rwp(kp,2))
        else
          etd = 0.
        end if

        rval = max(0., rsmn - rwp(kp,2))
        etd = min(rval,etd)              ! added 10/9/18
        rsmn = rsmn - etd                ! evapotranspiration
        rsma = rsma - etd
        
        rsm(kp,2) = rsmn
        reta(kp,2) = etd                 ! actual et

        rdp(kp,2) = rdp(kp,2)
        if(debug) write(*,*) '  pass y3i  ' 

        
        
        if (qopt) then  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>             
        ! update TDS and nitrogen +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! soil surface
        sTmass(kp,2) = sTmass(kp,2) + atmdep(1,1)              ! dry deposit, TDS
        sNmass(kp,2) = sNmass(kp,2) + atmdep(2,1)              ! dry deposit, N
                      if(debug) write(*,*) '  pass y3j  ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! loss to atmosphere
        qsNloss = sNmass(kp,2) * sNloss(1)
        qsTloss = qsNloss * 62./14.
        sTmass(kp,2) = sTmass(kp,2) - qsTloss
        sNmass(kp,2) = sNmass(kp,2) - qsNloss
        
        ! transfer of surface material to the rootzone
        sN2rN = sNmass(kp,2) * sNloss(2)
c        sT2rT = sN2rN * 62./14.
        sT2rT = sTmass(kp,2) * sNloss(2)
        sTmass(kp,2) = sTmass(kp,1) - sT2rT                                
        sNmass(kp,2) = sNmass(kp,1) - sN2rN         
                      if(debug) write(*,*) '  pass y3j1  '
        ! root zone ================================================================
        
        rTmass(kp,2) = rTmass(kp,2) + sT2rT  ! transfer from surface
        rNmass(kp,2) = rNmass(kp,2) + sN2rN
        
        TNtr(iha,1,1)= TNtr(iha,1,1) + sT2rT * dcg
        TNtr(iha,1,2)= TNtr(iha,1,2) + sN2rN * dcg
        
        if(rinf(kp,4).gt.0.) then           ! rainfall infiltration
          Ctran = rinf(kp,4) * atmdep(1,2) * cfact1      !(lb/ac)
          rTmass(kp,2) = rTmass(kp,2) + Ctran
          TNtr(iha,1,1)= TNtr(iha,1,1) + Ctran * dcg
          Ctran = rinf(kp,4) * atmdep(2,2) * cfact1      
          rNmass(kp,2) = rNmass(kp,2) + Ctran
          TNtr(iha,1,2)= TNtr(iha,1,2) + Ctran * dcg
        end if
        
        if(ilu.eq.20) then                          ! dairy wash water
          Ctran = dwflw * (cTawh(iha) + dwTinc) * cfact1
          rTmass(kp,2) = rTmass(kp,2) + Ctran
          TNtr(iha,2,1) = TNtr(iha,2,1) + Ctran * dcg               
          Ctran = dwflw * (cNawh(iha) + dwNinc) * cfact1
          rNmass(kp,2) = rNmass(kp,2) + Ctran
          TNtr(iha,2,2) = TNtr(iha,2,2) + Ctran * dcg                
        end if
        if(debug) write(*,*) '  pass y3j2  '        
        ! on non-irrigated area, no manure and fertilizer
        
        ! nitrogen uptake
        if(debug) write(*,*) '  pass y3j3  ' 
      
        ! nitrogen uptake
        cNuptake = min(rNmass(kp,2), rNuptake(klubar,icm))
        cTuptake = min(rTmass(kp,2), cNuptake*62./14.)  ! assume nitrogen salt loss
        
        rTmass(kp,2) = rTmass(kp,2) - cTuptake 
        rNmass(kp,2) = rNmass(kp,2) - cNuptake
        
        TNtr(iha,3,1) = TNtr(iha,3,1) + cTuptake * dcg
        TNtr(iha,3,2) = TNtr(iha,3,2) + cNuptake * dcg

    
        ! deep perc
        ! the last rTmass and rNmass is the total mass subject to leach;
        if(debug) write(*,*) '  pass y3j4  '  
        if(debug) write(*,'(4i5,e15.6)') kp,iha,ilu,ist,rsmn                      ! problem  rsmn is negative -----------------------------------

        rTleach(kp,2) = 0.
        rNleach(kp,2) = 0.
        if(rsma.gt.0.) then
          rTleach(kp,2) = rTmass(kp,2) * rdp(kp,2) / rsma   ! (lb/ac)
          rNleach(kp,2) = rNmass(kp,2) * rdp(kp,2) / rsma
        end if  
        if(debug) write(*,*) '  pass y3j5  '         
        rTmass(kp,2) = rTmass(kp,2) - rTleach(kp,2)
        rNmass(kp,2) = rNmass(kp,2) - rNleach(kp,2)

        TNtr(iha,4,1) = TNtr(iha,4,1) + rTleach(kp,2) * dcg
        TNtr(iha,4,2) = TNtr(iha,4,2) + rNleach(kp,2) * dcg     


        if(debug) write(*,*) '  pass y3k  ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end if ! qopt <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>     
        
        ! monthly summary =========================================================================

        xminf(iha) = xminf(iha) + rinf(kp,4) * dcg12                ! total infiltration by HA
        xmpc(iha) = xmpc(iha) + rdp(kp,2) * dcg12   

        if(debug) write(*,'(i5,2i3,2i5,a5,2(3e10.3,10x))') icy,icm,
     *      icd,iha,kp,'   nr',xmpc(iha),rdp(kp,2),rdpr(kp,2),dcg12
     
        xmpcr(iha) = xmpcr(iha) + rdpr(kp,2) * dcg12
        xmet(iha) = xmet(iha) + reta(kp,2) * dcg12                  ! ET by HA
        if(debug) then
               write(*,*) '  pass y3k1  ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        write(*,*) ctxmpc(iha),rtleach(kp,2),cNxmpc(iha),rNleach(kp,2)
        end if
        if(qopt) then
          cTxmpc(iha) = cTxmpc(iha) +  rTleach(kp,2) * dcg       ! (lbs)
          cNxmpc(iha) = cNxmpc(iha) +  rNleach(kp,2) * dcg
        end if
                      if(debug) write(*,*) '  pass y3k2  ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
 195    continue
 
!dl3
c        if(ndebug.eq.2) then
c          rv1 = rdp(kp,1) * dac(kp,1)
c          rv2 = rdp(kp,2) * dac(kp,2)
c          rv3 = rdpr(kp,1)+rdpa(kp,2)
c          write(44,'(i4,2i2.2,4i4,5f7.4,f7.2,4f7.4,f7.2)')
c     &      icy,icm,icd,iha,ilu,ist,k,
c     &      rinf(kp,3),rapw(kp,1),reta(kp,1),rv3,rsm(kp,1),rv1,
c     &      rinf(kp,4),reta(kp,2),rdpr(kp,2),rsm(kp,2),rv2
c        end if
         if(debug) write(*,*) '  pass y3m  ' 


 200  continue ! kp - HA/LU/ST loop
c        if(debug .and.
c     &     icy.ge.idsy .and. icy.le.idey)
c     &       write(94,'(i4,2i2,10f8.2)') icy,icm,icd,
c     &                rinf(jpg,3),
c     &           rapw(jpg,1),reta(jpg,1),rdp(jpg,1),rsm(jpg,1),
c     &                rinf(jpg,4),
c     &           rapw(jpg,2),reta(jpg,2),rdp(jpg,2),rsm(jpg,2)
c  
 300  continue ! day loop


      ! rootzone soil moisture storage
      do kp=1,kpg
        iha = ihls(kp,1)
        ilu = ihls(kp,2)
        ist = ihls(kp,3)
        if(ist.gt.0) then  ! ignore impervious area
          xmsto(iha) = xmsto(iha) +
     &        (rsm(kp,1) * dac(kp,1) + rsm(kp,2) * dac(kp,2)) / 12.
        end if
      end do

      do i=1,nhapr                                      ! added to eliminate truncation error blow-up
        if(xmpc(ihapt(i)).lt.0.1) xmpc(ihapt(i))= 0.    ! added to eliminate truncation error blow-up
	  end do                                            ! added to eliminate truncation error blow-up
	  
      ! print monthly variables
      write(11,fmym82) icy,icm,(xminf(ihapt(i)),i=1,nhapr)
      write(12,fmym82) icy,icm,(xmet(ihapt(i)),i=1,nhapr)
      write(13,fmym81) icy,icm,(xmpc(ihapt(i)),i=1,nhapr)
      write(14,fmym81) icy,icm,(xmaw(ihapt(i),1),i=1,nhapr)
      write(15,fmym81) icy,icm,(xmaw(ihapt(i),2),i=1,nhapr)
      if(qopt) write(16,fmym81) icy,icm,(xmaw(ihapt(i),3),i=1,nhapr)
      write(17,fmym82) icy,icm,(xmaw(ihapt(i),4),i=1,nhapr)

      write(21,fmym81) icy,icm,(xmmws(iwsa),iwsa=1,nwsa)
cremove      write(22,fmym82) icy,icm,(xmgws(ihapt(i)),i=1,nhapr)
cremove      cremoveif(qopt) write(23,fmym82) icy,icm,(xmrws(irwz),irwz=1,nrwz)
cremove      write(24,fmym82) icy,icm,xmdws
      write(25,fmym81) icy,icm,(xmawl(ilu), ilu=1,nlu)
      write(26,fmym82) icy,icm,(xmsto(ihapt(i)),i=1,nhapr)

      write(91,'(i4,i2.2,25f8.1)') icy,icm,(xmpcl(ilu),ilu=1,nlu)
      write(92,'(i4,i2.2,25f8.1)') icy,icm,(xmpclr(ilu),ilu=1,nlu)
      write(93,'(i4,i2.2,25f8.1)') icy,icm,(xmpcla(ilu),ilu=1,nlu)  
    

      if (qopt) then  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
      do iha=1,nha 
        if(xmpc(iha).gt.0.) then
          cTxmpc(iha) = cTxmpc(iha) / xmpc(iha) * cfact2
          cNxmpc(iha) = cNxmpc(iha) / xmpc(iha) * cfact2
          if(cTxmpc(iha).lt.80.0) cTxmpc(iha) = 80.0   ! for leaching of rainfall
          if(cNxmpc(iha).lt.0.5) cNxmpc(iha) = 0.5          
        else
           cTxmpc(iha) = 0.
           cNxmpc(iha) = 0.
        end if
      end do
      write(33,fmym81) icy,icm,(cTxmpc(ihapt(i)),i=1,nhapr)
      write(53,fmym82) icy,icm,(cNxmpc(ihapt(i)),i=1,nhapr)
         if(debug) write(*,*) '  pass y3n  '       



      write(81,'(i4,i2.2,500f10.0)') icy,icm,
     &                            (TNtr(ihapt(i),1,1),i=1,nhapr)
      write(82,fmym80) icy,icm,(TNtr(ihapt(i),1,2),i=1,nhapr)
      write(83,'(i4,i2.2,500f10.0)') icy,icm,
     &                            (TNtr(ihapt(i),2,1),i=1,nhapr)
      write(84,fmym80) icy,icm,(TNtr(ihapt(i),2,2),i=1,nhapr)
      write(85,'(i4,i2.2,500f10.0)') icy,icm,
     &                            (TNtr(ihapt(i),3,1),i=1,nhapr)
      write(86,fmym80) icy,icm,(TNtr(ihapt(i),3,2),i=1,nhapr)
      write(87,'(i4,i2.2,500f10.0)') icy,icm,
     &                            (TNtr(ihapt(i),4,1),i=1,nhapr)
      write(88,fmym80) icy,icm,(TNtr(ihapt(i),4,2),i=1,nhapr)
      end if ! qopt <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>    
 400  continue ! month loop

 500  continue ! year loop

 990  close(1)
      close(3)
      close(4)
      close(6)
      
      xmet = 0. 
      do kp=1,kpg                                             !!!!delete
        iha = ihls(kp,1)                                      !!!!delete
        xmet(iha) = xmet(iha) + dac(kp,1) + dac(kp,2)         !!!!delete
      end do                                                  !!!!delete
           
      write(13,'(/,a6,500f8.0)') 'dacsum', (xmet(ihapt(i)),i=1,nhapr)     ! delete=====================
      
      
      
      do iu=11,26 
        close(iu)
      end do

 8081 format(i4,i2.2,25f8.2,/,(6x,25f8.2))

      stop
      end

