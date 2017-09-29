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

! This routine reads and setup link/node system

      subroutine LNSetup()

      include 'RO_DIM_4.MAX'

c---- general model control properties
      include 'GenParm.VAR'

c---- link/node properties
      include 'LKprpt.var'

c---- ha property
      include 'HAprpt.var'

c---- reservoir property
      include 'RVprpt.VAR'

c---- diversion properties
      include 'DVprpt.VAR'

c---- table properties
      include 'TBprpt.VAR'

      character*8  ch8(100), bl8,chend
      real         rtemp(MXHA)


      bl8 = '        '
      chend = 'END     '       
      velocity = 1.0
c
c     define node properties
c
      write(*,'(/,a)') ' Reading node data'
      j=0
 11   j=j+1
 12   read(3,'(a8,2i8)') ndname(j),nha(j)
      if(ndname(j).eq.bl8) goto 12
      call ctrim8r(ndname(j))
      if(ldebug) write(*,'(i10,2x,a8)') j,ndname(j)
      if(ndname(j).eq.chend) then
        nnodes = j-1
        write(*,'(a,i5)') ' Number of Nodes = ',nnodes
        goto 20
      end if

      if(nha(j).le.0) go to 11
      if(nha(j).gt.MXHN) then
        write(*,'(a,i4,a10,a)') '  Number of loads at node',j,ndname(j),
     &                         ' is greater than MXHN'
        stop
      end if
      read(3,1044) (ch8(k),k=1,nha(j))
      do k=1,nha(j)
        call ctrim8r(ch8(k))
      end do
      ! read multiplier
      ! this is for small mountain watershed
      read(3,'(10f8.0)') (rtemp(k),k=1,nha(j))

c      write(6,'(/,a)') ' Checking source of water for nodes ========'
      do k=1,nha(j)
        do kk=1,nhat2
          if(ch8(k).eq.haname(kk)) then
            nodha(j,k) = kk
            rmult(j,k) = rtemp(k)
            go to 13
          end if
        end do
        write(6,'(a,i2,1x,a8,a,i4,1x,a8,a)') ' HA name ',k,ch8(k),
     &           ' for node #',j,ndname(j),' is not available'
        ierror = .true.
        jerror = jerror + 1
 13     continue
      end do
      go to 11

 20   close(3)

 1001 format(10f8.0)
 1044 format(10a8)
 1045 format(' Match up of nodes and hydrologic areas, node= ',a8,
     1    ' node index=',i4,'  ha= ', a8,' column in runoff file=', i8)
 1047 format(3i9,2a8)

c
c     read link data
c
      write(*,'(/,a)') ' Reading link data'
      write(6,'(/,a)') ' Reading link data ==========='
      nres = 0   ! reservoir counter
      ndiv = 0   ! diversion counter
      j = 0      ! link counter
 21   j = j+1
 22   read(2,'(3a8,5i8)') lkname(j),ch8(1),ch8(2),convtype(j)
      if(lkname(j).eq.bl8) goto 22
      call ctrim8r(lkname(j))
      if(ldebug) write(*,'(i10,2x,a8)') j,lkname(j)
      if(lkname(j).eq.chend) then
        write(*,'(a)') ' Finished reading link file'
        goto 30
      end if

      call ctrim8r(ch8(1))
      call ctrim8r(ch8(2))
      do in=1,nnodes
        if(ch8(1).eq.ndname(in)) then
          usn(j) = in
          goto 23
        end if
      end do
      write(6,'(a,a,a,a,a)')  ' Upstream node   ',ch8(1), ' for link ',
     &          lkname(j),' is not defined.'
      ierror = .true.
      jerror = jerror + 1
 23   do in=1,nnodes
        if(ch8(2).eq.ndname(in)) then
          dsn(j) = in
          goto 24
        end if
      end do
      write(6,'(a,a,a,a,a)')  ' Downstream node ',ch8(2), ' for link ',
     &          lkname(j),' is not defined.'
      ierror = .true.
      jerror = jerror + 1
 24   continue

c
c    convtype = 1   conveyance/percolation
c                   Link ID
c                   Manning     - open channel with Manning's equaiton
c                   Blank       - just move water to downstream
c                   Table name  - flow-width table name
c             = 2   conveyance only, like enclosed pipe
c             = 3   storage and conveyance
c             = 4   diversion
c             = 5   dummy link

      goto (210,220,230,245,250) convtype(j)
        write(6,'(a,i3,a)') ' Link type ',convtype(j),
     &                      ' is not available.  Stop execution.'
        write(6,'(5x,3a8,5i8)') lkname(j),ch8(1),ch8(2),convtype(j)
        ierror = .true.
        jerror = jerror + 1
        goto 250

 210  read(2,'(a8,10f8.0)') ch8(1),(rtemp(k),k=1,9)
      if(rtemp(8).ne.0) rxn(j) = rtemp(8)
      call ctrim8r(ch8(1))
      if(ch8(1).eq.bl8) then
        slop(j)    = 0.
        mannn(j)   = 0.
        b(j)       = 0.
        z(j)       = 0.
        length(j)  = 0.
        dpb(j)     = 0.
        dpz(j)     = 0.
      elseif(ch8(1).eq.'MANNING') then
        slop(j)    = rtemp(1)
        mannn(j)   = rtemp(2)
        b(j)       = rtemp(3)
        z(j)       = rtemp(4)
        length(j)  = rtemp(5)
        lkrch(j)   = rtemp(6)
        dpb(j)     = rvrch(lkrch(j),1)
        dpz(j)     = rvrch(lkrch(j),2)
        if(dpz(j) .eq. 0.) dpz(j) = dpb(j)
      else  ! this could be flow-width rating option channel
        do ich=1,nch
          if(ch8(1).eq.chname(ich)) then
            convtype(j) = 6
            ptrch(j) = ich
            length(j) = rtemp(5)
c           dpb(j)    = rtemp(6)
            lkrch(j)  = rtemp(6)
            dpb(j)    = rvrch(lkrch(j),1)
            dpz(j)    = rvrch(lkrch(j),2)
            if(dpz(j) .eq. 0.) dpz(j) = dpb(j)
            goto 215
          end if
        end do
        write(*,'(a,1x,a8,a)') ' CH name ',ch8(1),
     &           ' is not available'
        ierror = .true.
        jerror = jerror + 1
      end if

 215  rxn(j) = rxn(j)*length(j)/(velocity*86400.)

      if(rtemp(9).gt.0. .and. rtemp(9).le.MXRH) linkrch(j) = rtemp(9)

      ! unit conversion
      dpb(j) = dpb(j)/86400.               ! (ft/day) to (ft/sec)
      dpz(j) = dpz(j)/86400.
      go to 250
 220  read(2,'(8x,10f8.0)') slop(j),mannn(j),b(j)
      go to 250

c
c     Link is a reservoir
c
 230  nres=nres+1
      ptrrv(j)=nres
      read(2,'(a40)') rvname(nres)
      write(*,'(16x,i3,1x,a)') nres,rvname(nres)
      read(2,'(i8)') npts(nres)

      slimit(nres,1) = 0.                  ! dead storage
      slimit(nres,2) = 0.                  ! storage at spillway starts
      slimit(nres,3) = 0.                  ! maximum storage
      do kk=1,npts(nres)
        read(2,1001) rvelev(nres,kk),rvarea(nres,kk),rvstor(nres,kk),
     &               rvout(nres,kk,1), rvout(nres,kk,2),rvpr(nres,kk)
        if(rvout(nres,kk,1).eq.0. .and. rvout(nres,kk,2).eq.0)
     &               slimit(nres,1) = rvstor(nres,kk)
        if(                          rvout(nres,kk,2).eq.0)
     &               slimit(nres,2) = rvstor(nres,kk)
      end do  
      alimit(nres)   = rvarea(nres,npts(nres))
      slimit(nres,3) = rvstor(nres,npts(nres))

      ! check for storage increase
      do kk=1,npts(nres)-1
        if(rvstor(nres,kk).ge.rvstor(nres,kk+1)) then
          write(*,'(a)') ' Problem with reservoir storage'
          write(*,'(a,i3)') ' Reservoir ',nres
          k=nres
          write(*,'(6f8.2)') rvelev(k,kk),rvarea(k,kk),rvstor(k,kk)
          write(*,'(6f8.2)') rvelev(k,kk+1),rvarea(k,kk+1),
     &                       rvstor(k,kk+1)
          stop
        end if
      end do

c     read downstream link for each outflow
      read(2,'(24x,2a8)') rdlink(nres,1), rdlink(nres,2)
      call ctrim8r(rdlink(nres,1))
      call ctrim8r(rdlink(nres,2))
      write(*,'(10x,2a8)')  rdlink(nres,1), rdlink(nres,2) 
      lkatrv(nres) = j

      read(2,1001) (evwts(nres,kk),kk=1,newts)
      read(2,'(10f8.0)') (rtemp(kk),kk=1,10)
      rxn(j) = rtemp(9)
      if(rtemp(10).gt.0..and.rtemp(10).le.MXRH) linkrch(j) = rtemp(10)
      krvtype(nres) = rtemp(2)
      ffflow(nres) = rtemp(4)  ! this is the first flush to be wasted on ffmode day

      slimit(nres,4) = rtemp(3)

      if(rtemp(3).eq.0.) then ! operating storage is not specified
        if(krvtype(nres).eq.1)  then
          slimit(nres,4) = slimit(nres,2)
        elseif (krvtype(nres).eq.2)  then
          slimit(nres,4) = slimit(nres,1)
        else
          slimit(nres,4) = (slimit(nres,1) + slimit(nres,2)) / 2.
        end if
      end if
      ! krvtype() = 1 --- conservation basin
      ! krvtype() = 2 --- flood control basin
      ! krvtype() = 3 --- multipurpose basin
      ! for conservation basins, outlet gates are set
      !   for best conservation
      ! for flood control basins, outlet gates should be set
      !   to drain reservoir within 24 hours
      ! for multipurpose basins, all gates should be set for conservation
      !   when a storm is coming, all water should be just released
      !   to the first downstream link.

      do kk=1,npts(nres)
        if(rvpr(nres,kk).eq.0.) rvpr(nres,kk) = rtemp(1)
      end do

      ! setup drain flow tables
      k = nres
      call EASSET(k)
      goto 250
c
c     diversion link
c
 245  ndiv = ndiv+1
      ptrdv(j)=ndiv
      read(2,'(a40)') dvname(ndiv)
      read(2,'(i8)') ndpts(ndiv)
      do kk=1,ndpts(ndiv)
        read(2,1001) dflow(ndiv,kk),divout(ndiv,kk)
      end do
c
c     read downstream links
c       linkdv(ndiv,1) is the main link
c       linkdv(ndiv,2) is the diversion link
c
      read(2,'(2a8)') ddlink(ndiv,2),ddlink(ndiv,1)
      call ctrim8r(ddlink(ndiv,1))
      call ctrim8r(ddlink(ndiv,2))
      lkatdv(ndiv) = j
       ! *** diversion can not be made to upstream link

 250  go to 21
 30   nlinks = j-1
      write(*,'(a,i5)') ' Number of Links = ',nlinks
      close(2)


c
c     assign link indices to diversion destination links
c
      write(*,'(a)') ' Assigning diversion links'
      do k=1,ndiv
        do i2=1,2
          linkdv(k,i2)=0
          do jk=lkatdv(k)+1,nlinks          ! diversion should not be made to a upstream link
           if(ddlink(k,i2).eq.lkname(jk)) then
             linkdv(k,i2)=jk
             go to 40
           end if
          end do
 40       continue
        end do

        if(linkdv(k,1).eq.0 .or. linkdv(k,2).eq.0) then
          write (*,'(a,i3,a)') ' Outflow links for diversion',k,
     &                            ' are not properly specified'
          ierror = .true.
          jerror = jerror + 1
        end if
      end do
c
c     assign link indices to reservoir outlet links
c
      write(*,'(a)') ' Assigning reservoir outlet links'
      do k=1,nres
        do i2=1,2
          linkrv(k,i2)=0
          if(rdlink(k,i2).ne.bl8) then
            do jk=lkatrv(k)+1,nlinks
             if(rdlink(k,i2).eq.lkname(jk)) then
               linkrv(k,i2)=jk
               go to 41
             end if
            end do
 41         continue
          end if
        end do

        if( (rdlink(k,1).ne.bl8 .and. linkrv(k,1).eq.0) .or.
     &      (rdlink(k,2).ne.bl8 .and. linkrv(k,2).eq.0) ) then
          write (*,'(a,i3,a)') ' Outflow links for reservoir',k,
     &                         ' is not properly specified'
          ierror = .true.
          jerror = jerror + 1
        end if
      end do

c     ==================================================================
c
c     Configure Link-Node system
c
      write(*,'(a)') ' Configuring link-node system'
      ! top link does not have incoming link
c      nincl(1) = 0
      do j=1,nlinks
        kk=0
        do jj=1,nlinks
          if(usn(j).eq.dsn(jj)) then
            kk=kk+1
            incoml(usn(j),kk)=jj
          end if
        end do
        nincl(usn(j))=kk
        if(ldebug) write(*,'(2(i5,1x,a8),i5, 10(i5,1x,a8))') 
     &                   j, lkname(j),usn(j), ndname(usn(j)),
     &                   nincl(usn(j)),
     &               (incoml(usn(j),k),lkname(incoml(usn(j),k)),
     &                    k=1,nincl(usn(j)))
      end do

c     ==================================================================
c
c     write link-node configuration to main output file
c
      write(6,'(//,a,/)') '   *** Link - Node Configuration ***'
      write(6,'(a43,a21,a)')
     &        'Link Name    Connecting Nodes        Type  ',
     &        'Node Name   NINCL NHA',
     &        '  Incomming Links'
      write(6,'(13x,a26,27x,a)')
     &        'Upstream     Downstream   ',
     &        'HA Column # and Names'

      do j=1,nlinks
        ! upstream node
        iun = usn(j)
        if(nincl(iun).eq.0) then
          write(6,'(/,43x,i4,1x,a8,2i4,2x,a)')
     &             iun,ndname(iun),nincl(iun),nha(iun),'Headwater'
        else
          write(6,'(43x,i4,1x,a8,2i4,2x,5(i4,1x,a8))')
     &             iun,ndname(iun),nincl(iun),nha(iun),
     &             (incoml(iun,kk),
     &              lkname(incoml(iun,kk)),kk=1,nincl(iun))
        end if
        if(nha(iun).gt.0) then
          write(6,'(66x,10(i4,1x,a8))')
     &             (nodha(iun,kk),haname(nodha(iun,kk)),kk=1,nha(iun))
          write(6,'(66x,10(5x,f8.3))')
     &             (rmult(iun,kk),kk=1,nha(iun))
        end if

        ! write link information
        write(6,'(3(i4,1x,a8),i2)') j,lkname(j),
     &          usn(j),ndname(usn(j)),dsn(j),ndname(dsn(j)),
     &          convtype(j)
        if(convtype(j).eq.3) then
          ires = ptrrv(j)
          write(6,'(5x,a40)') rvname(ires)
          write(6,'(5x,a,2(i4,1x,a8))') 'Downlink',
     &            linkrv(ires,1),lkname(linkrv(ires,1)),
     &            linkrv(ires,2),lkname(linkrv(ires,2))
        elseif(convtype(j).eq.4) then
          idiv = ptrdv(j)
          write(6,'(5x,a40)') dvname(idiv)
          write(6,'(5x,a,2(i4,1x,a8))') 'Downlink',
     &            linkdv(idiv,2),lkname(linkdv(idiv,2)),
     &            linkdv(idiv,1),lkname(linkdv(idiv,1))
        end if
      end do ! nlinks

      ! last node
      idn = dsn(nlinks)
      write(6,'(43x,i4,1x,a8,2i4,2x,5(i4,1x,a8))')
     &         idn,ndname(idn),nincl(idn),nha(idn),
     &         (incoml(idn,kk),lkname(incoml(idn,kk)),kk=1,nincl(idn))
      if(nha(idn).gt.0) then
        write(6,'(66x,10(i4,1x,a8))')
     &           (nodha(idn,kk),haname(nodha(idn,kk)),kk=1,nha(idn))
        write(6,'(66x,10(5x,f8.3))')
     &           (rmult(idn,kk),kk=1,nha(idn))
      end if
      write(6,'(/,a,/)') 'End of Link-Node Configuration'
      write(6,'(a,i5)') 'Number of Links, nlinks = ',nlinks
      write(6,'(a,i5)') 'Number of Nodes, nnodes = ',nnodes


      write(6,'(/,a)') 'Diversions in the system'
      do k=1,ndiv
        write(6,'(2i5,2x,a,2f10.1)') k,lkatdv(k),dvname(k)
      end do
      write(6,'(/,a)') 'Reservoirs in the system'
      do k=1,nres
        write(6,'(2i5,a10,a,5f10.1)') k,lkatrv(k),lkname(lkatrv(k)),
     &               rvname(k),alimit(k),slimit(k,1),slimit(k,2),
     &               slimit(k,3)
      end do

      return
      end
