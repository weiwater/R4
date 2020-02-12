c Updated 2/6/2018
c Modified from rdoutp1.for   
c to simulate planned diversion from reservoir
c  
      subroutine rdoutp(file0)

      include 'RO_DIM_3.MAX'

      include 'GenParm.VAR'
      include 'OVprpt.VAR'
      include 'lkprpt2.VAR'

      character*8 ch8(MXPL)
      character*50 ofile,bl50,file0

      bl50 = '                                                  '
      
c     read location IDs for daily flow output

      read(1,'(i8)') nplk
      write(*,'(a,i5)') ' NPLK =',nplk
      if(nplk.gt.MXPL) then
        write(*,'(a,i5)')   ' Specified nplk = ',nplk
        goto 990
      end if
 
      if(nplk.gt.0) then
        write(6,'(/,a,i3,a)') 'Daily flow at ',nplk,
     &                        ' locations will be printed'
        write(6,'(a)')   '    Name  Number'
        do i=1,nplk
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do j=1,nlinks
            if(ch8(i).eq.lkname(j)) then
              iplk(i,1) = j
              iplk(i,2) = 1     ! location is link
              goto 10
            end if
          end do
          do j=1,nnodes
            if(ch8(i).eq.ndname(j)) then
              iplk(i,1) = j
              iplk(i,2) = 2     ! locaiton is node
              goto 10            
            end if
          end do

          write(*,'(a,a8,a)')' Link or node ',ch8(i),' is not available'
          ierror = .true.
 10       continue
          write(6,'(a8,2i4)') ch8(i),iplk(i,1),iplk(i,2)
        end do

c       read the file name to print flow at iplk(i, ) and open the file
c
        call OpnOF2(1,61,ofile,'* Daily flow output (cfs)',file0)
c        write(61,'(a)') 
        write(61,'(i8)') nplk
        write(61,'(8x,100a8)') (ch8(i), i=1,nplk)
        write(61,'(8x,100i8)') (iplk(i,1),i=1,nplk)
        write(61,'(8x,100i8)') (iplk(i,2),i=1,nplk)
        if(qopt) then
          call OpnOF2(1,62,ofile,
     &               '* TDS concentration in daily flow (mg/L)',file0)
          write(62,'(i8)') nplk
          write(62,'(8x,100a8)') (ch8(i), i=1,nplk)
          write(62,'(8x,100i8)') (iplk(i,1),i=1,nplk)
          write(62,'(8x,100i8)') (iplk(i,2),i=1,nplk)
          call OpnOF2(1,63,ofile,
     &            '* Nitrate concentration in daily flow (mg/L)',file0)
          write(63,'(i8)') nplk
          write(63,'(8x,100a8)') (ch8(i), i=1,nplk)
          write(63,'(8x,100i8)') (iplk(i,1),i=1,nplk)
          write(63,'(8x,100i8)') (iplk(i,2),i=1,nplk)
        end if
      end if !nplk

c     read location IDs for monthly flow output file

      read(1,'(i8)') npnd
      write(*,'(a,i5)') ' NPND =',npnd

      if(npnd.gt.MXPL) then
        write(*,'(a,i5)')   ' Specified npnd = ',npnd
        goto 990
      end if

      if(npnd.gt.0) then
        write(6,'(/,a,i3,a)') 'Monthly flow at', npnd,
     &                        ' locations will be printed'
        write(6,'(a)')   'NodeName  Number'
        do i=1,npnd
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do in=1,nnodes
            if(ch8(i).eq.ndname(in)) then
              ipnd(i) = in
              goto 20
            end if
          end do
          write(*,'(a,a8,a)')' Node ',ch8(i),' is not available'
          ierror = .true.
 20       continue
          write(6,'(a8,2i4)') ch8(i),ipnd(i)
        end do

c       read file name to print flow at ipnd(i) and open file
c       print at each link using
c
        call OpnOF2(1,64,ofile,
     &              '* Total monthly flow passing the nodes',file0)
        write(64,'(a)') '* Units (ac-ft)'
        write(64,'(i6)') npnd
        write(64,'(6x,100a8)') (ch8(i), i=1,npnd)
        write(64,'(1x)')
        if(qopt) then
          call OpnOF2(1,65,ofile,
     &               '* Total monthly TDS passing the nodes',file0)
c          write(65,'(a)') '* Units (pounds)'
          write(65,'(a)') '* Units (mg/L)'
          write(65,'(i6)') npnd
          write(65,'(6x,100a8)') (ch8(i), i=1,npnd)
          call OpnOF2(1,66,ofile,
     &               '* Total monthly TIN passing the nodes',file0)
c          write(66,'(a)') '* Units (pounds)'
          write(66,'(a)') '* Units (mg/L)'
          write(66,'(i6)') npnd
          write(66,'(6x,100a8)') (ch8(i), i=1,npnd)
        end if
      end if !npnd

c     read reach names to print stream reach deep percolation data
      read(1,'(i8)') nprch
      write(*,'(a,i5)') ' NPRCH =',nprch

      if (nprch.gt.MXRH) then
        write(*,'(a)') ' nprch > MXRH'
        write(*,'(a)') ' Stop execution'
        stop
      end if
  
      if(nprch.gt.0) then
        write(6,'(/,a,i3,a)') 'Monthly percolation to',nprch,
     &                        ' reaches will be printed'
        write(6,'(a)')   'Number Reach Name'

        do i=1,nprch
          read(1,'(a8)') rchname(i)
          call ctrim8r(rchname(i))
          write(6,'(i4,2x,a8)') i,rchname(i)
        end do

        call OpnOF2(1,67,ofile,
     &      '* Total monthly percolation to the river reaches',file0)
        write(67,'(a)') '* Units|(ac-ft)'
        write(67,'(i5)') nprch
        write(67,'(6x,(50a8))') (rchname(i), i=1,nprch)
c        write(67,'(1x)')

        if(qopt) then
          call OpnOF2(1,68,ofile,
     &               '* Total TDS monthly percolation',file0)
c          write(68,'(a)') '* Units|(pounds)'
          write(68,'(a)') '* Units|(mg/L)'
          write(68,'(i5)') nprch
          write(68,'(6x,(50a8))') (rchname(i), i=1,nprch)
          call OpnOF2(1,69,ofile,
     &               '* Total TIN monthly percolation',file0)
c          write(69,'(a)') '* Units|(pounds)'
          write(69,'(a)') '* Units|(mg/L)'
          write(69,'(i5)') nprch
          write(69,'(6x,(50a8))') (rchname(i), i=1,nprch)
        end if
      end if

c     read recharge basins IDs for monthly percolation and evaporation
      read(1,'(i8)') nqls
      write(*,'(a,i5)') ' NQLS =',nqls

      if(nqls.gt.0) then
        write(6,'(/,a,i3,a)') 'Monthly percolation to',nqls,
     &              ' recharge basins will be printed.'
        write(6,'(a)')   'LinkName  Number'
        do i=1,nqls
          read(1,'(a8)') ch8(i)
          call ctrim8r(ch8(i))      ! trim the name and convert to upper case
          do il=1,nlinks
            if(ch8(i).eq.lkname(il)) then
              ipql(i) = il
              goto 30
            end if
          end do
          write(*,'(a,a8,a)')' Link ',ch8(i),' is not available'
          ierror = .true.
          xqevp(i) = 0.
          cxqloss(i,1) = 0.
          cxqloss(i,2) = 0.
 30       xqloss(i) = 0.
          write(6,'(a8,2i4)') ch8(i),ipql(i)
        end do
        
        call OpnOF2(1,81,ofile,
     &             '* Total monthly percolation',file0)
        write(81,'(a)') '* Units|(ac-ft)'
        write(81,'(i8)') nqls
        write(81,'(8x,100a8)') (ch8(i),i=1,nqls)        
        write(81,'(8x,100i8)') (ipql(i),i=1,nqls)
        if(qopt) then
          call OpnOF2(1,82,ofile,
     &               '* Total TDS in monthly percolation',file0)
c          write(82,'(a)') '* Units|(pounds)'
          write(82,'(a)') '* Units|(mg/L)'
          write(82,'(i8)') nqls
          write(82,'(8x,100a8)') (ch8(i),i=1,nqls)        
          write(82,'(8x,100i8)') (ipql(i),i=1,nqls)
          
          call OpnOF2(1,83,ofile,
     &               '* Total nitrate in monthly percolation',file0)
c          write(83,'(a)') '* Units|(pounds)'
          write(83,'(a)') '* Units|(mg/L)'
          write(83,'(i8)') nqls
          write(83,'(8x,100a8)') (ch8(i),i=1,nqls)        
          write(83,'(8x,100i8)') (ipql(i),i=1,nqls)
        end if

        call OpnOF2(1,84,ofile,
     &             '* Total monthly evaporation',file0)
        write(84,'(a)') '* Units|(ac-ft)'
        write(84,'(i8)') nqls
        write(84,'(8x,100a8)') (ch8(i),i=1,nqls)        
        write(84,'(8x,100i8)') (ipql(i),i=1,nqls)        
        
        call OpnOF2(1,85,ofile,
     &             '* End-of-day Storage',file0)
        write(85,'(a)') '* Units|(ac-ft)'
        write(85,'(i8)') nqls
        write(85,'(8x,100a8)') (ch8(i),i=1,nqls)        
        write(85,'(8x,100i8)') (ipql(i),i=1,nqls)        

      end if ! nqls
      

      if(npdd.gt.0) then
        read(1,'(a50)',end=40) ofile                              ! 7/21/14 *.ADV file.
        if (ofile.ne.bl50) then                            
          open(91,file=ofile,status='UNKNOWN')             
          write(*,'(1x,a40,a)') ofile,' is open'
          write(91,'(a)') '* Echo print of predefined diversion'
          write(91,'(a)') '* Units|(ac-ft)'                 
          write(91,'(i5)') npdd
        end if  
      end if    
      
      goto 40
 
 990  write(*,'(/,a,i5)') ' Maximum number = ',MXPL
      write(*,'(a)') ' Stop execution'
      stop
 40   return      
      end


