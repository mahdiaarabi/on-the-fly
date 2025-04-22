!Updates in this version:
!1. Can now declare symmetric coordinates to skip negative values in that dimension
!
!Readme notes:
!1. The gin.log file initially should contain the energy at the equilibrium geometry,
!with whatever the convention specified within the crystal code.
!Currently, the gin.log can be replaced by something like the following line*:
! SCF Done:  E(RTPSSh) =  -4444.4
!
!*Note that this file will change on execution,
!so it must be re-newed each time before submitting the code
!
!2. To run the file in Frontera after step 1, 
!enter the following 2 lines (from the same folder):
!ifort crystal.f90
!sbatch submit.sh
!
!3. The output file will contain cordinates and energies: "outputLattices.txt"
!
!
!Comments: 
!For now, the change in energy is used to determine if Gaussian calculations have completed. 
!Please make sure every symmetric cordinate q is set as symmetric(q)=1 inside the code. 


program main
implicit none
real, dimension (:,:),allocatable:: gog,nogo,arrae
real, dimension(:),allocatable:: qq,qq2, mcPoint
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real,dimension(:),allocatable:: pot, pot2
real::Ecut0,vbar,start,finish,time,step,dumdu,now,addIt, vbarTMP, ein, einOld
integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij
integer:: msize,p,s,chki,ibar2,z,dump,grth,leth,el,t,y,xx,same,energyChecked

character (len=1024) :: text
character (len=1024) :: text2,text3
character (len=1024) :: xxx
character:: reader
!character(len=1024)::H1XPs, H1YPs

real:: H1X,H1Y,nscal
real:: H1Xp,H1Yp,oldH1xp,oldH1yp
real:: scal

!writ(*,*) "checking"
scal=1.0
d=2  !!--------> This should be the dimensions of normal mode coordinates 

allocate(symmetric(d)) 
symmetric=0
!symmetric(3) =1  !!--------> Set = 1 for dimensions that are symmetric. NOTE: Fix this
                 !!Crystal skips the negative direction of these coordinates


Open(UNIT = 12, file="outputlattices.txt")      !Output lattices go here
write(12,*) 
Close(12)

dump=0

msize=9999999  ! max basis size (for allocating memory).
Ecut0=-2638.09     !!! cutoff energy for accepting lattices. (6000 cm-1)
nscal=1.0      ! lattice spacing."

! Emin=-2638.11738478 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

open(unit=180, file="radiux.txt") !        <<<<<<<<<<<<<<<
                                  
Allocate(rad(d))                  !                      |
                                  !                      |
!defined below for simplicial growth. otherwise, this array can be externally replaced by the desired replacement of growth type. In that case, comment out this part.

ct=d*2+1

rad=0
write(180,*)rad
do iii=1,d
    rad=0
    do i=1,d
      if (i==iii) then
      rad(iii)=1
      end if
    end do
    write(180,*)rad
end do

do iii=1,d
  if(symmetric(iii)==0) then
    rad=0
    do i=1,d
      if (i==iii) then
      rad(iii)=1
      end if
    end do
      write(180,*)rad*(-1)
  else
    ct=ct-1
  end if    
end do
rewind(180)


Allocate(arrae(ct,d))     !Now reading the above file
do i=1,ct
read(180,*) (arrae(i,j),j=1,d)
end do

allocate (gog(d,msize))              
allocate (nogo(d,msize))               
allocate (le(msize))
allocate (ri(msize))
allocate (leNo(msize))
allocate (riNo(msize))
allocate(qq(d))
allocate(qq2(d))
allocate(mcPoint(d))
allocate(pot(msize))
allocate(pot2(msize))

ci=1
cf=1
ibar=1
ibar2=1
le=-1
ri=-1
leNo=-1
riNo=-1
step=0
same=2
energyChecked=2
text="abc"

!                       Defining the first seed point                   
do i=1,d                                                                                
gog(i,1)=0      
nogo(i,1)=99999                                                                         
end do                                                                                  

ein=-99999   !Insert equilibrium energy here.
vbarTMP=ein

 H1X= (0.996292908516088)
 H1Y= (1.2820849950736215)

!do while (step<=0) !-->this loop is for modified terminations
ci=1
cf=ibar

do while (cf>=ci)
!!
if (ibar>msize/8) then !To avoid unwanted symmetry breaks
write(4,*) "killing code due to the msize/5 allocation", ibar
write(*,*) "killing code due to the msize/5 allocation", ibar
dump=1
exit
end if
!!

dumdu=0

do i=ci,cf      



  do chki=1,d
  qq(chki)=gog(chki,i)
  qq2(chki)=gog(chki,i)
  end do

  do t=2,ct
      do chki=1,d
        qq(chki)=gog(chki,i)+arrae(t,chki)*1.0
      end do

      do chki=1,d
       if (qq(chki) .eq. qq2(chki)) then
         mcPoint(chki)=qq(chki)
       else
         mcPoint(chki)=qq(chki) !!! Choose either qq(chki) or (qq(chki)+qq2(chki))/2.0
       end if
      end do


      call sort(mcPoint,d,nogo,leNo,riNo,msize,i,unique,iii)
      if(unique.ne.0) then
        dumdu=unique
        call sort(qq,d,gog,le,ri,msize,i,unique,iii)    
        if(unique.ne.0) then

          !!!!!!!!!!! Transformation Matrix goes here!!!!!!!!!!!!
          

          H1X= 0.996292908516088 + (0.08091089939954188*qq(1))
          H1Y= 1.2820849950736215 + (0.06647614720468225*qq(2))

          open(unit=24, file="gin.gjf")
          write(24,*) "%NProcShared=56"
          write(24,*) "%mem=64GB"
          write(24,*) "#P GFInput IOP(6/7=3) 5D"
          write(24,*) "Gen TPSSh Nosym scf(maxcycle=500)"
          write(24,*)
          write(24,*) "DFT_2D"
          write(24,*)
          write(24,*) "+1, 1"
          write(24,*)"Fe 	0.0		 0.0		0.0"
          WRITE(24,*)"H" , "    ", H1X, "     ", H1Y, "       ", 0.00000
          write(24,*)"H 	1.4861807219207144	 0.581403851120372	0.0"
          write(24,*)"H 	0.996292908516088	-1.149132700664257	0.0"
          write(24,*)"H 	-1.2797790665095365	-2.620657920322011	1.0664087784621854"
          write(24,*)"H 	-2.768416240838588	-1.527111079009803	0.0"
          write(24,*)"H 	-1.2797790665095365	-2.620657920322011	-1.0664087784621854"
          write(24,*)"H 	-0.5175962886994335	0.6405308818173794	3.0904580087596476"
          write(24,*)"H 	0.08994225990261373	-1.384068462242426	2.8107581920712565"
          write(24,*)"H 	1.528737951421494	0.192460529568697	2.694330638467"
          write(24,*)"H 	-2.9111353767303982	1.2164807772857158	0.0"
          write(24,*)"H 	-1.6147440955681684	2.5206850612186917	-1.0617140996494856"
          write(24,*)"H 	-1.6147440955681684	2.5206850612186917	1.0617140996494856"
          write(24,*)"H 	1.528737951421494	0.192460529568697	-2.694330638467"
          write(24,*)"H 	0.08994225990261373	-1.384068462242426	-2.8107581920712565"
          write(24,*)"H 	-0.5175962886994335	0.6405308818173794	-3.0904580087596476"
          write(24,*)"P 	0.25447691775962866	-0.13007503855615057	-2.1969981742522284"
          write(24,*)"P 	0.25447691775962866	-0.13007503855615057	2.1969981742522284"
          write(24,*)"P 	-1.5542045006274163	1.597393919014446	0.0"
          write(24,*)"P 	-1.3692629015444848	-1.708066275716068	0.0"
          write(24,*)" "
          write(24,*)"Fe P H 0"
          write(24,*)"def2TZVPP"
          write(24,*)"****"
          write(24,*)" "
          write(24,*)" "
          write(24,*)
          close(24)


          same=2
          energyChecked=2
!          write(*,*) "abc"

          do zzz=1,90000000

            open(unit=14, file="gw.log")
                      
            do
!              write(*,*) "reading line by line"

              read (14,"(a)",iostat=irr) text ! Read line into character variable
              
              IF (irr/=0) then
               write(*,*) "End of file. exiting gin.log"
               close(14)
               EXIT
              end if

              if(index(text, "Charge =  0 Multiplicity = 1").ne.0) then
!                write(*,*) "reading"
!                READ (14, *) reader, H1XP, H1YP
               
                 read (14,"(a)",iostat=irr2) text2



!                 character::string = reader, H1XPs, H1YPs
!                 integer:: n, iarray(100), H1XP,
!                 n = count(transfer(text2, 'a', len(text2)) == ",")
!                 read(text2,*) reader, H1XPs, H1YPs 
!                 write(*,*) "read succesfully, now writing"
!                 read(H1XPs,*) H1XP
!                 read(H1YPs,*) H1YP

                 oldH1xp=h1xp
                 oldH1yp=h1yp

                 if(irr2/=0) then
                   write(*,*) "crash due to reading blank averted"
                   close(14)
                   exit
                 end if
                 read(text2,*) reader, H1XP, H1YP
!                 write(*,*) "read text succesfully, now writing"  

                write(*,*) "checking H"
                if((Abs(h1xp-h1x)<0.001) .and. (Abs(h1yp-h1y)<0.001)) then
                  same=1
               
                end if

!                  Write(*,*) "same? ", same

              end if
             
              if(index(text, "SCF Done").ne.0 .and. same==1) then
                energyChecked=1
                write(*,*) "gaussian energy output found"
                first= index(text, "=")
                last= index(text,"A.U.")
                xxx= text(first+1:last-1)
                close(14)
                exit
              else
!                WRITE(*,*) "New but incomplete gin.log. Waiting to complete"
                energyChecked=2
              end if

              if(index(text, "Problem with the distance matrix").ne.0 .and. same==1) then
                energyChecked=3
                vbarTMP=Ecut0+100000
                exit
              end if

            end do

              if (energyChecked==3) then
                write(*,*) "Problem with the distance matrix"
                exit
              end if

              if (energyChecked==1) then
                write(*,*) "saving energy"
                read(xxx,*) ein                   ! Convert xxx into integer
                write(*,*) "energy saved"
              end if
            close(14)


            if(energyChecked==2) then
              call sleep(1)
            else
              vbarTMP=ein
              write(*,*) ein
!              close(14)    ! Should probably add this line
              exit
            end if
          end do

 
          if(vbarTMP<=Ecut0 .and. step==0 .and. ein.ne.9999999) then
            if (unique==1) then
            !unique is 0 if the vector is already in the list
            !Else it is 1 or -1. Depends on if it goes to left or right branch at the end.
              call addTo(ibar,ri,gog,iii,qq,d,pot,vbartmp,msize)
              !adds vector to gog, and adds vector's label to ri(iii)
              
              open(12, file="outputlattices.txt", status="old", position="append", action="write")
              write(12,*) qq, vbartmp
              close(12)
              
            else if (unique==-1) then
              call addTo(ibar,le,gog,iii,qq,d,pot,vbartmp,msize)
              
              open(12, file="outputlattices.txt", status="old", position="append", action="write")
              write(12,*) qq, vbartmp
              close(12)
              
            end if
          else if(vbarTMP>Ecut0 .and. step==0)then
            if (dumdu==1) then
              call addTo2(ibar2,riNo,nogo,iii,mcPoint,d,pot2,vbartmp,msize)
            else if (dumdu==-1) then
              call addTo2(ibar2,leNo,nogo,iii,mcPoint,d,pot2,vbartmp,msize)
            end if
          end if
          
          
          
        end if
      end if
      
    if (dump==1) then   
      exit
    end if
  
  end do
  
  if (dump==1) then   
    exit
  end if
  
end do

ci=cf+1                                                                                 !
cf=ibar 
call cpu_time(now)
!write(4,*) ibar, now-start
!write(4,*) "ibar,ci,cf: ", ibar, ci, cf

!!if (step>=1) then  !steps where you want to terminate after one iteration
!!exit
!!end if

end do                                                                                  !
!_______________________________________________!

!write(4,*) ecut0, "   ", ibar ,"   ",now-start

!do i=1,ibar
!  write(4,*) (gog(chki,i),chki=1,d)
!end do
!do i=1,ibar2
!  write(4,*) (nogo(chki,i),chki=1,d)
!end do


call cpu_time(finish)
time=finish-start
!write(*,*) time

deallocate (gog)              
deallocate (nogo)               
deallocate (le)
deallocate (ri)
deallocate (leNo)
deallocate (riNo)
deallocate(qq)
deallocate(qq2)
deallocate(mcPoint)
deallocate(pot)
deallocate(pot2)
deallocate(rad)                    !                   |
deallocate(arrae)

rewind(180)

write(*,*) "Code completed"

!end do
end program


function grth(d,qq,i,iii,msize,gog)                     !returns 1 if greater, else 0
integer:: d,i,iii,a,grth,z,xx,b,msize
real, dimension(d,msize):: gog
real, dimension(d):: qq
xx=0
do z=1,d
  if(qq(z)>gog(z,iii))then
    xx=1
    exit
  else if(qq(z)<gog(z,iii))then
    exit
  end if
end do
grth=xx
return
end function


function leth(d,qq,i,iii,msize,gog)                     !returns 1 if lesser, else 0
integer:: d,i,iii,a,leth,z,xx,b,msize
real, dimension(d,msize):: gog
real, dimension(d):: qq
xx=0
do z=1,d
  if(qq(z)<gog(z,iii))then
    xx=1
    exit
  else if(qq(z)>gog(z,iii))then
    exit
  end if
end do
leth=xx
return
end function


subroutine sort(qq,d,gogo,le,ri,msize,i,unique,iii) !unique is (0)/(1)/(-1) if the vector is
!(already in list)/(on right side of lowest element in the binary-sort tree)/(on left side)
integer:: i,d,msize,iii,unique,grth,leth
real,dimension(d)::q,qq
real,dimension(d,msize)::gogo
integer, dimension(msize):: le,ri
         unique=0                             
         iii=1                                                                  
         do                                                                             
             if ((grth(d,qq,i,iii,msize,gogo))==1) then !If (the vector being sorted is greater 
                                                      !than the iii-th in the binary tree) then
                 if (ri(iii)>0) then  !If (there is an element to the right of iii)
                   iii = ri(iii) !Then, the next time when this the loop repeats, 
                     ! compare the vector with the element on the right side of iii
                 else                                                   
                   unique=1   !Else report to the main code that the vector is a new element. (iii and unique are the outputs of this section)
                 exit
                 end if                                                 
             else if (leth(d,qq,i,iii,msize,gogo)==1) then      
                 if (le(iii)>0) then                    
                   iii = le(iii)                                
                 else
                   unique=-1
                 exit
                 end if                                                 
             else                                                               
                 exit    !Unique remains 0 if the vector is same as iii
             end if                                                             
         end do                                                                 
end subroutine


subroutine addTo(iba,ri,gogo,iii,qq,d,pot,vbartmp,msize)
integer:: iba,iii,el,d
integer, dimension(msize):: ri
real,dimension(d,msize)::gogo
real,dimension(d)::qq
real,dimension(msize)::pot
iba=iba+1       
ri(iii)=iba                                             !
do el=1,d
  gogo(el,iba)=qq(el)
end do
pot(iba)=vbartmp
end subroutine

subroutine addTo2(iba,ri,gogo,iii,qq,d,pot2,dimen,msize)
integer:: iba,iii,el,d
integer, dimension(msize):: ri
real,dimension(d,msize)::gogo
real,dimension(d)::qq
real,dimension(msize)::pot2
iba=iba+1       
ri(iii)=iba                                             !
do el=1,d
  gogo(el,iba)=qq(el)
end do
pot2(iba)=dimen
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
