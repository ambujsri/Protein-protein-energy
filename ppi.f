c	Program to calculate the inter-molecular energy in 
c	protein-protein complexes
c	file: ~gromiha/backup/work/protein-protein/ppi.f
c	date: 07 Feb 2008

	implicit real*8(a-h,o-z)

	character*3 filein*99
	character*15 fileout*15

	character*5 pdbfile1
	character*5 pdbfile2
	character*9 pdbout

	dimension natom(9999),filein(9999),ame1(9999),no1(9999),
     $    x1(3,9999),nat1(9999),atname1(9999),ch1(9999)
	dimension ame2(9999),no2(9999),
     $    x2(3,9999),nat2(9999),atname2(9999),ch2(9999)
	dimension atnameaa(9999),ameaa(9999),charge(9999),vdwrad(9999),
     $    vdweps(9999)
	dimension charge1(9999),vdwrad1(9999),vdweps1(9999),
     $    charge2(9999),vdwrad2(9999),vdweps2(9999)
	dimension energy(9999,9999)

	open(unit=1,file="pdb_list.dat")
	open(unit=31,file="aa20_vdrch.dat")
	open(unit=21,file="test1.out")
	open(unit=41,file="energy.out")

	print *, " give number of data, ndat "
	read(*,*)ndat

c	*** Reading the list of PDB codes ***
	do 10 i=1,ndat
	read(1,2)natom(i),filein(i)
	write(*,2)natom(i),filein(i)
 2	format(8x,i6,1x,a5)
 10	continue

c	*** Reading the charge and vdW radii for all 167 atoms in 20 AA ***
	do 501 i=1,167
	read(31,31)atnameaa(i),ameaa(i),charge(i),vdwrad(i),vdweps(i)
c	write(21,31)atnameaa(i),ameaa(i),charge(i),vdwrad(i),vdweps(i)
 31	format(12x,a4,1x,a3,8x,f9.4,3x,f7.4,3x,f7.4)
 501	continue

	do 20 i1=1,ndat,2
 	nata=natom(i1)
	nbta=natom(i1+1)
c	na=500	
	pdbfile1=filein(i1)
	pdbfile2=filein(i1+1)
c	write(*,4)nata,pdbfile1,nbta,pdbfile2
	write(21,4)nata,pdbfile1,nbta,pdbfile2
 4	format(2(1x,i6,2x,a5))

c	*** Reading PDB files ***
	open(unit=2,file=pdbfile1)

	read(2,40)(nat1(j),atname1(j),AME1(j),ch1(j),NO1(j),
     $             X1(1,j),X1(2,j),X1(3,j),j=1,nata)
c	write(21,40)(nat1(j),atname1(j),AME1(j),ch1(j),NO1(j),
c     $             X1(1,j),X1(2,j),X1(3,j),j=1,nata)
 40     FORMAT(6x,i5,1x,a4,1x,A3,1x,a1,1x,I3,4X,3F8.3)

	open(unit=4,file=pdbfile2)

	read(4,40)(nat2(j),atname2(j),AME2(j),ch2(j),NO2(j),
     $             X2(1,j),X2(2,j),X2(3,j),j=1,nbta)
c	write(21,40)(nat2(j),atname2(j),AME2(j),ch2(j),NO2(j),
c     $             X2(1,j),X2(2,j),X2(3,j),j=1,nbta)

c	pdbout=pdbfile1//".dist"
c	write(*,3)pdbfile1,pdbout
c 3	format(1x,a5,1x,9a)
c
c	open(unit=3,file=pdbout)

c	*** Assigning charges and vdW parameters ***

	do 411 i=1,nata
	 do 412 j=1,167
	  if((atname1(i).eq.atnameaa(j)).and.
     $                      (ame1(i).eq.ameaa(j))) then
	   charge1(i)=charge(j)
	   vdwrad1(i)=vdwrad(j)
	   vdweps1(i)=vdweps(j)
	go to 411
c	write(21,31)atname1(i),ame1(i),charge1(i),vdwrad1(i),vdweps1(i)
	   else
	   charge1(i)=0.0
	   vdwrad1(i)=0.0
	   vdweps1(i)=0.0
	   endif
 412	continue
 411	continue

	do 421 i=1,nbta
	 do 422 j=1,167
	  if((atname2(i).eq.atnameaa(j)).and.
     $                      (ame2(i).eq.ameaa(j))) then
	   charge2(i)=charge(j)
	   vdwrad2(i)=vdwrad(j)
	   vdweps2(i)=vdweps(j)
	go to 421
c	write(21,31)atname2(i),ame2(i),charge2(i),vdwrad2(i),vdweps2(i)
	   else
	   charge2(i)=0.0
	   vdwrad2(i)=0.0
	   vdweps2(i)=0.0
	   endif
 422	continue
 421	continue
	
c	do 701 i=1,999
c	do 701 j=1,999
c 701	energy(i,j)=0.0

c	nres1=1
c	nres2=1


           do 801 i=1,nata
		energyres=0.0
c		if(no1(i).ne.no1(i+1)) nres1=nres1+1
              do 901 j=1,nbta
c		if(no2(j).ne.no2(j+1)) nres2=nres2+1
c		if(no2(j).ne.no2(j-1)) then
c		  write(41,*)nres1,nres2,energy(nres1,nres2)
c		endif
                 eps=vdweps1(i)*vdweps2(j)
                 vdwr=vdwrad1(i)+vdwrad2(j)
                 chrg=charge1(i)*charge2(j)

                 vdwr2=vdwr*vdwr
                 vdwr6=vdwr2*vdwr2*vdwr2
                 vdwr12=vdwr6*vdwr6
                 epssqr=dsqrt(eps)

                 aec12ab=epssqr*vdwr12
                 aec6ab=2.0d0*epssqr*vdwr6

                 rxij=x2(1,j)-x1(1,i)
                 ryij=x2(2,j)-x1(2,i)
                 rzij=x2(3,j)-x1(3,i)
                 rijs=rxij**2+ryij**2+rzij**2
                 rs=1.0d0/rijs
                 rij2=rs
                 rij6=rs*rs*rs
                 rij12=rij6*rij6
                 v1=aec12ab*rij12-aec6ab*rij6
                 v2=chrg*rij2             
                 v12=v1+v2
		 energyres=energyres+v12
c                 energy(nres1,nres2)=energy(nres1,nres2)+v12
 901          continue
	write(41,51)atname1(i),AME1(i),NO1(i),energyres
 51	format(1x,a4,1x,A3,1x,I3,1x,f13.7)
 801       continue

c	close(unit=41)
	close(unit=4)
	close(unit=2)

 20	continue

	stop
	end
 
