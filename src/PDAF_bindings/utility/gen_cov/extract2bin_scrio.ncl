; This ncl script is to extract schism output for PDAF eof ensemble
; It will produce 6 fortran binary files & dim.dat for preparation of analysis.
;
; usage example: ncl extract2bin_scrio.ncl
; fs&fe is scribe file start & end number
; user can also specify irec loop index for extraction

system("rm -f ?.bin dim.dat") ; remove files first

fs=32
fe=396

dp=new(3,integer)
icount=0
do nof=fs,fe,2
  f1=addfile("./data/"+"out2d_"+nof+".nc","r")
  f2=addfile("./data/"+"horizontalVelX_"+nof+".nc","r")
  f3=addfile("./data/"+"horizontalVelY_"+nof+".nc","r")
  f4=addfile("./data/"+"salinity_"+nof+".nc","r")
  f5=addfile("./data/"+"temperature_"+nof+".nc","r")
  tt=f1->time
  nrec=dimsizes(tt)
  elev=f1->elevation
  hvelX=f2->horizontalVelX
  hvelY=f3->horizontalVelY
  temp=f5->temperature
  salt=f4->salinity
  ;ww=f1->vertical_velocity
  ww=hvelX ; We don't have outputs here
  ww=0.
  if (nof.eq.fs) then
     node=dimsizes(temp)
     ;printVarSummary(node)
     ;print(node)
     dp(0:1)=toint(node(1:2))
  end if
  do irec=0,0 ;nrec-1,12
   ; output records
     fbinrecwrite("z.bin",-1,todouble(elev(irec,:)))
     fbinrecwrite("w.bin",-1,todouble(ww(irec,:,:)))
     fbinrecwrite("u.bin",-1,todouble(hvelX(irec,:,:)))
     fbinrecwrite("v.bin",-1,todouble(hvelY(irec,:,:)))
     fbinrecwrite("t.bin",-1,todouble(temp(irec,:,:)))
     fbinrecwrite("s.bin",-1,todouble(salt(irec,:,:)))
     icount=icount+1
     print(icount)
     system("date")
  end do
 ; close file
   delete(f1)
   delete(f2)
   delete(f3)
   delete(f4)
   delete(f5)
   delete(tt)
   delete(nrec)
   delete(elev)
   delete(ww)
   delete(hvelX)
   delete(hvelY)
   delete(temp)
   delete(salt)
end do
dp(2)=icount
asciiwrite("dim.dat",dp)


