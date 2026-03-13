; This ncl script is to extract schism output for PDAF eof ensemble
; It will produce 6 fortran binary files & dim.dat for preparation of analysis.
;
; usage axample: ncl extract2bin.ncl
; fs&fe is schout file start & end number
; user can also specify irec loop index for extraction

system("rm -f ?.bin dim.dat") ; remove files first

fs=1
fe=5

dp=new(3,integer)
icount=0
do nof=fs,fe
  f1=addfile("schout_"+nof+".nc","r")
  tt=f1->time
  nrec=dimsizes(tt)
  elev=f1->elev
  ww=f1->vertical_velocity
  hvel=f1->hvel
  temp=f1->temp
  salt=f1->salt
  if (nof.eq.fs) then
     node=dimsizes(temp)
     dp(0:1)=toint(node(1:2))
  end if
  do irec=0,nrec-1
   ; output records
     fbinrecwrite("z.bin",-1,todouble(elev(irec,:)))
     fbinrecwrite("w.bin",-1,todouble(ww(irec,:,:)))
     fbinrecwrite("u.bin",-1,todouble(hvel(irec,:,:,0)))
     fbinrecwrite("v.bin",-1,todouble(hvel(irec,:,:,1)))
     fbinrecwrite("t.bin",-1,todouble(temp(irec,:,:)))
     fbinrecwrite("s.bin",-1,todouble(salt(irec,:,:)))
     icount=icount+1
     print(icount)
  end do
 ; close file
   delete(f1)
   delete(tt)
   delete(nrec)
   delete(elev)
   delete(ww)
   delete(hvel)
   delete(temp)
   delete(salt)
end do
dp(2)=icount
asciiwrite("dim.dat",dp)


