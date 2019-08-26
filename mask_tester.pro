Function mask_tester, DataHeader, DataData, outputRes,ScanProto, SENSEFactor
; Receives as input the header of the cpx, the data, outputRes, the slices, the information of the Scan protocol, the SENSE factor, the mask and the time mask
;SensesitivitiesTime=SYSTIME(/SECONDS ) ; Time spent until now
;slice=0
dir='H:\CORONARIAS\'
;ScanProto=ReadProtocol(dir+'CORONARIAS+EDEMA.txt') ; Same for the Scan protocol
;maskKspace=Generate_KSpace_mask(maskKspace)
nlines = (file_lines(DataHeader))(0) ; Number of lines of the file
openr, 1, DataHeader ; Open the header, index 1
header=strarr(21, nlines) ; String array of 21 characters
rstring=''
for i=0, 300 do begin
  readf,1,rstring
  if strmatch(rstring, '*NOI*') and (strmatch(rstring, '*#*') ne 1) then begin ; Look for lines beginning by NOI but not beginning by #
     header=long((STRSPLIT(rstring, /extract))(1:*)) ; Split the lines into separate strings when finding a 1
  endif
endfor
close,/all
openr, 1, DataHeader
openr, 2, DataData
VectorPosition=0LL
KzPos=0
KyPos=0
CoilChan=0
for i=0., nlines-1 do begin
    readf,1,rstring
    if strmatch(rstring, '*#*') ne 1 then begin ; Look for lines not starting by #
       if strmatch(rstring, '*STD*') then begin;or strmatch(rstring, '*REJ*') ; Look for lines starting by STD
            header=long((STRSPLIT(rstring, /extract))(1:*))
            header(8)=temporary(header(8))*SENSEFactor(0) ; k-space coordinate in 1st encoding direction (Y), affected by the SENSE factor in Y
            header(9)=temporary(header(9))*SENSEFactor(1) ; k-space coordinate in 2nd encoding direction (Z), affected by the SENSE factor in Z

            KyPos=[temporary(KyPos),header(8)] ; Vector with the Y positions in the k-space
            KzPos=[temporary(KzPos),header(9)] ; Vector with the Z positions in the k-space
            CoilChan=[temporary(CoilChan),header(5)] ; Vector with the synchronization channel numbers
            VectorPosition=[temporary(VectorPosition),header(19)] ; Vector position determined by the data vector offset (in bytes)
       endif
     endif
endfor
    close,/all
    VectorPosition=temporary(VectorPosition(1:*)) ; Saves the number of line where the file is being read
    KyPos=temporary(KyPos(1:*)) ; Vector with the Y dimensions in the k-space
    ind=where(KyPos lt 0, count) ; Indexes where the Y positions in the k-space are 0, they are also counted
    if count gt 0 then KyPos(ind)=outputRes(0)+temporary(KyPos(ind)) ; If there are positions which are 0 in the Y dimension in the k-space, they are assigned the value of the number of rows of the final image
    ; Same in Z dimension of the k-space
    KzPos=temporary(KzPos(1:*))
    ind=where(KzPos lt 0, count)
    if count gt 0 then KzPos(ind)=outputRes(1)+temporary(KzPos(ind))
    CoilChan=temporary(CoilChan(1:*)) ; Vector with the synchronization channel numbers
    CoilChan=CoilChan-min(CoilChan) ; The first synchronization channel number is set to 0 in case it is not 0
    nlines=n_elements(VectorPosition) ; Number of lines of the file
    ;n_channels=max(coilchan)-min(coilchan)+1 ; The number of channels used is the difference between the maximum and minimum channel synchronization numbers
    mask=fltarr(outputRes(0),outputRes(1),outputRes(2),1) ; Floating array of the same dimensions than "data" but only in one slice
    for i=0L, (nlines-1) do begin
       mask(KyPos(i), KzPos(i), CoilChan(i),*)=temporary(mask(KyPos(i), KzPos(i), CoilChan(i),*))+1.0 ; Mask is set to 1
    endfor
close, /all
return, mask
end