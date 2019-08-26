Function data_tester ,DataHeader, DataData, outputRes, slice,ScanProto, SENSEFactor, mask ;TFEFactor=TFEFactor,
; Receives as input the header of the cpx, the data, outputRes, the slices, the information of the Scan protocol, the SENSE factor, the mask and the time mask
;SensesitivitiesTime=SYSTIME(/SECONDS ) ; Time spent until now
dir='H:\CORONARIAS\'
ScanProto=ReadProtocol(temporary(dir)+'CORONARIAS+EDEMA.txt') ; Same for the Scan protocol
;maskKspace=Generate_KSpace_mask(maskKspace)
print,'Reading Raw Data'
nlines = (file_lines(DataHeader))(0) ; Number of lines of the file
openr, 1, DataHeader ; Open the header, index 1
header=strarr(21, nlines) ; String array of 21 characters
rstring=''


MinCoil=0 ; Minimum coil
for i=0, 300 do begin
  readf,1,rstring
  if strmatch(rstring, '*NOI*') and (strmatch(rstring, '*#*') ne 1) then begin ; Look for lines beginning by NOI but not beginning by #
     header=long((STRSPLIT(rstring, /extract))(1:*)) ; Split the lines into separate strings when finding a 1
     ; RR=[RR,header(17)]
     MinCoil=[temporary(MinCoil),header(5)] ; Information for the minimum coil is in the 5th character of each line
  endif
endfor

MinCoil=min(temporary(MinCoil(1:*)))  ; Minimum value of the vector with the coils with the minimum values


close,/all
openr, 1, DataHeader
openr, 2, DataData
VectorPosition=0LL
KzPos=0
KyPos=0
CoilChan=0
SENSEFactor=fltarr(2)
SENSEFactor=ScanProto.SENSE
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
       endif else begin
       if strmatch(rstring, '*NOI*') then begin
            header=long((STRSPLIT(rstring, /extract))(1:*))
            dummy=complexarr(round(header(18)/8.)) ; Complex array of length in bits
            point_lun, 2, header(19)
            readu, 2, dummy ; Dummy receives the information from data
            if keyword_set(noise) ne 1 then begin ; Determines if "noise" has been defined
               noise=complexarr(ncoils, round(header(18)/8)) ; Double-precision, complex array of size the number of coils x length of "dummy" vector
               noise(header(5)-MinCoil, *)=dummy ; Noise matrix is filled with the noises from all the coils used
            endif else noise(header(5)-MinCoil, *)=dummy ; If "noise" already exists, it receives the same value than it receives in the case it does not exist
            endif
            if strmatch(rstring, '*kx_range*') then begin
                xdim=fix((strsplit(rstring, /extract))(7))-fix((strsplit(rstring, /extract))(6))+1 ; X dimensions in the k-space
            endif
            if strmatch(rstring, '*ky_range*') then ydim=fix((strsplit(rstring, /extract))(7))*SENSEFactor(0); Y dimensions, affected by the first component of the SENSE factor
            if strmatch(rstring, '*kz_range*') then zdim=fix((strsplit(rstring, /extract))(7))*SENSEFactor(1); Z dimensions, affected by the second component of the SENSE factor
            if strmatch(rstring, '*number_of_locations*') then nlocs=fix((strsplit(rstring, /extract))(6)) ; Keeps an integer value of the number of locations
            if strmatch(rstring, '*0  number of coil channels*') then if nlocs eq 1 then  ncoils=fix((strsplit(rstring, /extract))(9)) else ncoils=fix((strsplit(rstring, /extract))(9)+1)
            ; If the number of locations is 1, the number of coil channels is the one stated in the file, otherwise it is increased by 1
            if strmatch(rstring, '*0  kx_oversample_factor*') then FH_OversamplingFactor=float(((strsplit(rstring, /extract)))(6)) ; Oversampling factor
            if strmatch(rstring, '*0  ky_oversample_factor*') then ScanProto.FOV(1)=ScanProto.FOV(1)*float(((strsplit(rstring, /extract)))(6)) ; Y component of the FOV
            if strmatch(rstring, '*0  kz_oversample_factor*') then ScanProto.FOV(2)=ScanProto.FOV(2)*float(((strsplit(rstring, /extract)))(6)) ; Z component of the FOV
          endelse
       endif
    endfor
    close,/all
    VectorPosition=temporary(VectorPosition(1:*)) ; Saves the number of line where the file is being read

    outputRes(0)= ceil((2.*ydim)/(SenseFactor(0)))*(SenseFactor(0)) ; Y component of the matrix size (double of the dimensions previously calculated)
    outputRes(1)= ceil((2.*zdim)/(SenseFactor(1)))*(SenseFactor(1)) ; Z component of the matrix size (double of the dimensions previously calculated)
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
    n_channels=max(coilchan)-min(coilchan)+1 ; The number of channels used is the difference between the maximum and minimum channel synchronization numbers
    DataVector=complexarr(xdim) ; Complex vector with the length of the X dimensions in the k-space
    Data=complexarr(outputRes(0),outputRes(1),outputRes(2),outputRes(3)) ; Final output with the matrix size, the number of cardiac phases, the number of slices and the number of coils used. It is saved in a complex array
    ;mask=fltarr(outputRes(0),outputRes(1),outputRes(2),1) ; Floating array of the same dimensions than "data" but only in one slice
    filter=exp(-((findgen(xdim)-xdim/2)/xdim)^2) ; Gaussian filter along dimensions of X. It is a 1D vector with length the number of slices

    openr, 2, DataData ; Open the data file, index=2
    filter=fltarr(xdim)
    filter=exp(-(findgen(xdim)-xdim/2.)^2./20000)
    filter(0:4)=0.0 ; First 5 elements of the filter, which will be used in the edges, will be 0
    filter(xdim-6:*)=0.0 ; Last 5 elements of the filter, which will be used in the edges, will be 0, too
    DataPoint=complexarr(outputRes(3))
    filter=convol(filter, FilterGen1D(30,30),/EDGE_t) ; Final filter is convolution of the original 1D filter with a Gaussian filter of 30 samples and width 30
    for i=0L, (nlines-1) do begin
    aux=complexarr(xdim)
       point_lun, 2, VectorPosition(i) ; Moves the pointer to the beginning of the file
       readu, 2, DataVector ; Reads each line the data file
       DataVector=temporary(DataVector)*filter ; Filtered positions in the X dimension
       DataVector=shift(fft(shift(temporary(DataVector), -xdim/2),/inverse, /OVERWRITE ), xdim/2)
       ; Shift the filtered positions in the X dimension, compute the IFFT and undo the shifting
       DataPoint=DataVector(xdim/2-floor(outputRes(3)/2.)+slice:xdim/2-floor(outputRes(3)/2.)+outputRes(3)-1+slice)
       Data(KyPos(i), KzPos(i), CoilChan(i),*) = temporary(Data(KyPos(i), KzPos(i), CoilChan(i),*))+DataPoint ; Data in each point are actualized to the previous data calculated for each point
       ;mask(KyPos(i), KzPos(i), CoilChan(i),*) +=1.0 ; Mask is set to 1


    endfor
    close,/all
    ind=where(mask gt 0, count) ; Indexes where mask is positive
    ;print, max(mask)   ; Print the maximum value of the mask
;    if count gt 0 then TimeMask(ind)=TimeMask(ind)/mask(ind)
maskData=rebin(mask, outputRes(0),outputRes(1),outputRes(2),outputRes(3))
ind=where(maskData gt 0, count)
if count gt 0 then data(ind)=temporary(data(ind))/maskData(ind)

; FILTERING
;maskKspace=complexarr(outputRes(0), outputRes(1))
;maskKspace=mask(*,*,0)
;maskKspace=Generate_KSpace_mask(maskKspace) ; Obtain a mask to select where in the k-space we are going to work
;vectorX=(where(maskKspace(*,0) eq 0)) ; Indexes where the rows of the mask are 0
;vectorY=(where(maskKspace(0,*) eq 0)) ; Indexes where the columns of the mask are 0
;radio=float([vectorX(0),vectorY(0)])+1 ; Contains the first indexes of rows and columns with zeros. A one is added
;LowerRadio=radio-[3.,3.] ; The lower radio is 3 elements shorter than the initial radio
;HalfScan=[vectorX(n_elements(vectorX)-1),vectorY(n_elements(vectorY)-1)] ; Halfscan contains the last elements of the rows and columns of the mask which are 0
;LowerHalfScan=HalfScan+[3.,3.] ; Lower half scan is 3 elements higher than the original half scan
;maskKspace=fltarr(outputRes(0),outputRes(1)) ; Float array with dimensions the matrix size of the final image
;filter=fltarr(outputRes(0),outputRes(1)) ; Float array with dimensions the matrix size of the final image
;for j=0, outputRes(1)-1 do begin
;   for i=0, outputRes(0)-1 do begin
;        radius=((i-outputRes(0)/2.)/radio(0))^2.+((j-outputRes(1)/2.)/radio(1))^2. ; Distance to the central pixel in the k-space normalized by the first elements where the mask is 0
;        if radius le 1 then maskKspace(i,j)=1.0 ; Looks if the radius is lower than 1, assigning then the mask a 1 --> LP filter in the k-space
;        LowerRadius=((i-outputRes(0)/2.)/LowerRadio(0))^2.+((j-outputRes(1)/2.)/LowerRadio(1))^2. ; Distance to the central pixel in the k-space normalized by the first elements where the mask is 0, minus 3
;        if LowerRadius le 1 then filter(i,j)=1.0 ; Looks if the lower radius is lower than 1, assigning then the filter a 1 --> LP filter in the k-space
;    endfor
;endfor
;
;maskKspace=shift(temporary(maskKspace),outputRes(0)/2.,outputRes(1)/2.) ; Shift the mask with half of the dimensions of the matrix size
;maskKspace(radio(0):HalfScan(0),*)=0 ; The mask between the first and last indexes of rows and columns where it was 0 is now 0
;maskKspace(*,radio(1):HalfScan(1))=0
;
;filter=shift(temporary(filter),outputRes(0)/2.,outputRes(1)/2.) ; Shift the filter with half of the dimensions of the matrix size
;filter(LowerRadio(0):LowerHalfScan(0),*)=0 ; The filter between the first indexes of rows and columns where it was 0, minus 3, and the last indexes of rows and columns where it was 0, plus 3, is now 0
;filter(*,LowerRadio(1):LowerHalfScan(1))=0
;
;HomodineFilter=maskKspace ; The final mask is a LP Homodyne filter
;
;filter=shift(temporary(filter),-outputRes(0)/2.,-outputRes(1)/2.) ; The filter is unshifted
;filter=(convol(temporary(filter), FilterGen(25,19),/EDGE_t)) ; The filter is now convolved with a Gaussian filter of 25 samples and a width of 19, convolving also in the edges
;
;maskKspace=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; The mask is unshifted
;
;weight=fltarr(1,outputRes(1)) ; Vector of weights with the dimensions of Z
;weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./40000) ; Weights are normally distributed
;weight=rebin(temporary(weight), outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image
;
;a=fltarr(outputRes(0),1) ; Vector of weights with the dimensions of Y
;a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./40000) ; Weights are normally distributed
;a=rebin(a, outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image
;
;weight = temporary(weight)*a ; Weights in Y and Z directions are now multiplied between them
;
;filter=temporary(filter)*shift(HomodineFilter, -outputRes(0)/2.,-outputRes(1)/2.)*weight
;; Homodyne filter (the initial filter) is shifted with half of the dimensions of the final matrix size and multiplied by the convolved version of the filter and by the Gaussian weights
;weight=0
;
;filter=reform(shift(temporary(filter), -outputRes(0)/2.,-outputRes(1)/2.),  outputRes(0),outputRes(1),1,1) ; The filter is shifted again and adapted to the matrix size of the final image
;filter=rebin(temporary(filter),outputRes(0),outputRes(1),1,outputRes(3)) ; The filter is adapted to the dimensions of the final matrix size and the slice number
;
;for i=0, outputRes(2)-1 do begin
;    data(*,*,i,*)=temporary(data(*,*,i,*))*filter ; Filtered data
;endfor

return, data
end