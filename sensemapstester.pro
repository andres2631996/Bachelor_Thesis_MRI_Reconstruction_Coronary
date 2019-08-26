Function SenseMapsTester, Coils, QBodyCoil, RefProto, ScanProto, outputRes, FOVCorrection, slice, mask, ExportMask
;Fractional SENSE correction
QBThreshold=0.05
RealImageFOV=[ScanProto.FOV(1),ScanProto.FOV(2),ScanProto.FOV(0)] ; Image FOV (YxZ and X is represented as slices)
; LINE DOWN: SPECIFIC FOR FRACTIONAL SENSE FACTORS
ScanProto.FOV(1:2)=temporary(ScanProto.FOV(1:2))*ceil(ScanProto.SENSE)/ScanProto.SENSE ; FOV modified by the SENSE factors in different directions
TotalSlices  = round(ScanProto.FOV(0)/ScanProto.VOXELSIZE(0)) ; Slices are obtained by dividing the FOV along X over the voxel size along X
print, 'SENSE Maps Generation'
  ;if keyword_set(QBThreshold) ne 1 then QBThreshold=1e-2 ; If the Q-body threshold is not specified, it is 10^(-2)
  if keyword_set(FOVCorrection) then ScanProto.FOV=temporary(ScanProto.FOV)/FOVCorrection ; If the FOVCorrection is not specified, the FOV is divided by that correction
  dim=size(Coils,/dimension) ; Dimensions of the array containing the information about the coils
  if ScanProto.FOV(2) gt RefProto.FOV(2) then begin
  ; Compare the Z component of the FOV between the scan and ref protocols
      FOVInc=ceil((ScanProto.FOV(2)-RefProto.FOV(2))/RefProto.VOXELSIZE(2)/2.)>0 ; Quotient between the difference in both FOVs and the half of the Z component of the voxel size
      inter=complexarr(dim(0), dim(1)+2*FOVInc,dim(2)) ; Complex array of dimensions the ones of Coils, but as first dimension, the first dimension of Coils + twice the previous calculated quotient
      inter(*,FOVInc:FOVInc+dim(1)-1,*)=QBodyCoil ; The complex array above for the first dimension between the quotient and the quotient + the first dimension of Coils equals the information from the Q-body coil
      QBodyCoil=inter ; The information saved in inter is now saved in Q-body coil
      inter=complexarr(dim(0), dim(1)+2*FOVInc,dim(2),dim(3)) ; Inter is now a complex array of the dimension of Coils, including its third dimension, but in the first dimension, twice the quotient between FOV differences is also added
      inter(*,FOVInc:FOVInc+dim(1)-1,*,*)=Coils ; The complex array between the quotient and the quotient + the first dimension of Coils equals Coils
      Coils=inter ; Coils is saved in the previous complex array
      inter=0
      dim(1)=dim(1)+2*FOVInc ; The Z dimension is increased by twice the quotient of the difference between FOVs
      RefProto.FOV(2)=dim(1)*RefProto.VOXELSIZE(2) ; The new FOV Z component is calculated, too, multiplying by the Z component of the voxel size
  endif
  if ScanProto.FOV(1) gt RefProto.FOV(1) then begin
    ; Compare the Y component of the FOV between the scan and ref protocols. Same process than for Z, except the FOV quotient
      FOVInc=ceil((ScanProto.FOV(1)-RefProto.FOV(1))/RefProto.ReconVoxel/2.)>0 ; Same quotient than for Z, but now dividing the difference between the size of the reconstruction voxel
      inter=complexarr(dim(0)+2*FOVInc, dim(1), dim(2))
      inter(FOVInc:FOVInc+dim(0)-1,*,*)=QBodyCoil
      QBodyCoil=inter
      inter=complexarr(dim(0)+2*FOVInc, dim(1),dim(2),dim(3))
      inter(FOVInc:FOVInc+dim(0)-1,*,*,*)=Coils
      Coils=inter
      inter=0
      dim(0)=dim(0)+2*FOVInc
      RefProto.FOV(0)=dim(0)*RefProto.ReconVoxel
  endif
  final_res=[ScanProto.FOV(1)/outputRes(0),ScanProto.FOV(2)/outputRes(1),ScanProto.FOV(0)/TotalSlices] ; Voxel size on the scan protocol (in mm)
;  OriProtoRes=([RefProto.ReconVoxel, RefProto.VOXELSIZE(2), RefProto.ReconVoxel]) ; Voxel size on the Ref protocol (mm)
;  offCenter=refproto.offcenter-scanproto.offcenter ; Difference between Off centers in Ref and Scan protocols
;  offCenter=[offCenter(1),offCenter(2),offCenter(0)]+OriProtoRes ; The difference between off centers is added now the voxel size on the ref protocol
;  relative_shift=[0,0,0];
  ;print, round(offCenter/final_res) ; Quotient between the difference between off centers and the voxel size in the scan protocol (mm)
  SMDim=round([RefProto.FOV(1),RefProto.FOV(2),RefProto.FOV(0)]/final_res) ; Quotient between the FOV in the ref protocol and the voxel size in the scan protocol (number of voxels)
;ind=where(abs(QBodyCoil) gt max(abs(QBodyCoil))*QBThreshold) ; Find where the absolute value of the Q-body coil is higher than the maximum value multiplied by a threshold
;  if keyword_set(mask) then begin ; Look if mask has been defined
;    mask2=fltarr(outputRes(0),outputRes(1),TotalSlices) ; Float array of the dimensions of the matrix size and the number of slices
;    ExportMask=fltarr(outputRes(0),outputRes(1),TotalSlices) ; Float array of the dimensions of the matrix size and the number of slices
;    mask2(ind)=1. ; If the Q-body coil values are higher than the maximum scaled by a threshold, they are set to 1
;    for i=0, TotalSlices-1 do begin ; Go slice by slice
;      a=MORPH_CLOSE(reform(mask2(*,*,i)),replicate(1,50.,50.))+mask2(*,*,i) ; Mask 2 is expressed as a matrix for each slice and closed with a kernel of 50x50 full of ones added to the original mask2
;      ind=where(a gt 0, count) ; Find where the closed mask + the original one is positive and substitute it by ones ---> BINARY MASK
;      if count gt 0 then begin
;      	a(ind)=1
;      	mask2(*,*,i)=a
;      	a=shift(float(dilate(reform(mask2(*,*,i)), replicate(1,10.,10.))),0,1) ; Dilate the mask with a kernel of 10x10 full of ones and then shift it 1 unit in Z direction
;      	a=convol(a,FilterGen(15,15),/EDGE_t); Convolve the dilated and shifted mask with a Gaussian filter of 15 samples and width 15
;      	b=shift(float(dilate(reform(mask2(*,*,i)), replicate(1,5.,5.))),0,1) ; Dilate the mask with a kernel of 5x5 full of ones and then shift it 1 unit in Z direction
;      	b=(convol(b,FilterGen(15,15),/EDGE_t)) ; Convolve the dilated and shifted mask with a Gaussian filter of 15 samples and width 15
;      	mask2(*,*,i)=a ; The dilated mask with a 10x10 kernel and shifted is saved in mask2
;      	ExportMask(*,*,i)=b ; The dilated mask with a 5x5 kernel and shifted is saved in ExportMask
;      endif else begin
;      	ExportMask(*,*,i)=make_array(outputRes(0),outputRes(1),value=1)
;      endelse
;    endfor
;  endif
SenseMaps=fltarr(outputRes(0),outputRes(1),outputRes(3),outputRes(2)) ; Complex array of dimensions: Matrix Size X Number of slices X Number of coils
for i=0, outputRes(2)-1 do begin ; For each coil
      am=complexarr(SMDim(0),SMDim(1),SMDim(2)) ; Complex array with dimensions the ones of the quotient between the FOV in the ref protocol and the voxel size in the scan protocol
      coil=Coils(*,*,*,i) ; Keeps the information for each coil
      coil=shift(fft(shift(temporary(coil),-Dim(0)/2,-Dim(1)/2,-Dim(2)/2), 1),Dim(0)/2,Dim(1)/2,Dim(2)/2)
      ; Shift along the dimensions of Coils, perform the IFFT of coil, and then unshift
      am(round((SMDim(0)-Dim(0))/2):round((SMDim(0)+Dim(0))/2)-1,$
        round((SMDim(1)-Dim(1))/2):round((SMDim(1)+Dim(1))/2)-1,$
        round((SMDim(2)-Dim(2))/2):round((SMDim(2)+Dim(2))/2)-1)=coil
      ; a is a vector whose values in the intermediate indexes are assigned to coil
      am=shift(fft(shift(temporary(am),-SMDim(0)/2,-SMDim(1)/2,-SMDim(2)/2),-1),SMDim(0)/2,SMDim(1)/2,SMDim(2)/2)
      ; a is now shifted to the dimensions obtained from the quotient between the FOV and the voxel size, its IFFT is obtained and then it is unshifted
      ;am=shift(am,relative_shift(0),relative_shift(1),relative_shift(2)) ; Shift "a" with respect to the quotient between the difference between off centers and the voxel size in the scan protocol (mm)
      am=temporary(am(round((SMDim(0)-outputRes(0))/2):round((SMDim(0)+outputRes(0))/2)-1, $
          round((SMDim(1)-outputRes(1))/2):round((SMDim(1)+outputRes(1))/2)-1, $
          round((SMDim(2)-outputRes(3))/2)+slice:round((SMDim(2)+outputRes(3))/2)-1+slice))
      ; a just keeps the values of the intermediate indexes
      ;b=complexarr(outputRes(0),outputRes(1),TotalSlices) ; Complex array of dimensions: Matrix Size X Number of slices
      if keyword_set(mask) then begin ; See if "mask" is defined
          SenseMaps(*,*,*,i)=am/QBodyCoil ; b is the normalized version of a by the Q-body coil information
          SenseMaps(*,*,*,i)=temporary(SenseMaps(*,*,*,i))*ExportMask;mask2 ; b is filtered by the mask that was previously dilated and then shifted with a 5x5 kernel
      endif else SenseMaps(*,*,*,i)=am/QBodyCoil ; If mask has already been defined, do the same as if it had not been defined
      ;SenseMaps(*,*,*,i)=b ; b is the sensitivity map for each coil
endfor
;mask=ExportMask
; Q-BODY COIL LINEAR INTERPOLATION
; The Q-body coil images that are obtained need to be interpolated in order to align them in the Z axis with images obtained in the final reconstructed volume
return, SenseMaps
end