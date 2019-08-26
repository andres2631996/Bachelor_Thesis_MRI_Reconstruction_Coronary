Function QBodyCoil_tester, Coils, mask, FOVCorrection, outputRes, QBodyCoil
QBThreshold=5e-2
dir='H:\CORONARIAS\'
ScanProto=ReadProtocol(dir+'CORONARIAS+EDEMA.txt')
RefProto=ReadProtocol(dir+'Ref_XL_Torso.txt')
  if keyword_set(QBThreshold) ne 1 then QBThreshold=1e-2 ; If the Q-body threshold is not specified, it is 10^(-2)
  ;     Ky          Kz            nCoil       Slices      Phases
  ;outputRes(0),outputRes(1),outputRes(2),outputRes(3),outputRes(4));
  if keyword_set(FOVCorrection) then ScanProto.FOV=ScanProto.FOV/FOVCorrection ; If the FOVCorrection is not specified, the FOV is divided by that correction
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

  final_res=[ScanProto.FOV(1)/outputRes(0),ScanProto.FOV(2)/outputRes(1),ScanProto.FOV(0)/outputRes(3)] ; Voxel size on the scan protocol (in mm)
  ;                 RL                      AP                    FH
  OriProtoRes=([RefProto.ReconVoxel, RefProto.VOXELSIZE(2), RefProto.ReconVoxel]) ; Voxel size on the Ref protocol (mm)

  offCenter=refproto.offcenter-scanproto.offcenter ; Difference between Off centers in Ref and Scan protocols
  offCenter=[offCenter(1),offCenter(2),offCenter(0)]+OriProtoRes ; The difference between off centers is added now the voxel size on the ref protocol

  relative_shift=round(offCenter/final_res) ; Quotient between the difference between off centers and the voxel size in the scan protocol (mm)
  SMDim=round([RefProto.FOV(1),RefProto.FOV(2),RefProto.FOV(0)]/final_res) ; Quotient between the FOV in the ref protocol and the voxel size in the scan protocol (number of voxels)

  if keyword_set(FOVCorrection) then ScanProto.FOV=ScanProto.FOV*FOVCorrection
  ; See if FOVCorrection has been defined. In that case, the FOV in the scan protocol is scaled by the FOV correction
  ;a=complexarr(SMDim(0),SMDim(1),SMDim(2)/2) ; Complex array with dimensions the ones of the quotient between the FOV in the ref protocol and the voxel size in the scan protocol
  QBodyCoil=shift(fft(shift(QBodyCoil,-Dim(0)/2,-Dim(1)/2,-Dim(2)/2), 1),Dim(0)/2,Dim(1)/2,Dim(2)/2)
  ; The parameters of the Q-body coil are shifted around half of the dimension of Coils, and then the IFFT is computed
  ; Then, the IFFT is unshifted
;  a(round((SMDim(0)-outputRes(0))/2):round((SMDim(0)+outputRes(0))/2)-1,$
;    round((SMDim(1)-outputRes(1))/2):round((SMDim(1)+outputRes(1))/2)-1,$
;    round((SMDim(2)-outputRes(3))/2):round((SMDim(2)+outputRes(3))/2)-1)=QBodyCoil
  ; The central indexes of "a" are equal to the information transmitted by the Q-body coil, so the information of the Q-body coil is saved there
  QBodyCoil=fft(shift(QBodyCoil,-SMDim(0)/2,-SMDim(1)/2,-SMDim(2)/2),-1) ; The Q-body coil information is shifted, performed a FFT, and unshifted
  ;a=shift(a,relative_shift(0),relative_shift(1),relative_shift(2)) ; Shift "a" with respect to the quotient between the difference between off centers and the voxel size in the scan protocol (mm)

;  QBodyCoil=a(round((SMDim(0)-outputRes(0))/2):round((SMDim(0)+outputRes(0))/2)-1,$
;              round((SMDim(1)-outputRes(1))/2):round((SMDim(1)+outputRes(1))/2)-1,$
;              round((SMDim(2)-outputRes(3))/2):round((SMDim(2)+outputRes(3))/2)-1)
  ; The Q-body coil is set to the values in the intermediate indexes
;  ind=where(abs(QBodyCoil) gt max(abs(QBodyCoil))*QBThreshold) ; Find where the absolute value of the Q-body coil is higher than the maximum value multiplied by a threshold
;  if keyword_set(mask) then begin ; Look if mask has been defined
;    mask2=fltarr(outputRes(0),outputRes(1),outputRes(3)) ; Float array of the dimensions of the matrix size and the number of slices
;    ExportMask=fltarr(outputRes(0),outputRes(1),outputRes(3)) ; Float array of the dimensions of the matrix size and the number of slices
;    mask2(ind)=1. ; If the Q-body coil values are higher than the maximum scaled by a threshold, they are set to 1
;    for i=0, outputRes(3)-1 do begin ; Go slice by slice
;      a=MORPH_CLOSE(reform(mask2(*,*,i)),replicate(1,50.,50.))+mask2(*,*,i) ; Mask 2 is expressed as a matrix for each slice and closed with a kernel of 50x50 full of ones added to the original mask2
;      ind=where(a gt 0, count) ; Find where the closed mask + the original one is positive and substitute it by ones ---> BINARY MASK
;      a(ind)=1
;      mask2(*,*,i)=a
;      a=shift(float(dilate(reform(mask2(*,*,i)), replicate(1,10.,10.))),0,1) ; Dilate the mask with a kernel of 10x10 full of ones and then shift it 1 unit in Z direction
;      a=convol(a,FilterGen(15,15),/EDGE_t); Convolve the dilated and shifted mask with a Gaussian filter of 15 samples and width 15
;      b=shift(float(dilate(reform(mask2(*,*,i)), replicate(1,5.,5.))),0,1) ; Dilate the mask with a kernel of 5x5 full of ones and then shift it 1 unit in Z direction
;      b=(convol(b,FilterGen(15,15),/EDGE_t)) ; Convolve the dilated and shifted mask with a Gaussian filter of 15 samples and width 15
;      mask2(*,*,i)=a ; The dilated mask with a 10x10 kernel and shifted is saved in mask2
;      ExportMask(*,*,i)=b ; The dilated mask with a 5x5 kernel and shifted is saved in ExportMask
;    end
;  endif
;  SenseMaps=complexarr(outputRes(0), outputRes(1),outputRes(3),outputRes(2)) ; Complex array of dimensions: Matrix Size X Number of slices X Number of coils
;  for i=0, outputRes(2)-1 do begin ; For each coil
;      am=complexarr(SMDim(0),SMDim(1),SMDim(2)*2) ; Complex array with dimensions the ones of the quotient between the FOV in the ref protocol and the voxel size in the scan protocol
;      coil=Coils(*,*,*,i) ; Keeps the information for each coil
;      coil=shift(fft(shift(coil,-Dim(0)/2,-Dim(1)/2,-Dim(2)/2), 1),Dim(0)/2,Dim(1)/2,Dim(2)/2)
;      ; Shift along the dimensions of Coils, perform the IFFT of coil, and then unshift
;      am(round((SMDim(0)-Dim(0))/2):round((SMDim(0)+Dim(0))/2)-1,$
;        round((SMDim(1)-Dim(1))/2):round((SMDim(1)+Dim(1))/2)-1,$
;        20+round((SMDim(2)-Dim(2))/2):round((SMDim(2)+Dim(2))/2)+19)=coil
;      ; a is a vector whose values in the intermediate indexes are assigned to coil
;      am=shift(fft(shift(am,-SMDim(0)/2,-SMDim(1)/2,-SMDim(2)/2),-1),SMDim(0)/2,SMDim(1)/2,SMDim(2)/2)
;      ; a is now shifted to the dimensions obtained from the quotient between the FOV and the voxel size, its IFFT is obtained and then it is unshifted
;      am=shift(am,relative_shift(0),relative_shift(1),relative_shift(2)) ; Shift "a" with respect to the quotient between the difference between off centers and the voxel size in the scan protocol (mm)
;      am=am(round((SMDim(0)-outputRes(0))/2):round((SMDim(0)+outputRes(0))/2)-1, $
;          round((SMDim(1)-outputRes(1))/2):round((SMDim(1)+outputRes(1))/2)-1, $
;          round((SMDim(2)-outputRes(3))/2):round((SMDim(2)+outputRes(3))/2)-1)
;      ; a just keeps the values of the intermediate indexes
;      b=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Complex array of dimensions: Matrix Size X Number of slices
;      if keyword_set(mask) then begin ; See if "mask" is defined
;          b=am/QBodyCoil ; b is the normalized version of a by the Q-body coil information
;          b=b*ExportMask;mask2 ; b is filtered by the mask that was previously dilated and then shifted with a 5x5 kernel
;      endif else b=am/QBodyCoil ; If mask has already been defined, do the same as if it had not been defined
;      SenseMaps(*,*,*,i)=b ; b is the sensitivity map for each coil
;  endfor
  return, QBodyCoil
  end