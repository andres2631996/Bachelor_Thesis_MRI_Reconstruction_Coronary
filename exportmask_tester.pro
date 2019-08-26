function ExportMask_tester, QBodyCoil, QBThreshold, outputRes
;TotalSlices = round(ScanProto.FOV(0)/ScanProto.VOXELSIZE(0))
QBThreshold=0.05
ind1=where(abs(QBodyCoil) gt max(abs(QBodyCoil)*QBThreshold))
; Find where the absolute value of the Q-body coil is higher than the maximum value multiplied by a threshold
  ;if keyword_set(mask) then begin ; Look if mask has been defined
    mask2=fltarr(outputRes(0),outputRes(1),outputRes(3)) ; Float array of the dimensions of the matrix size and the number of slices
    ExportMask=fltarr(outputRes(0),outputRes(1),outputRes(3)) ; Float array of the dimensions of the matrix size and the number of slices
    mask2(ind1)=1. ; If the Q-body coil values are higher than the maximum scaled by a threshold, they are set to 1
    for i=0, outputRes(3)-1 do begin ; Go slice by slice
      a=MORPH_CLOSE(reform(mask2(*,*,i)),replicate(1,50.,50.))+mask2(*,*,i) ; Mask 2 is expressed as a matrix for each slice and closed with a kernel of 50x50 full of ones added to the original mask2
      ind=where(a gt 0, count) ; Find where the closed mask + the original one is positive and substitute it by ones ---> BINARY MASK
      a(ind)=1
      mask2(*,*,i)=a
      a=shift(float(dilate(reform(mask2(*,*,i)), replicate(1,10.,10.))),0,1) ; Dilate the mask with a kernel of 10x10 full of ones and then shift it 1 unit in Z direction
      a=convol(a,FilterGen(15,15),/EDGE_t);smooth(a,8, /edge_t) ; Convolve the dilated and shifted mask with a Gaussian filter of 15 samples and width 15
      b=shift(float(dilate(reform(mask2(*,*,i)), replicate(1,5.,5.))),0,1) ; Dilate the mask with a kernel of 5x5 full of ones and then shift it 1 unit in Z direction
      b=(convol(b,FilterGen(15,15),/EDGE_t));smooth(b,4, /edge_t) ; Convolve the dilated and shifted mask with a Gaussian filter of 15 samples and width 15
      mask2(*,*,i)=a;fltarr(outputRes(0), outputRes(1))+1;a; ; The dilated mask with a 10x10 kernel and shifted is saved in mask2
      ExportMask(*,*,i)=b;fltarr(outputRes(0), outputRes(1))+1;b; ; The dilated mask with a 5x5 kernel and shifted is saved in ExportMask
    end
  ;endif
  return, ExportMask
  end