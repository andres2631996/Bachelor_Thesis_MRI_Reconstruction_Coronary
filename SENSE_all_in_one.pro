function FilterGen, nSamples, FilterWith
; Gaussian filter of dimensions nSamples x nSamples and width FilterWith
    Filter=fltarr(nSamples,nSamples)
    for i=0, nSamples-1 do begin
      for j=0, nSamples-1 do begin
          Filter(i,j)= exp(-((i-floor(nSamples/2.))^2.+(j-floor(nSamples/2.))^2.)/FilterWith)
      endfor
    endfor
    Filter=Filter/total(Filter)
    return, Filter
end


Function MedianFilter, data, XYZwidth
; Function that receives as inputs the data and the kernel dimensions
  dim=size(data, /dimen) ; Matrix size X Number of coils X Number of slices X Number of cardiac phases
  TotalLength=XYZwidth(0)*XYZwidth(1)*XYZwidth(2) ; The total length is the product of the kernel dimensions in the 3 directions
  Center=round(XYZwidth/2) ; Center of the kernel
  for k=XYZwidth(2), dim(2)-XYZwidth(2)-1 do begin ; Go through the X dimension
      for j=XYZwidth(1), dim(1)-XYZwidth(1)-1 do begin ; Go through the Z dimension
          for i=XYZwidth(0), dim(0)-XYZwidth(0)-1 do begin ; Go through the Y dimension
              inter=data(i:i+XYZwidth(0)-1,j:j+XYZwidth(1)-1, k:k+XYZwidth(2)-1) ; Zone of the image that is going to be analysed
              inter=reform(inter,TotalLength ) ; Convert inter to a 1D vector
              ind=sort(inter) ; Sort inter in ascending order
              data(i+Center(0),j+Center(1),k+Center(2))=inter(ind(TotalLength/2)) ; The central voxel is assigned the value in the middle of the sorted vector
          endfor
      endfor
  endfor
  return, data ; The function returns data based on neighboring values in order to remove noise
end

Function AppliedVolumeSENSE, data, SenseMatrix, SENSEFactor=SENSEFactor,outputRes=outputRes,ReducedFOV=ReducedFOV, QBodyCoil=QBodyCoil, XYShift=XYShift;(0), -XYShift(1))
; Function receiving as inputs the original data in the k-space, the Sensitivity Matrix, the matrix size, the number of slices, the number of coils, the number of cardiac phases,
; the reduced FOV dimensions, the SENSE factor, the information about the Q-body coil and the XY shift
  print, "Applying Volume SENSE reconstruction"
  volumen=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Final reconstructed image by SENSE, with dimensions Matrix size (Y x Z) and the number of slices (X)
  inter=data ; Auxiliary variable that keeps the original data
  inter=fft(fft(inter,-1,dimension=1, /double),-1, dimension=2, /double) ; The original data are FFT in 1D, then they are FFT in 2D
  steps=[outputRes(0)/SENSEFactor(0), outputRes(1)/SENSEFactor(1)] ; Number of times we have to calculate sets of unfolded pixels in the horizontal and vertical direction
;  if keyword_set(ReducedFOV) eq 0 then window, 3, xsize=8*outputRes(0), ysize=4*outputRes(1)
;sacar_image=0
;stop
  for i=0, outputRes(2)-1 do begin ; Go through all the coils
    volumenInter=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Auxiliary variable of the size of our final reconstructed image
    for j=0, (SENSEFactor(0)-1) do volumenInter(j*steps(0):(j+1)*steps(0)-1, 0:steps(1)-1,*,*)=reform(inter(*,*,i,*)) ; Assign the FFT of the original data in each slice to each of the aliased pixels in the horizontal direction
    for j=0, (SENSEFactor(1)-1) do volumenInter(*, j*steps(1):(j+1)*steps(1)-1,*,*)=volumenInter(*, 0:steps(1)-1,*) ; Assign the aliased pixels in the vertical direction to the variable steps
    if keyword_set(ReducedFOV) eq 0 then volumenInter=shift(volumenInter, 0, outputRes(1)/4, 0)
    ; If the reduced FOV has not been defined, shift volumenInter a quarter of the number of columns down
;    if sacar_image ne 0 then tvscl,(abs(volumenInter(*,*,50))), i
;    if keyword_set(ReducedFOV) eq 0 then tvscl,(abs(volumenInter(*,*,60))), i
    volumen +=volumenInter*SenseMatrix(*,*,*,i) ; The final volume is added to the aliased images in the k-space multiplied by the sensitivity matrix
  endfor
  volumen=volumen/(SENSEFactor(0)*SENSEFactor(1)) ; The final volume is reduced by the SENSE factor in the horizontal and vertical directions
  if keyword_set(QBodyCoil) then volumen=volumen*abs(QBodyCoil) ; If there is available an image obtained with the Q-body coil with a priori information, its mask is multiplied by the final volume, in order to enhance SNR and CNR
  return, volumen ; The function returns the final volume
end


Function BuildSENSEMatrix, SenseMaps,outputRes=outputRes, SENSEFactor, SVDthreshold=SVDthreshold, QBodyCoil=QBodyCoil
   print, 'Build SENSEMatrix Reconstruction'
   ; The function builds the sensitivity matrix starting from the sensitivity maps, the matrix size, the number of slices, the number of coils and the cardiac phase number,
   ; the SENSE factor, the SVD threshold and the Q-body coil information
   if keyword_set(SVDthreshold) ne 1 then SVDthreshold=1e-5 ; If no SVD threshold is stated, it is set to 10^(-5)

   SENSEmatrix=complexarr(outputRes(0),outputRes(1),outputRes(3),outputRes(2))
   ; The sensitivity matrix is a complex array of dimensions: matrix size X number of slices X number of channels
   for k=0, outputRes(3)-1 do begin ;Loop over slices
       SenseMapsInter=reform(SenseMaps(*,*,*,k)) ; Auxiliary variable which is assigned the values of the sensitivity map of each slice

       if keyword_set(QBodyCoil) then for NCoils=0, 19 do SenseMapsInter(*,*,NCoils)=SenseMapsInter(*,*,NCoils)*reform(abs(QBodyCoil(*,*,k)))
       ; If there is an image with information a priori from the Q-body coil, the sensitivity maps are multiplied by the mask obtained with the Q-body coil, in each coil
       matrix=complexarr(SENSEFactor(0)*SENSEFactor(1),outputRes(2)) ; Complex array with dimensions the product of the SENSE factors in both directions as rows and the number of channels as columns
       for j=0, round(outputRes(1)/SENSEFactor(1))-1 do begin ; Go for every aliased voxel in the vertical direction
           for i=0, round(outputRes(0)/SENSEFactor(0))-1 do begin ; Go for every aliased voxel in the horizontal direction
               for l=0, SENSEFactor(1)-1 do begin
                   for p=0, SENSEFactor(0)-1 do begin
                       matrix(SENSEFactor(0)*l+p,*)=SenseMapsInter(p*round(outputRes(0)/SENSEFactor(0))+i,l*round(outputRes(1)/SENSEFactor(1))+j,*) ; Assign the matrix the different values of the sensitivity maps in the indexes of the aliased voxels
                   endfor
               endfor
               la_svd, matrix,w,u,v, /double ; Apply the SVD in order to reduce noise and enhance SNR and CNR
               ind=where(w gt 1e-5, count) ; Indexes where the SVs are higher than SVDthreshold--> TRUNCATED SVD
               if count gt 0 then begin ; If there SVs above the threshold, in u, v and w keep just the values in the indexes whose values are above the threshold
                   w=w(ind)
                   u=u(ind,*)
                   v=v(ind,*)
                   w=w/(w^2.+(SVDthreshold*w(0)^2.0)) ; Build a matrix for Tikhonov regularization, basing on the truncated values. Alpha= SVD threshold*Maximum SV
                   inv_matrix=v##diag_matrix(w)## transpose(conj(u)) ; Inverse matrix= V·Tikhonov matrix·U*
                   for l=0, SENSEFactor(1)-1 do begin ; Go through the SENSE factors in the horizontal and vertical directions
                       for p=0, SENSEFactor(0)-1 do begin
                            SENSEmatrix(p*round(outputRes(0)/SENSEFactor(0))+i,l*round(outputRes(1)/SENSEFactor(1))+j,k,*)=inv_matrix(*,SENSEFactor(0)*l+p)
                            ; Assign the voxels in the aliased indexes the columns of the inverse matrix, which will be the final values of the unaliased voxels in the final reconstructed image
                       endfor
                   endfor
               endif
            endfor
        endfor
   endfor
   SENSEmatrix=SENSEmatrix*SENSEFactor(0)*SENSEFactor(1); to compensate by 0 filling in non acquired samples

   return, SENSEmatrix
end

function SENSE_all_in_one, Data, SenseMaps,outputRes=outputRes, SENSEFactor=SENSEFactor,SVDthreshold=SVDthreshold,QBodyCoil=QBodyCoil,fullSamplingMask=fullSamplingMask
   ;     Ky          Kz            nCoil       Slices      Phases
   ;outputRes(0),outputRes(1),outputRes(2),outputRes(3),outputRes(4));
  Final3DCine=complexarr(outputRes(0),outputRes(1),outputRes(3),outputRes(4)) ; Final reconstructed image of dimensions: Matrix size X Number of slices X Cardiac phase number
  SenseMatrix=BuildSENSEMatrix(SenseMaps,outputRes=outputRes, SENSEFactor=SENSEFactor, SVDthreshold=SVDthreshold, QBodyCoil=QBodyCoil););  Obtain the Sensitivity Matrix

  if keyword_set(fullSamplingMask) then ind=where(fullSamplingMask gt 0)
  ; Look if the full sampling mask (the image with the full FOV) is defined. If this happens, look where the sampling mask is positive
  for i=0,outputRes(4)-1 do begin ; Go through every cardiac phase
      DynData=reform(Data(*,*,*,*,i)) ; Create a "matrix" per cardiac phase
      DynData=AppliedVolumeSENSE(DynData, SenseMatrix,outputRes=outputRes, SENSEFactor=SENSEFactor, QBodyCoil=QBodyCoil) ; Obtain the final reconstructed volume
      if keyword_set(fullSamplingMask) then begin ; Look if the full sampling mask (the image with the full FOV) is defined.
          a=dcomplexarr(outputRes(0), outputRes(1), outputRes(3)) ; Complex array with dimensions the Matrix Size X Number of slices
          a(ind)=DynData(ind) ; Where the full sampling mask is positive, the values of a are assigned to the values of the final reconstructed volume
          DynData=a*fullSamplingMask ; Where the full sampling mask is 0, the reconstructed volume is set to 0 and where the full sampling mask is positive, the reconstructed volume is conserved
          a=0
      endif
      Final3DCine(*,*,*,i)=DynData ; Assign the final reconstructed volume with the mask to the final result
  endfor
  return, Final3DCine ; The function returns the final reconstructed 3D volume
end

