function final_reconstruction, SVDthreshold, outputRes, data, SENSEFactor, SenseMaps, QBodyCoil
SENSEmatrix=complexarr(outputRes(0),outputRes(1),outputRes(3),outputRes(2))
volumen=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Final reconstructed image by SENSE, with dimensions Matrix size (Y x Z) and the number of slices (X)
inter=data ; Auxiliary variable that keeps the original data
inter=fft(fft(temporary(inter),-1,dimension=1, /double),-1, dimension=2, /double) ; The original data are FFT in 1D, then they are FFT in 2D
steps=[round(outputRes(0)/SENSEFactor(0)), round(outputRes(1)/SENSEFactor(1))]
volumenInter=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Auxiliary variable of the size of our final reconstructed image
SenseMapsAux=complexarr(outputRes(0),outputRes(1),outputRes(3),outputRes(2))
; The sensitivity matrix is a complex array of dimensions: matrix size X number of slices X number of channels
for k=0, outputRes(3)-1 do begin ; Loop over slices
   SenseMapsInter=reform(SenseMaps(*,*,k,*)) ; Auxiliary variable which is assigned the values of the sensitivity map of each slice
   for NCoils=0, outputRes(2)-1 do SenseMapsInter(*,*,NCoils)=temporary(SenseMapsInter(*,*,NCoils))*reform(abs(QBodyCoil(*,*,k)))
   SenseMapsAux(*,*,k,*)=SenseMapsInter
   ; If there is an image with information a priori from the Q-body coil, the sensitivity maps are multiplied by the mask obtained with the Q-body coil, in each slice
   matrix=complexarr(SENSEFactor(0)*SENSEFactor(1),outputRes(2)) ; Complex array with dimensions the product of the SENSE factors in both directions as rows and the number of channels as columns
   for j=0, steps(1)-1 do begin ; Go for every aliased voxel in the vertical direction
       for i=0, steps(0)-1 do begin ; Go for every aliased voxel in the horizontal direction
            for l=0, SENSEFactor(1)-1 do begin
               for p=0, SENSEFactor(0)-1 do begin
                   matrix(SENSEFactor(0)*l+p,*)=SenseMapsInter(p*round(outputRes(0)/SENSEFactor(0))+i,l*round(outputRes(1)/SENSEFactor(1))+j,*) ; Assign the matrix the different values of the sensitivity maps in the indexes of the aliased voxels
               endfor
            endfor
            la_svd, matrix,w,u,v, /double ; Apply the SVD in order to reduce noise and enhance SNR and CNR
            ind=where(w gt 1e-5, count) ; Indexes where the SVs are higher than SVDthreshold--> TRUNCATED SVD
            if count gt 0 then begin ; If there SVs above the threshold, in u, v and w keep just the values in the indexes whose values are above the threshold
                w=temporary(w(ind))
                u=temporary(u(ind,*))
                v=temporary(v(ind,*))
                w=w/(w^2.+SVDthreshold*w(0)^2.0) ; Build a matrix for Tikhonov regularization, basing on the truncated values. Alpha= SVD threshold*Maximum SV
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
SENSEmatrix=temporary(SENSEmatrix)*SENSEFactor(0)*SENSEFactor(1); to compensate by 0 filling in non acquired samples
; OBTAIN THE FINAL RECONSTRUCTED VOLUME
for i=0, outputRes(2)-1 do begin ; Go through all the coils
	for j=0, (SENSEFactor(0)-1) do volumenInter((j*steps(0)):((j+1)*steps(0))-1, 0:steps(1)-1,*,*)=reform(inter(*,*,i,*)) ; Assign the FFT of the original data in each slice to each of the aliased pixels in the horizontal direction
    for j=0, (SENSEFactor(1)-1) do volumenInter(*, j*steps(1):(j+1)*steps(1)-1,*,*)=volumenInter(*, 0:steps(1)-1,*) ; Assign the aliased pixels in the vertical direction to the variable steps
    if keyword_set(ReducedFOV) eq 0 then volumenInter=shift(temporary(volumenInter), 0, outputRes(1)/4, 0)
    ; If the reduced FOV has not been defined, shift volumenInter a quarter of the number of columns down
    volumen=temporary(volumen)+volumenInter*SenseMatrix(*,*,*,i) ; The final volume is added to the aliased images in the k-space multiplied by the sensitivity matrix
endfor
volumen=temporary(volumen)/(SENSEFactor(0)*SENSEFactor(1)) ; The final volume is reduced by the SENSE factor in the horizontal and vertical directions
if keyword_set(QBodyCoil) then volumen=temporary(volumen)*abs(QBodyCoil) ; If there is available an image obtained with the Q-body coil with a priori information, its mask is multiplied by the final volume, in order to enhance SNR and CNR
stop
return, volumen
end