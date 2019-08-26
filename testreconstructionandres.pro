pro testReconstructionAndres
; Obtain the images with the raw data, the sensitivity maps and the Q-body coil with a priori information
outputRes=[324.,196.,16.,20.,1.] ; Vector telling us the matrix size (324 x 176), the number of coils, the number of slices and the number of cardiac phases (just 1)
data=complexarr(outputRes(0), outputRes(1), outputRes(2), outputRes(3))
SenseMaps=complexarr(outputRes(0), outputRes(1), outputRes(2), outputRes(3))
QBodyCoil=complexarr(outputRes(0), outputRes(1), outputRes(3))
dir='H:\CORONARIAS\' ; Main directory
dir1=dir+'Data\' ; Directory where to read the raw data
dir2=dir+'SenseMaps\' ; Directory where to read the sensitivity maps
dir3=dir+'QBodyCoil\' ; Directory where to read the a priori information from the Q-Body coil
dir4=dir+'FinalVolumeAndres\'
for i=0, outputRes(3)-1 do begin ; Go through all the slices
	QBodyCoil(*,*,i)= read_png(dir3+'QBodyCoilSlice'+string(i)+'.png')
	for j= 0, outputRes(2)-1 do begin ; Go through all the coils
		data(*, *, j, i)= float(read_png(dir1+'SavedDataSlice'+string(i)+'Coil'+string(j)+'.png'))
		SenseMaps(*, *, j, i)= read_png(dir2+'SenseMapSavedSlice'+string(i)+'Coil'+string(j)+'.png')
	endfor
endfor
SENSEFactor=[3.0, 2.0] ; SENSE factors in the horizontal and vertical directions
; Data aliasing
data=data(SENSEFactor(0)*findgen(outputRes(0)/SENSEFactor(0)),*,*,*)
data=data(*,SENSEFactor(1)*findgen(outputRes(1)/SENSEFactor(1)),*,*)
for i= 0, outputRes(2)-1 do begin
	SenseMaps(*,*,i,*)*=QBodyCoil; Multiply all the sensitivity maps by the a priori information obtained with the Q-body coil
endfor
; Obtain the sensitivity matrix (32 x 6)
steps=[outputRes(0)/SENSEFactor(0), outputRes(1)/SENSEFactor(1)]
for i=0, outputRes(3)-1 do begin
	for j=0, outputRes(2)-1 do begin
	 	data(*,*,j,i)=fft(data(*,*,j,i),1) ; Aliased image in the spatial domain
	 	data(*,*,j,i)=shift(data(*,*,j,i), steps(0), steps(1))
	 endfor
endfor
; Sensitivity matrix and its hermitian
SenseMatrix=complexarr(SENSEFactor(0)*SENSEFactor(1), outputRes(2))
SenseMatrixH=complexarr(outputRes(2),SENSEFactor(0)*SENSEFactor(1))
; Threshold for the truncated SVD
SVDthreshold=1e-5
; Unfolding matrix
Unf=complexarr(SENSEFactor(0)*SENSEFactor(1), outputRes(2))
; Vector with the raw data
a= complexarr(outputRes(2))
; Vector with the final values of the image
fin=complexarr(SENSEFactor(0)*SENSEFactor(1))
finalVolume=complexarr(outputRes(0), outputRes(1), outputRes(3))
for k=0, outputRes(3)-1 do begin ; Go through all slices
	for i=0, steps(1)-1 do begin ; Go in the vertical direction
		for j=0, steps(0)-1 do begin ; Go in the horizontal direction
			for n=0, SENSEFactor(0)-1 do begin ; Horizontal SENSE factor
				for m=0, SENSEFactor(1)-1 do begin ; Vertical SENSE factor
					SENSEMatrix(SENSEFactor(0)*m+n,*)=SenseMaps(j+n*steps(0),i+m*steps(1),*,k)
				endfor
			endfor
		endfor
	endfor
	SENSEMatrixH=transpose(conj(SenseMatrix)) ; Hermitian
	la_svd, SENSEMatrixH##SENSEMatrix, W, U, V ; Singular value decomposition for the regularization
	UH=transpose(conj(U))
	for q=0, SENSEFactor(0)*SENSEFactor(1)-1 do begin
		if W(q) ge SVDthreshold then W(q)=W(q)/(W(q)^2.+(W(0)^2.0)*SVDthreshold) else W(q)=0 ; Tikhonov regularization
		diag=diag_matrix(W)
	endfor
	Unf=(V##diag##UH)##SENSEMatrixH
	for r=0, steps(1)-1 do begin ; Go in the vertical direction
		for o=0, steps(0)-1 do begin ; Go in the horizontal direction
			for p= 0, outputRes(2)-1 do begin
				a(p)=data(o,r,p,k)
			endfor
			fin=Unf##a
		    for n=0, SENSEFactor(0)-1 do begin ; Horizontal SENSE factor
				for m=0, SENSEFactor(1)-1 do begin ; Vertical SENSE factor
					finalVolume(o+n*steps(0), r+m*steps(1), k)=fin(SENSEFactor(0)*m+n)
				endfor
			endfor
		endfor
	endfor
	finalVolume(*,*,k)=abs(finalVolume(*,*,k))
	finalVolume(*,*,k)*=QBodyCoil(*,*,k)
	finalVolume(*,*,k)=finalVolume(*,*,k)*255/max(finalVolume(*,*,k))
	write_png, dir4+'FinalVolumeAndresSlice'+string(k)+'.png', finalVolume(*,*,k)
endfor
stop
end
