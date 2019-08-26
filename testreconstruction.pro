pro TestReconstruction
; The procedure tests the data in the cpx and tries to reconstruct the images from the data
dir='H:\CORONARIAS\' ; Main directory
dir1=dir+'Data\' ; Directory where to read the raw data
dir2=dir+'SenseMaps\' ; Directory where to read the sensitivity maps
dir3=dir+'QBodyCoil\' ; Directory where to read the a priori information from the Q-Body coil
dir4=dir+'FinalVolume\'
RefProto=ReadProtocol(dir+'Ref_XL_Torso.txt') ; Obtain the FOV, the voxel size, the TFE factor, the TFE shot duration and shot interval, the angulation and other aspects from the Reference protocol
ScanProto=ReadProtocol(dir+'CORONARIAS+EDEMA.txt') ; Same for the Scan protocol
;Fractional SENSE correction
RealImageFOV=[ScanProto.FOV(1),ScanProto.FOV(2),ScanProto.FOV(0)] ; Image FOV (YxZ and X is represented as slices)
; LINE DOWN: SPECIFIC FOR FRACTIONAL SENSE FACTORS
ScanProto.FOV(1:2)=ScanProto.FOV(1:2)*ceil(ScanProto.SENSE)/ScanProto.SENSE ; FOV modified by the SENSE factors in different directions. IN CASE OF FRACTIONAL SENSE FACTORS
SENSEfactor= ScanProto.SENSE
TotalSlices  = round(ScanProto.FOV(0)/ScanProto.VOXELSIZE(0)) ; Slices are obtained by dividing the FOV along X over the voxel size along X
SlicesCenter = 0 ; positive is foot direction negative head direction

;RR_interval=((60./float(ScanProto.CardFrequency))*1000.)/float(ScanProto.CardPhases) ; Time between two Rs in an ECG, in ms. Inverse of the cardiac frequency
; Search coil and data files in the same directory
CoilHeaderCPX=file_search(dir+'cpx*.list')
CoilDataCPX=file_search(dir+'cpx*.data')
DataHeader=file_search(dir+'raw*.list')
DataData=file_search(dir+'raw*.data')

Coils=read_cpx(CoilHeaderCPX, CoilDataCPX, QBodyCoil, NCoils)
outputRes=float([324,196, NCoils, 20, 1]) ; Vector with the matrix size of the final image, the number of coils, the number of slices and the number of cardiac phases
; FOV shift due to different FOV sampling scheme for fractional SENSE factors
; FOLLOWING TWO LINES ARE SPECIFIC FOR FRACTIONAL SENSE FACTORS
FOVShift=(ceil(ScanProto.SENSE)-ScanProto.SENSE)/ceil(ScanProto.SENSE) ; Shift of the FOV: depends on the SENSE factor
ScanProto.SENSE=ceil(ScanProto.SENSE) ; The SENSE factor is approximated to the closer and higher integer
FOVCorrection=scanproto.FOV ; Final FOV
data=complexarr(outputRes(0), outputRes(1), outputRes(2), outputRes(3))
SenseMaps=complexarr(outputRes(0), outputRes(1), outputRes(2), outputRes(3))
; Read the already obtained data, sensitivity maps and QBodyCoil information
QBodyCoil= complexarr(outputRes(0), outputRes(1), outputRes(3))
for j=0, outputRes(3)-1 do begin ; Go through all the slices
QBodyCoil(*,*,j)= read_png(dir3+'QBodyCoilSlice'+string(j)+'.png')
	for i=0, NCoils-1 do begin ; Go through all the coils
		data(*, *, i, j)= read_png(dir1+'SavedDataSlice'+string(j)+'Coil'+string(i)+'.png')
		SenseMaps(*, *, i, j)= read_png(dir2+'SenseMapSavedSlice'+string(j)+'Coil'+string(i)+'.png')
	endfor
endfor
aux=data
for i=0, NCoils-1 do begin
	for j=0, outputRes(3)-1 do begin
		data=shift(abs(aux(*,*,i,j)),outputRes(0)/2,outputRes(1)/2)
	endfor
endfor
data=aux
; FILTERING
;maskKspace=complexarr(outputRes(0), outputRes(1))
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
;maskKspace=shift(maskKspace,outputRes(0)/2.,outputRes(1)/2.) ; Shift the mask with half of the dimensions of the matrix size
;maskKspace(radio(0):HalfScan(0),*)=0 ; The mask between the first and last indexes of rows and columns where it was 0 is now 0
;maskKspace(*,radio(1):HalfScan(1))=0
;
;filter=shift(filter,outputRes(0)/2.,outputRes(1)/2.) ; Shift the filter with half of the dimensions of the matrix size
;filter(LowerRadio(0):LowerHalfScan(0),*)=0 ; The filter between the first indexes of rows and columns where it was 0, minus 3, and the last indexes of rows and columns where it was 0, plus 3, is now 0
;filter(*,LowerRadio(1):LowerHalfScan(1))=0
;
;HomodineFilter=maskKspace ; The final mask is a LP Homodyne filter
;
;filter=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; The filter is unshifted
;filter=(convol(filter, FilterGen(25,19),/EDGE_t)) ; The filter is now convolved with a Gaussian filter of 25 samples and a width of 19, convolving also in the edges
;
;maskKspace=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; The mask is unshifted
;
;weight=fltarr(1,outputRes(1)) ; Vector of weights with the dimensions of Z
;weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./40000) ; Weights are normally distributed
;weight=rebin(weight, outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image
;
;a=fltarr(outputRes(0),1) ; Vector of weights with the dimensions of Y
;a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./40000) ; Weights are normally distributed
;a=rebin(a, outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image
;
;weight *= a ; Weights in Y and Z directions are now multiplied between them
;
;filter=filter*shift(HomodineFilter, -outputRes(0)/2.,-outputRes(1)/2.)*weight
;; Homodyne filter (the initial filter) is shifted with half of the dimensions of the final matrix size and multiplied by the convolved version of the filter and by the Gaussian weights
;weight=0
;
;filter=reform(shift(filter, -outputRes(0)/2.,-outputRes(1)/2.),  outputRes(0),outputRes(1),1,1) ; The filter is shifted again and adapted to the matrix size of the final image
;filter=rebin(filter,outputRes(0),outputRes(1),1,outputRes(3)) ; The filter is adapted to the dimensions of the final matrix size and the slice number
;
;for i=0, outputRes(2)-1 do begin
;    data(*,*,i,*)=data(*,*,i,*)*filter ; Filtered data
;endfor
;filter=reform(filter) ; Dimensions with just one value are removed
;
;
;
;mask=fltarr(outputRes(0),outputRes(1),outputRes(2),1,outputRes(4)) ; Float array analyzed in just one slice
;ind=where(abs(data(*,*,*,0)), count) ; Indexes of the absolute value of our data
;mask(ind)=1. ; The mask is one in all the indexes of our data
;
;mask=total(mask,5) ; Sum of all the elements of data for the same cardiac phase
;ind=where(abs(mask), count) ; Keep the mask indexes
; Obtain the aliased data
; DATA SAMPLING
Data=Data(ScanProto.SENSE(0)*findgen(outputRes(0)/(ScanProto.SENSE(0))),*,*,*) ; Data is restricted in Y and Z directions by the SENSE factors in both directions
Data=Data(*,ScanProto.SENSE(1)*findgen(outputRes(1)/(ScanProto.SENSE(1))),*,*)
; COIL SELECTION
ChannelSNR=total(total(abs(data(0,0,*,*)),2),2)/outputRes(3) ; Sum of the absolute value of the first element (Row 0, Column 0) in all the slices, divided by the number of slices and multiplied by the number of cardiac phases
ind=reverse(sort(ChannelSNR)) ; Absolute values of the first elements of data in each slice sorted in descending order
ChannelSelection=where(ChannelSNR gt (mean(ChannelSNR(ind))-2.0*stddev(ChannelSNR(ind))), NumSelectedChannels)
; Indexes where the absolute values of the first elements of the data in each slice are higher than the sorted absolute values of the first elements (from the first slice to the slice 2/3 of the of the total number of coils)
; minus twice the standard deviation of the sorted absolute values (from the first element to the slice whose number is equal to half of the coils used
; The number of times this happens is saved in "NumSelectedChannels"
if NumSelectedChannels lt 8 then NumSelectedChannels=outputRes(2)
; See if NumSelectedChannels is lower than 8. In that case, the number of selected channels equals the number of coils
if NumSelectedChannels lt outputRes(2) then begin
	;print, ChannelSelection ; Print the channel indexes if the number of selected channels is lower than the number of coils
	ChannelSNR=0
	data=data(*,*,ChannelSelection,*) ; The data is now just kept in the channels we selected
	SenseMaps=temporary(SenseMaps(*,*,ChannelSelection,*)) ; Sensitivity Maps are saved temporarily in the selected channels
    outputRes(2)=NumSelectedChannels ; The number of channels in outputRes is equal to the number of selected channels
endif
; Build the sensitivity matrix (Ncoils X NCoils) and apply SENSE
SVDthreshold=1e-5
SENSEmatrix=complexarr(outputRes(0),outputRes(1),outputRes(3),outputRes(2))
volumen=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Final reconstructed image by SENSE, with dimensions Matrix size (Y x Z) and the number of slices (X)
inter=data ; Auxiliary variable that keeps the original data
inter=fft(fft(inter,-1,dimension=1, /double),-1, dimension=2, /double) ; The original data are FFT in 1D, then they are FFT in 2D
steps=[round(outputRes(0)/SENSEFactor(0)), round(outputRes(1)/SENSEFactor(1))]
volumenInter=complexarr(outputRes(0),outputRes(1),outputRes(3)) ; Auxiliary variable of the size of our final reconstructed image
; The sensitivity matrix is a complex array of dimensions: matrix size X number of slices X number of channels
for k=0, outputRes(3)-1 do begin ; Loop over slices
   SenseMapsInter=reform(SenseMaps(*,*,*,k)) ; Auxiliary variable which is assigned the values of the sensitivity map of each slice
   for NCoils=0, outputRes(2)-1 do SenseMapsInter(*,*,NCoils)=SenseMapsInter(*,*,NCoils)*reform(abs(QBodyCoil(*,*,k)))
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
SENSEmatrix*=SENSEFactor(0)*SENSEFactor(1); to compensate by 0 filling in non acquired samples
; OBTAIN THE FINAL RECONSTRUCTED VOLUME
for i=0, outputRes(2)-1 do begin ; Go through all the coils
	for j=0, (SENSEFactor(0)-1) do volumenInter((j*steps(0)):((j+1)*steps(0))-1, 0:steps(1)-1,*,*)=reform(inter(*,*,i,*)) ; Assign the FFT of the original data in each slice to each of the aliased pixels in the horizontal direction
    for j=0, (SENSEFactor(1)-1) do volumenInter(*, j*steps(1):(j+1)*steps(1)-1,*,*)=volumenInter(*, 0:steps(1)-1,*) ; Assign the aliased pixels in the vertical direction to the variable steps
    if keyword_set(ReducedFOV) eq 0 then volumenInter=shift(volumenInter, 0, outputRes(1)/4, 0)
    ; If the reduced FOV has not been defined, shift volumenInter a quarter of the number of columns down
    volumen +=volumenInter*SenseMatrix(*,*,*,i) ; The final volume is added to the aliased images in the k-space multiplied by the sensitivity matrix
endfor
volumen=volumen/(SENSEFactor(0)*SENSEFactor(1)) ; The final volume is reduced by the SENSE factor in the horizontal and vertical directions
if keyword_set(QBodyCoil) then volumen=volumen*abs(QBodyCoil) ; If there is available an image obtained with the Q-body coil with a priori information, its mask is multiplied by the final volume, in order to enhance SNR and CNR
for i=0, outputRes(3)-1 do begin
	;volumen(*,*,i)=abs(volumen(*,*,i))
	;volumen(*,*,i)=volumen(*,*,i)*255/max(volumen(*,*,i))
	write_png, dir4+'FinalVolumeSlice'+string(i)+'.png', volumen(*,*,i)
endfor
stop
end