;pro HomodyneRecon, data,outputRes=outputRes, filter=filter
;; HOMODYNE FILTER
;; Inputs:
;; data obtained from the scan in the k-space
;; outputRes: matrix size of the reconstructed image
;; filter: filter with the matrix size of outputRes
;print, "Homodyne filtering"
;  filter=float(filter) ; Filter is a floating point
;  filter=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; Filter is shifted along the matrix size
;  mask=filter ; Filter is saved in variable "mask"
;  ;float(erode(filter, replicate(1.,10,10)))
;  a=shift(reverse(shift(reverse(filter,1),1,0),2),0,1)
;  ; Filter is reversed in its columns and shifted 1 column to the right. Then, it is reversed along its rows and shifted one row down
;  ; Changes are saved in variable "a", float
;  b=filter+a ; Changes applied in "a" are added to the original filter, the sum is saved in variable "b"
;  ind=where(b eq 2, count) ; Find the indexes where "b" equals 2 and count the times this happens
;    a=fltarr(outputRes(0),outputRes(1))
;    a(ind)=1 ; In the indexes where b=2, a=1
;    a=a+shift(a,10,0)
;    a=a+shift(a,0,10)
;    ; Shift "a" 10 columns to the right and 10 rows down. Add all this to the previous "a"
;    ind=where(a gt 0) ; Indexes where "a" is positive
;  filter(ind)=0.7 ; If "a" is positive, filter=0.7 at those indexes
;  filter=filter*mask ; Multiply the filter by the initial filter saved in "mask"
;  filter=convol(filter, FilterGen(30,19),/EDGE_t) ; Convolution of filter with a gaussian kernel of 30x30 with a width of 19.
;  ; The values at the edges are repeated outside in order to be able to convolve in the edges
;     weight=fltarr(1,outputRes(1)) ; Vector with the length of the number of columns of the matrix size. Contains float
;     weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./50000) ; Contains exponential weights which are 1 in the center and minimum in the beginning and in the end
;     weight=rebin(weight, outputRes(0), outputRes(1))
;     ; With the weights computed, the vector "weight" is resized into a matrix with the dimensions of the original matrix size
;     a=fltarr(outputRes(0),1)
;     a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./50000)
;     a=rebin(a, outputRes(0), outputRes(1))
;     ; "a" is the same than "weight", but with the number of rows
;     weight *= a ; "weight" and "a" are multiplied among them to generate a matrix of weights
;
;  filter=filter*weight ; The convoluted filter is multiplied by the weight matrix, having now a Gaussian weight
;  tvscl, filter ; Display the weighted filter in screen
;  filter=shift(filter,outputRes(0)/2.,outputRes(1)/2.) ; Filter is unshifted
;  filter=rebin(filter, outputRes(0),outputRes(1),outputRes(3)) ; Filter is resized to the original number of rows and columns and to the slice number
;
;  data=fft(fft(data, 1,dimension=1), 1, dimension=2); pasa al espacio-K
;  ; Compute the IFFT along only 1 dimension of the data. Then, compute the IFFT of the first IFFT along two dimensions
;  for i=0, outputRes(4)-1 do begin
;    inter=data(*,*,*,i)*filter ; Filtered data saved in "inter"
;    data(*,*,*,i)=inter
;  endfor
;  data=fft(fft(data, -1,dimension=1), -1, dimension=2); pasa al espacio real
;  ; FFT of the filtered data in 1D. Then, a FFT is applied again in 2D
;end

Function Generate_KSpace_mask, maskKspace
; K SPACE MASK
    maskKspace=total(total(total(maskKspace,3),3),3) ; Sum of all the elements along the 3rd dimension of mask space
    dim=size(maskKspace, /dimension) ; Size of the mask
    KSpace_mask=shift(maskKspace, -dim(0)/2.,-dim(1)/2.) ; Shifted mask
    KSpace_mask=smooth(KSpace_mask,5) ; Smooth the original mask with a width of 5
    ind=where(KSpace_mask gt 0) ; Indexes where the mask is positive
    KSpace_mask(ind)=1. ; In the indexes where the mask is positive, now the mask is 1
    KSpace_mask=erode(KSpace_mask,replicate(1.,5,5)) ; Using a kernel of 5x5, the mask is eroded (binary operation that reduces the size of the object)
    KSpace_mask=shift(KSpace_mask, dim(0)/2.,dim(1)/2.) ; Unshift
    return, KSpace_mask ; Final mask in the k-space, which has been shifted, smoothened, made binary, eroded and unshifted again
end

pro D3SENSE
  argument=command_line_args(count=count)

  dir=argument(0)
;
dir='H:\CORONARIAS\' ; Directory where to read the raw data. IT MAY HAVE CHANGED!!!
;  stop
; READ DIRECTORY, PROTOCOLS AND CPX: OBTAIN DATA
loadct,0 ; Load color table number 0
  if file_test(dir) then begin ; Look for the file specified in the directory
     genVideo=1
     FOVReduction=0
     Recursive=0
     mask=1
     fullSamplingMask=1
     WriteData=1
     RingingFilterOrder=50
     XYshift=[10,-5];X=Negative Left must be multiplos de 6; Y=Negative Anterior must be multiplos de 4

     if file_test(dir+'ReconConfig.txt') then begin
          params=((read_ascii(dir+'ReconConfig.txt')).field1)(*,1)
          ; Save the initial values as an array in a struct called "params"
          genVideo=params(0)
          FOVReduction=params(1)
          Recursive=params(2)
          mask=params(3)
          fullSamplingMask=params(4)
          XYshift=params(5:6)
     endif

     ReadingTime=0.0
     ReconTime=0.0
     SensesitivitiesTime=0.0

     ReadingTime=SYSTIME(/SECONDS ) ; Time spent until here reading the script
     RefProto=ReadProtocol(dir+'Ref_XL_Torso.txt') ; Obtain the FOV, the voxel size, the TFE factor, the TFE shot duration and shot interval, the angulation and other aspects from the Reference protocol
     ;dialog_pickfile(filter='*.txt'))
     ScanProto=ReadProtocol(dir+'CORONARIAS+EDEMA.txt') ; Same for the Scan protocol

     ;Fractional SENSE correction
     RealImageFOV=[ScanProto.FOV(1),ScanProto.FOV(2),ScanProto.FOV(0)] ; Image FOV (YxZ and X is represented as slices)
     ; LINE DOWN: SPECIFIC FOR FRACTIONAL SENSE FACTORS
     ScanProto.FOV(1:2)=ScanProto.FOV(1:2)*ceil(ScanProto.SENSE)/ScanProto.SENSE ; FOV modified by the SENSE factors in different directions

     TotalSlices  = round(ScanProto.FOV(0)/ScanProto.VOXELSIZE(0)) ; Slices are obtained by dividing the FOV along X over the voxel size along X
     SlicesCenter = 0 ; positive is foot direction negative head direction

     ;RR_interval=((60./float(ScanProto.CardFrequency))*1000.)/float(ScanProto.CardPhases) ; Time between two Rs in an ECG, in ms. Inverse of the cardiac frequency
     ; Search coil and data files in the same directory
     CoilHeaderCPX=file_search(dir+'cpx*.list')
     CoilDataCPX=file_search(dir+'cpx*.data')

     DataHeader=file_search(dir+'raw*.list')
     DataData=file_search(dir+'raw*.data')

     Coils=read_cpx(CoilHeaderCPX, CoilDataCPX, QBodyCoil, NCoils)
     outputRes=float([640,640, NCoils, 20, 1]) ; Vector with the matrix size of the final image, the number of coils, the number of slices and the number of cardiac phases
     ; FOV shift due to different FOV sampling scheme for fractional SENSE factors
     ; FOLLOWING TWO LINES ARE SPECIFIC FOR FRACTIONAL SENSE FACTORS
     FOVShift=(ceil(ScanProto.SENSE)-ScanProto.SENSE)/ceil(ScanProto.SENSE) ; Shift of the FOV: depends on the SENSE factor
     ScanProto.SENSE=ceil(ScanProto.SENSE) ; The SENSE factor is approximated to the closer and higher integer
     ReadingTime=SYSTIME(/SECONDS )-ReadingTime ; Time spent between the first call to ReadingTime and now
     ; DATA OBTENTION
     SensesitivitiesTime=SYSTIME(/SECONDS ) ; Time spent until now
     FOVCorrection=scanproto.FOV ; Final FOV
     data=read_raw_list_by_TFEFactor(DataHeader, DataData, outputRes=outputRes, slice=SlicesCenter,ScanProto=ScanProto, noise,SENSEFactor=ScanProto.SENSE, mask=maskKspace,TimeMask=TimeMask);RR_interval=RR_interval,
     ;prepare data for reconstruction
     FOVCorrection=scanproto.FOV/FOVCorrection ; Correction of FOV

     ;window,1, xsize=outputRes(0)*4., ysize=outputRes(1)*4., title=dir ; Open a window with 4 times the matrix sizeof outputRes and with the name of the directory in the title
     ;for i=0, outputRes(4)-1 do tvscl,shift(alog(abs(data(*,*,0,outputRes(3)/2,i))),outputRes(0)/2,outputRes(1)/2) ,i
     ; Show the shifted logarithm of the middle slice in the k-space, with matrix size half of the one introduced in outputRes, in all cardiac phases
     if total(ScanProto.SENSE) lt 2.0 then ScanProto.SENSE= [3.0,2.0] ; If no SENSE is applied in

;NoiseDeco,SenseMaps,data, noise, outputRes=outputRes

;    dummy=NormData(data, maskKspace, outputRes=outputRes,TimeMask=TimeMask)
     ; K-SPACE MASK AND FILTERING
;     maskKspace=complexarr(outputRes(0),outputRes(1))
;     maskKspace=Generate_KSpace_mask(maskKspace) ; Obtain a mask to select where in the k-space we are going to work
;     vectorX=(where(maskKspace(*,0) eq 0)) ; Indexes where the rows of the mask are 0
;     vectorY=(where(maskKspace(0,*) eq 0)) ; Indexes where the columns of the mask are 0
;     radio=float([vectorX(0),vectorY(0)])+1 ; Contains the first indexes of rows and columns with zeros. A one is added
;     LowerRadio=radio-[3.,3.] ; The lower radio is 3 elements shorter than the initial radio
;     HalfScan=[vectorX(n_elements(vectorX)-1),vectorY(n_elements(vectorY)-1)] ; Halfscan contains the last elements of the rows and columns of the mask which are 0
;     LowerHalfScan=HalfScan+[3.,3.] ; Lower half scan is 3 elements higher than the original half scan
;
;     maskKspace=fltarr(outputRes(0),outputRes(1)) ; Float array with dimensions the matrix size of the final image
;     filter=fltarr(outputRes(0),outputRes(1)) ; Float array with dimensions the matrix size of the final image
;     for j=0, outputRes(1)-1 do begin
;       for i=0, outputRes(0)-1 do begin
;            radius=((i-outputRes(0)/2.)/radio(0))^2.+((j-outputRes(1)/2.)/radio(1))^2. ; Distance to the central pixel in the k-space normalized by the first elements where the mask is 0
;            if radius le 1 then maskKspace(i,j)=1.0 ; Looks if the radius is lower than 1, assigning then the mask a 1 --> LP filter in the k-space
;            LowerRadius=((i-outputRes(0)/2.)/LowerRadio(0))^2.+((j-outputRes(1)/2.)/LowerRadio(1))^2. ; Distance to the central pixel in the k-space normalized by the first elements where the mask is 0, minus 3
;            if LowerRadius le 1 then filter(i,j)=1.0 ; Looks if the lower radius is lower than 1, assigning then the filter a 1 --> LP filter in the k-space
;       endfor
;     endfor
;
;     maskKspace=shift(maskKspace,outputRes(0)/2.,outputRes(1)/2.) ; Shift the mask with half of the dimensions of the matrix size
;     maskKspace(radio(0):HalfScan(0),*)=0 ; The mask between the first and last indexes of rows and columns where it was 0 is now 0
;     maskKspace(*,radio(1):HalfScan(1))=0
;
;     filter=shift(filter,outputRes(0)/2.,outputRes(1)/2.) ; Shift the filter with half of the dimensions of the matrix size
;     filter(LowerRadio(0):LowerHalfScan(0),*)=0 ; The filter between the first indexes of rows and columns where it was 0, minus 3, and the last indexes of rows and columns where it was 0, plus 3, is now 0
;     filter(*,LowerRadio(1):LowerHalfScan(1))=0
;
;     HomodineFilter=maskKspace ; The final mask is a LP Homodyne filter
;
;     filter=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; The filter is unshifted
;     filter=(convol(filter, FilterGen(25,19),/EDGE_t)) ; The filter is now convolved with a Gaussian filter of 25 samples and a width of 19, convolving also in the edges
;
;     maskKspace=shift(filter,-outputRes(0)/2.,-outputRes(1)/2.) ; The mask is unshifted
;
;     weight=fltarr(1,outputRes(1)) ; Vector of weights with the dimensions of Z
;     weight(0,*)=exp(-(findgen(outputRes(1))-outputRes(1)/2.)^2./40000) ; Weights are normally distributed
;     weight=rebin(weight, outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image
;
;     a=fltarr(outputRes(0),1) ; Vector of weights with the dimensions of Y
;     a(*,0)=exp(-(findgen(outputRes(0))-outputRes(0)/2.)^2./40000) ; Weights are normally distributed
;     a=rebin(a, outputRes(0), outputRes(1)) ; Weights are now distributed in a matrix with size the matrix size of our image
;
;     weight *= a ; Weights in Y and Z directions are now multiplied between them
;
;     filter=filter*shift(HomodineFilter, -outputRes(0)/2.,-outputRes(1)/2.)*weight
;     ; Homodyne filter (the initial filter) is shifted with half of the dimensions of the final matrix size and multiplied by the convolved version of the filter and by the Gaussian weights
;     weight=0
;
;     filter=reform(shift(filter, -outputRes(0)/2.,-outputRes(1)/2.),  outputRes(0),outputRes(1),1,1) ; The filter is shifted again and adapted to the matrix size of the final image
;     filter=rebin(filter,outputRes(0),outputRes(1),1,outputRes(3)) ; The filter is adapted to the dimensions of the final matrix size and the slice number
;
;     for i=0, outputRes(2)-1 do begin
;        for j=0, outputRes(4)-1 do begin
;          data(*,*,i,*,j)=data(*,*,i,*,j)*filter ; The original data is now filtered for each coil and cardiac phase
;        endfor
;     endfor
;     filter=reform(filter) ; Dimensions with just one value are removed
;
;
;
;     mask=fltarr(outputRes(0),outputRes(1),outputRes(2),1,outputRes(4)) ; Float array analyzed in just one slice
;     ind=where(abs(data(*,*,*,0,*)), count) ; Indexes of the absolute value of our data
;     mask(ind)=1. ; The mask is one in all the indexes of our data
;
;     mask=total(mask,5) ; Sum of all the elements of data for the same cardiac phase
;     ind=where(abs(mask), count) ; Keep the mask indexes
     ; STATIC DATA
;     StaticKSpaceData=total(data,5) ; Sum of all the elements of data for the same cardiac phase --> STATIC DATA: only one cardiac phase
;     for i=0, outputRes(3)-1 do begin ; Go slice by slice
;        inter=Reform(StaticKSpaceData(*,*,*,i)) ; Change the dimensions of the static data again to 4 (Matrix size X Number of coils X Number of Slices)
;        inter(ind)=inter(ind)/mask(ind) ; Normalize the values of the static data dividing by the mask
;        StaticKSpaceData(*,*,*,i)=reform(inter,outputRes(0),outputRes(1),outputRes(2),1) ; Set the static normalized mask to the following dimensions: Matrix size X Number of coils X Number of slices
;     end

    ; StaticOriginalKSpaceDataIndex=where(abs(reform(StaticKSpaceData(*,*,0,*))) gt 0, countOriginal) ; Indexes where the static normalized mask is positive
     ; DATA SAMPLING
     Data=Data(ScanProto.SENSE(0)*findgen(outputRes(0)/(ScanProto.SENSE(0))),*,*,*,*) ; Data is restricted in Y and Z directions by the SENSE factors in both directions
     Data=Data(*,ScanProto.SENSE(1)*findgen(outputRes(1)/(ScanProto.SENSE(1))),*,*,*)
     ;for i=0, outputRes(4)-1 do tvscl, shift(alog(abs(data(*,*,0,outputRes(3)/2,i))),-outputRes(0)/ScanProto.SENSE(0)/4. ,-outputRes(1)/ScanProto.SENSE(1)/4. ),i
     ; Display the logarithm of the shifted image in the k-space depending on the matrix size and the SENSE factors for every cardiac phase
     ;imagen=tvrd() ; The images shown are saved in this variable
     ;WRITE_PNG, dir+'samplingAlongPhases.png', Imagen ; Images previously saved in a variable are converted to PNG format and saved in our directory


;     TimeMask=TimeMask(2*ScanProto.SENSE(0)*findgen(outputRes(0)/(2.*ScanProto.SENSE(0))),*,*,*,*) ; Sampling of time masks in rows and columns
;     TimeMask=TimeMask(*,2.*ScanProto.SENSE(1)*findgen(outputRes(1)/(2.*ScanProto.SENSE(1))),*,*,*)
;
;     StaticKSpaceData=StaticKSpaceData(ScanProto.SENSE(0)*findgen(outputRes(0)/ScanProto.SENSE(0)),*,*,*) ; Sampling of static data in the k-space in rows and columns
;     StaticKSpaceData=StaticKSpaceData(*,ScanProto.SENSE(1)*findgen(outputRes(1)/ScanProto.SENSE(1)),*,*)

     ;StaticReducedKSpaceDataIndex=where(abs(reform(StaticKSpaceData(*,*,0,*))) gt 0, countReduced)
     ; Indexes where the static data in the first slice is positive
;     print,scanproto.FOV(1) ; Print the FOV along the Z direction
;     ; COIL SELECTION
;     ChannelSNR=total(total(abs(data(0,0,*,*)),2),2)/outputRes(3) ; Sum of the absolute value of the first element (Row 0, Column 0) in all the slices, divided by the number of slices and multiplied by the number of cardiac phases
;     ind=reverse(sort(ChannelSNR)) ; Absolute values of the first elements of data in each slice sorted in descending order
;     ChannelSelection=where(ChannelSNR gt mean(ChannelSNR(ind(0:2.*round(outputRes(2)/3.))))-2.0*stddev(ChannelSNR(ind(0:outputRes(2)/2.))), NumSelectedChannels)
;     ; Indexes where the absolute values of the first elements of the data in each slice are higher than the sorted absolute values of the first elements (from the first slice to the slice 2/3 of the of the total number of coils)
;     ; minus twice the standard deviation of the sorted absolute values (from the first element to the slice whose number is equal to half of the coils used
;     ; The number of times this happens is saved in "NumSelectedChannels"
;     if NumSelectedChannels lt 8 then NumSelectedChannels=outputRes(2)
;     ; See if NumSelectedChannels is lower than 8. In that case, the number of selected channels equals the number of coils
;     if NumSelectedChannels lt outputRes(2) then begin
;          print, ChannelSelection ; Print the channel indexes if the number of selected channels is lower than the number of coils
;          ChannelSNR=0
;          maskKspace=maskKspace(*,*,ChannelSelection,*,*) ; The mask is now just kept in the channels we selected
;          data=data(*,*,ChannelSelection,*,*) ; The data is now just kept in the channels we selected
;          StaticKSpaceData=StaticKSpaceData(*,*,ChannelSelection,*) ; The static data is now just kept in the channels we selected
;     endif
     ; OBTENTION OF SENSITIVITY MAPS
     fullSamplingMask=1
     SensesitivitiesTime=SYSTIME(/SECONDS ) ; Time that has passed until this line has been read
     SenseMaps=SenseMapsGen(Coils, QBodyCoil, RefProto, ScanProto, outputRes=outputRes, QBThreshold=5e-2, mask=fullSamplingMask,FOVCorrection=FOVCorrection) ; Obtain the Sensitivity Maps
     SensesitivitiesTime=SYSTIME(/SECONDS )-SensesitivitiesTime ; Time spent obtaining the sensitivity maps
;
;     if (NumSelectedChannels lt outputRes(2)) then begin ; See if the number of selected channels is lower than the number of channels in OutputRes
;        SenseMaps=temporary(SenseMaps(*,*,*,ChannelSelection)) ; Sensitivity Maps are saved temporarily in the selected channels
;        outputRes(2)=NumSelectedChannels ; The number of channels in outputRes is equal to the number of selected channels
;     end

     Coils=0
     QBodyCoil=abs(QBodyCoil) ; Take the absolute value of the information from the Q-body coils
     ; RECONSTRUCTION
     ReconTime=SYSTIME(/SECONDS) ; Current time
     SENSEFactor=fltarr(2)
     SENSEFactor=ScanProto.SENSE
     Data=SENSE_all_in_one(Data, SenseMaps,outputRes=outputRes,SENSEFactor,SVDthreshold=1e-3,QBodyCoil=QBodyCoil,fullSamplingMask=fullSamplingMask)
     ; Obtain the final reconstructed volume in the full FOV
     ReconTime=SYSTIME(/SECONDS)-ReconTime ; Time spent reconstructing
     for i=0, outputRes(2)-1 do begin
     	for j=0, outputRes(3)-1 do begin
     	 write_png, dir+'FinalVolumeJavierCode\FinalVolumeJavierSlice'+string(j)+'Coil'+string(i)+'.png', Data(*,*,i,j)
     	endfor
     endfor
stop
;     HomodyneRecon, data,outputRes=outputRes, filter=HomodineFilter
     ; RESIZE
     filename=dir+'DICOM\IM_0001'
     info=read_dicom_info(filename)
     help, info,/struct
     if keyword_set(WriteData) then begin ; See if WriteData has been defined and resize the image to the full FOV
            print, 'Resizing to the final image resolution'
            FinalImageFOV=[ScanProto.FOV(1), Scanproto.FOV(2), Scanproto.FOV(0)] ; Final FOV defined in the Scan protocol
            NewDim=round(FinalImageFOV/[info.pixelsize(0),info.pixelsize(2),info.pixelsize(1)]) ; Final matrix size and number of slices
            NewDim=[NewDim,fix(outputRes(4))] ; To the matrix size and the number of slices, the number of cardiac phases is also included
            Inter=fltarr(NewDim(0),NewDim(1),NewDim(2),outputRes(4)) ; Auxiliary variable with the final matrix size, number of slices and cardiac phases
            SlicePos=[round((NewDim(0)-outputRes(0))/2.),round((NewDim(1)-outputRes(1))/2.),round((NewDim(2)-outputRes(3))/2.)] ; Position of the slice
            for i=0, outputRes(4)-1 do begin ; Go through all cardiac phases
                a=complexarr(NewDim(0),NewDim(1),NewDim(2)) ; Auxiliary variable with the final matrix size and number of slices
                b=fft(reform(data(*,*,*,i)),-1) ; FFT of the aliased data
                b=shift(b,-round(outputRes(0)/2.),-round(outputRes(1)/2.),-round(outputRes(3)/2.)) ; FFT shifting
                a(SlicePos(0):SlicePos(0)+outputRes(0)-1,SlicePos(1):SlicePos(1)+outputRes(1)-1,$
                  SlicePos(2):SlicePos(2)+outputRes(3)-1)=b ; Assign array a to b
                ;a=a*Filter
                a=shift(a, round(NewDim(0)/2.),round(NewDim(1)/2.),round(NewDim(2)/2.)) ; Shift a
                a=fft(a,1) ; IFFT
                Inter(*,*,*,i)=abs(a) ; Assign the absolute value of a to the auxiliary variable with the final matrix size, number of slices and cardiac phases
                a=0
                b=0
            endfor

            data=fltarr(info.npixels(0), info.npixels(2), info.npixels(1),outputRes(4)) ; Float array with dimensions the matrix size, the number of slices and the cardiac phases

            if NewDim(0) lt info.npixels(0) then begin ; See if the dimensions along Y is lower than the number of pixels along that direction
             SlicePos=[round((info.npixels(0)-NewDim(0))/2.),round((NewDim(1)-info.npixels(2))/2.),round((info.npixels(1)-NewDim(2))/2.)] ; Slice position
               data(SlicePos(0):SlicePos(0)+NewDim(0)-1,*,SlicePos(2):SlicePos(2)+NewDim(2)-1,*)=$
                    Inter(*,SlicePos(1):SlicePos(1)+info.npixels(2)-1,*,*) ; Original data in Y and X dimensions is assigned to Inter
            endif else begin
              SlicePos=[round((NewDim(0)-info.npixels(0))/2.),round((NewDim(1)-info.npixels(2))/2.),round((info.npixels(1)-NewDim(2))/2.)] ; Slice position
              data(*,*,SlicePos(2):SlicePos(2)+NewDim(2)-1,*)=$
                    Inter(SlicePos(0):SlicePos(0)+info.npixels(0)-1,SlicePos(1):SlicePos(1)+info.npixels(2)-1,*,*) ; Original data in Y dimension is assigned to Inter
            endelse
            Inter=0
            data=abs(data) ; Obtain the absolute value of data
            maxDynImage=max(data(*,*,NewDim(2)/2-3:NewDim(2)/2+3,*)) ; Maximum value of the final reconstructed image

            data=data<(maxDynImage)

            data=((data*4000./(maxDynImage))<4000.)>0 ; Data normalization over its maximum value

            openw,1, dir+'time.txt' ; Open a file to write and save the times, file identifier=1
            printf,1,  "ReadingTime  SensesitivitiesTime   ReconTime"
            printf,1,  ReadingTime, SensesitivitiesTime, ReconTime ; Print in the file the reading time, the time spent to compute the sensitivity maps and the time spent reconstructing the image
            close,1
            help, data
;            write_dicom_volume, data, dir
            View3D, dir=dir, data=data ; 3D visualization of the final reconstructed images

       endif else begin
            data=abs(data) ; Take the absolute value of the aliased data
            maxDynImage=max(data(*,*,outputRes(3)/2-3:outputRes(3)/2+3,*)) ; Maximum value of the final reconstructed image

            data=data<(maxDynImage)

            data=((data*4000./maxDynImage)<4000.)>0 ; Data normalization over its maximum value
            ; Data is reversed along rows and columns
            data=reverse(data,1)
            data=reverse(data,2)
            View3D, dir=dir,data=data ; 3D visualization of the final reconstructed images

       endelse
  endif
end