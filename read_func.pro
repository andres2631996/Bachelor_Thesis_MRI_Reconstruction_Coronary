
Function ReadString, InString,i, TR=TR, NoNum=NoNum
; Receives as input a string, an index, TR and not to look for numbers
  if strmatch(InString, '*;*') ne 1 then begin ; Looks for lines in the file without ";"
    ArrayStr=strsplit(InString, /extract, escape=',') ; The "," is removed from the splitted strings
    readf,1,InString ; Read the final string
    i=i+1
  endif else ArrayStr=strsplit(InString, /extract, escape=';') ; If the string contains a ";", it separates the string from the following character
  if keyword_set(NoNum) then begin ; Looks if variable "NoNum" is defined, if there are no numbers
      return, (strsplit(ArrayStr(2), /extract, escape='"')) ; If no numbers are found, the second character in the string is splitted and "''''" are removed from the splitted strings
  endif else begin ; If the string contains numbers
    ;if keyword_set(TR) gt 0 then return, fltarr(strsplit(ArrayStr(4), /extract, escape='"')) else $
    ; If TR exists and is positive, then split the string into numbers and remove "''''", converting it into float. This is the output of the function
    si=size(ArrayStr,/dim)
    return, float(ArrayStr(si(0)-1))
    ; If TR is not specified as input, the whole string is taken and is converted into float. This is the output of the function
  endelse
end

Function ReadProtocol, name
; Receives as input the directory and the name of the protocol
  close,/all ; Closes all
  nlines = (file_lines(name))(0) ; Number of lines of the protocol
  rstring='' ; String used to read the protocol
  Protocol = { $ ; Struct that contains all the data from the protocol
              FOV:              fltarr(3),$; FH, RL, AP
              VoxelSize:        fltarr(3),$; FH, RL, AP
              OffCenter:        fltarr(3),$; FH, RL, AP in the protocol is the opposite AP, RL FH
              SENSE:            fltarr(2)+1.,$
              Angulation:       fltarr(3),$; FH, RL, AP in the protocol is the opposite AP, RL FH
              ReconVoxel:       0.0,      $
              CardPhases:       0,        $
              HalfScan:         fltarr(2)+1.,$
              CardFrequency:    0.,        $
              TFEShotDur:       0.,        $
              TriggerDevice:''   ,        $
              TR: 0.0     $
              }
  openr, 1,name ; Open the protocol

  for i=0, nlines-1 do begin
    readf,1,rstring ; Read the protocol
    if strmatch(rstring, '*FOV*') then begin ; Fill the info of the FOV
      Protocol.FOV(0)=ReadString(rstring,i) ; X component of FOV
      readf,1,rstring
      Protocol.FOV(1)=ReadString(rstring,i) ; Y component of FOV
      readf,1,rstring
      Protocol.FOV(2)=ReadString(rstring,i) ; Z component of FOV
      i += 2
    endif
    if strmatch(rstring, '*Voxel size*') then begin ; Fill the info of the voxel size
      Protocol.VoxelSize(0)=ReadString(rstring,i) ; X component of voxel size
      readf,1,rstring
      Protocol.VoxelSize(1)=ReadString(rstring,i) ; Y component of voxel size
      readf,1,rstring
      Protocol.VoxelSize(2)=ReadString(rstring,i) ; Z component of voxel size
      i += 2
    endif

    if strmatch(rstring, '*SENSE =*') then begin ; Fill the info of SENSE
      if strmatch(rstring, '*yes*') then begin
         readf,1,rstring
         Protocol.SENSE(0)=ReadString(rstring,i) ; P reduction
         readf,1,rstring
         readf,1,rstring
         Protocol.SENSE(1)=ReadString(rstring,i) ; S reduction
      endif
       i += 3
    endif

    if strmatch(rstring, '*TFE factor*') then Protocol.TFEShotDur=ReadString(rstring,i) ; Fill the info of the TFE factor
    if strmatch(rstring, '*Recon voxel size*') then Protocol.ReconVoxel=ReadString(rstring,i) ; Fill the info of the reconstructed voxel size
    if strmatch(rstring, '*heart phases*') then Protocol.CardPhases=ReadString(rstring,i) ; Fill the info of the number of heart phases
    if strmatch(rstring, '*Cardiac frequency*') then Protocol.CardFrequency=ReadString(rstring,i) ; Fill the info of the cardiac frequency
    ;if strmatch(rstring, '*Act. TR/TE (ms)*') then Protocol.TR=ReadString(rstring,i, /TR) ; Fill the info of the TR/TE
    if strmatch(rstring, '*device =*') then Protocol.TriggerDevice=ReadString(rstring,i, /NoNum) ; Fill the info of the triggering of the device
    if strmatch(rstring, '*Halfscan*') then begin
      a=strsplit(rstring, /extract, escape=';')   ; Splits the string in two whenever a ";" appears
      if a(n_elements(a)-1) eq '"yes"' then begin ; Determines if the last string splitted in "a" is "yes"
        readf,1,rstring
        Protocol.HalfScan(0)=ReadString(rstring,i) ; Fill the info of the factor Y of the halfscan
        readf,1,rstring
        Protocol.HalfScan(1)=ReadString(rstring,i) ; Fill the info of the factor Z of the halfscan
        i += 2
      endif
    endif

    if strmatch(rstring, '*TFE dur. shot*') then begin
        TFEShotDur=ReadString(rstring,i) ; Fill the info of the shot duration
        if TFEShotDur gt Protocol.TFEShotDur then Protocol.TFEShotDur=TFEShotDur
    endif
    if strmatch(rstring, '*TFE shot interval*') then begin
        TFEShotDur=ReadString(rstring,i) ; Fill the info of the shot interval
        if TFEShotDur gt Protocol.TFEShotDur then Protocol.TFEShotDur=TFEShotDur
    endif

    if strmatch(rstring, '*Slice Offc*') or strmatch(rstring, '*Stack Offc*') then begin ; Obtain the stack information
      Protocol.OffCenter(2)=ReadString(rstring,i) ; Z component
      readf,1,rstring
      Protocol.OffCenter(1)=ReadString(rstring,i) ; Y component
      readf,1,rstring
      Protocol.OffCenter(0)=ReadString(rstring,i) ; X component
;      Read also the angulation
      readf,1,rstring
      Protocol.Angulation(2)=ReadString(rstring,i) ; Z component
      readf,1,rstring
      Protocol.Angulation(1)=ReadString(rstring,i) ; Y component
      readf,1,rstring
      Protocol.Angulation(0)=ReadString(rstring,i) ; X component
      i += 5
    endif

  endfor
  close, /all
  return, Protocol ; The whole struct is the output of the function
end

Function read_cpx,CoilHeaderCPX, CoilDataCPX, QBodyCoil, NCoils
    print,'Reading SENSE ref scan'
    close,/all
    nlines = (file_lines(CoilHeaderCPX))(0) ; Number of lines of the cpx
    openr, 1, CoilHeaderCPX ; With index 1, open the coil header
    openr, 2, CoilDataCPX ; With index 2, open the coil data
    rstring=''
    for i=0., nlines-1 do begin
      readf,1,rstring
      if strmatch(rstring, '*#*') ne 1 then begin ; Look for lines not beginning by #
        if strmatch(rstring, '*STD*') then begin ; STD vector
          header=long((STRSPLIT(rstring, /extract))(1:*)) ; Vector keeping all the vector attributes associated to STD
          if header(4) eq 1 then begin ; See if the location number attribute is 1
            point_lun, 2, header(19) ; Drive the pointer until the end of the line
            readu, 2, vector ; Vector containing the coil data
            QBodyCoil(*,header(8),header(9))+=(vector) ; QBodyCoil is a vector with the coil data of the Q-body coil
            if header(5) gt 16 and NCoils gt 16 then coil_shift=0 else coil_shift=16 ; Determines whether the synchronization number attribute and number of coils are over 16
            ; If both conditions happen, the coils are not shifted, otherwise they are shifted by a factor of 16
          endif else begin
            point_lun, 2, header(19) ; Drive the pointer to the end of the line
            readu, 2, vector ; Vector containing the coil data
            coils(*,header(8),header(9), header(5)-coil_shift)+=(vector) ; Contains information for coils which are not Q-body, with a location number attribute different than 1
          endelse
        endif else begin
          if strmatch(rstring, '*X-resolution*') then xdim=fix((strsplit(rstring, /extract))(6)) ; Resolution in X
          if strmatch(rstring, '*Y-resolution*') then ydim=fix((strsplit(rstring, /extract))(6)) ; Resolution in Y
          if strmatch(rstring, '*Z-resolution*') then begin ; Resolution in Z
            zdim=fix((strsplit(rstring, /extract))(6))
            coils=complexarr(xdim,ydim,zdim,NCoils) ; Complex array with size the Resolutions in X, Y, and Z and the number of coils. It is an array that contains information for each coil
            QBodyCoil=complexarr(xdim,ydim,zdim) ; Complex array with size the Resolutions in X, Y and Z. It contains the information of the Q-body coil
            vector=complexarr(xdim) ; Complex array with size the Resolution in X
          endif
          if strmatch(rstring, '*0  number of coil channels*') then NCoils=fix((strsplit(rstring, /extract))(9)) ; Number of coil channels
        endelse
      endif
    endfor
    close,/all

    CopyCoil=complexarr(ydim,zdim,xdim, NCoils)
    CopyQBodyCoil=complexarr(ydim,zdim,xdim)
    for i=0, xdim-1 do begin
        CopyCoil(*,*,i,*)=reform(coils(xdim-1-i,*,*,*),ydim,zdim,1,NCoils) ; Change the dimensions of the vector with the information of the coils to YxZ, being X the direction with the slices
        CopyQBodyCoil(*,*,i)=reform(QBodyCoil(xdim-1-i,*,*)) ; Do the same with the Q-body coil
    endfor

    coils=CopyCoil
    QBodyCoil=CopyQBodyCoil

    CopyQBodyCoil=0
    CopyCoil=0

    return, coils ; The function returns the vector with the information of the coils in YxZ, being X the direction with the slices
end
Function read_raw_list_by_TFEFactor,DataHeader, DataData, outputRes=outputRes, slice=slice,ScanProto=ScanProto,  noise, SENSEFactor=SENSEFactor, mask=mask, TimeMask=TimeMask;TFEFactor=TFEFactor,
; Receives as input the header of the cpx, the data, outputRes, the slices, the information of the Scan protocol, the SENSE factor, the mask and the time mask
    print,'Reading Raw Data'
    close,/all
    nlines = (file_lines(DataHeader))(0) ; Number of lines of the file
    openr, 1, DataHeader ; Open the header, index 1
    header=strarr(21) ; String array of 21 characters
    rstring=''
;    RR=0
;    RR_interval=(60000.0/float(ScanProto.CardFrequency))/float(ScanProto.CardPhases); Time between Rs in the ECG. Inverse of the cardiac frequency
;    RR_interval=[RR_interval*0.6,RR_interval*1.4]; Vector with the minimum and maximum RR intervals

    MinCoil=0 ; Minimum coil
    for i=0, 300 do begin
      readf,1,rstring
      if strmatch(rstring, '*NOI*') and (strmatch(rstring, '*#*') ne 1) then begin ; Look for lines beginning by NOI but not beginning by #
        header=long((STRSPLIT(rstring, /extract))(1:*)) ; Split the lines into separate strings when finding a 1
       ; RR=[RR,header(17)]
        MinCoil=[MinCoil,header(5)] ; Information for the minimum coil is in the 5th character of each line
      endif
    endfor
   ; RR_interval=(mean(RR(1:*)))/outputRes(4)
    MinCoil=min(MinCoil(1:*))  ; Minimum value of the vector with the coils with the minimum values
    ;                   Ky          Kz            nCoil       Slices      Phases

    close,/all
    openr, 1, DataHeader
    openr, 2, DataData ; Open the data, index 2
    VectorPosition=0
    KzPos=0
    KyPos=0
    CoilChan=0
    RRInt=0
    RRTop=0

    for i=0., nlines-1 do begin
        readf,1,rstring
        if strmatch(rstring, '*#*') ne 1 then begin ; Look for lines not starting by #
          if strmatch(rstring, '*STD*') then begin;or strmatch(rstring, '*REJ*') ; Look for lines starting by STD
                header=long((STRSPLIT(rstring, /extract))(1:*)) ; Read line by line

                header(8)=header(8)*SENSEFactor(0) ; k-space coordinate in 1st encoding direction (Y), affected by the SENSE factor in Y
                header(9)=header(9)*SENSEFactor(1) ; k-space coordinate in 2nd encoding direction (Z), affected by the SENSE factor in Z

                KyPos=[KyPos,header(8)] ; Vector with the Y positions in the k-space
                KzPos=[KzPos,header(9)] ; Vector with the Z positions in the k-space
;                if header(8) ge 0 then KyPos=[KyPos,header(8)] else KyPos=[KyPos,outputRes(0)+header(8)]
;                if header(9) ge 0 then KzPos=[KzPos,header(9)] else KzPos=[KzPos,outputRes(1)+header(9)]

;                RRInt=[RRInt,header(16)] ; Vector with the RR interval lengths
;                RRTop=[RRTop,header(17)] ; Vector with the values of the R waves of the ECG
                CoilChan=[CoilChan,header(5)] ; Vector with the synchronization channel numbers
                VectorPosition=[VectorPosition,header(19)] ; Vector position determined by the data vector offset (in bytes)
          endif else begin
            if strmatch(rstring, '*NOI*') then begin
                header=long((STRSPLIT(rstring, /extract))(1:*))
                dummy=complexarr(round(header(18)/8.)) ; Complex array of length in bits
                point_lun, 2, header(19)
                readu, 2, dummy ; Dummy receives the information from data
                if keyword_set(noise) ne 1 then begin ; Determines if "noise" has been defined
                  noise=dcomplexarr(ncoils, round(header(18)/8)) ; Double-precision, complex array of size the number of coils x length of "dummy" vector
                  noise(header(5)-MinCoil, *)=dummy ; Noise matrix is filled with the noises from all the coils used
                endif else noise(header(5)-MinCoil, *)=dummy ; If "noise" already exists, it receives the same value than it receives in the case it does not exist
            endif
            if strmatch(rstring, '*kx_range*') then begin
                xdim=fix((strsplit(rstring, /extract))(7))-fix((strsplit(rstring, /extract))(6))+1 ; X dimensions in the k-space
                DataVector=complexarr(xdim) ; Complex vector with the length of the X dimensions in the k-space
            endif
            if strmatch(rstring, '*ky_range*') then ydim=fix((strsplit(rstring, /extract))(7))*SENSEFactor(0);-fix((strsplit(rstring, /extract))(6))+1 ; Y dimensions, affected by the first component of the SENSE factor
            if strmatch(rstring, '*kz_range*') then zdim=fix((strsplit(rstring, /extract))(7))*SENSEFactor(1);-fix((strsplit(rstring, /extract))(6))+1 ; Z dimensions, affected by the second component of the SENSE factor
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

    VectorPosition=VectorPosition(1:*) ; Saves the number of line where the file is being read

    outputRes(0)= ceil((ydim*2.)/(SenseFactor(0)*2.0))*(SenseFactor(0)*2.0) ; Y component of the matrix size (double of the dimensions previously calculated)
    outputRes(1)= ceil((zdim*2.)/(SenseFactor(1)*2.0))*(SenseFactor(1)*2.0) ; Z component of the matrix size (double of the dimensions previously calculated)

    KyPos=KyPos(1:*) ; Vector with the Y dimensions in the k-space
    ind=where(KyPos lt 0, count) ; Indexes where the Y positions in the k-space are 0, they are also counted
    if count gt 0 then KyPos(ind)=outputRes(0)+KyPos(ind) ; If there are positions which are 0 in the Y dimension in the k-space, they are assigned the value of the number of rows of the final image
    ; Same in Z dimension of the k-space
    KzPos=KzPos(1:*)
    ind=where(KzPos lt 0, count)
    if count gt 0 then KzPos(ind)=outputRes(1)+KzPos(ind)

    CoilChan=CoilChan(1:*) ; Vector with the synchronization channel numbers
    ;RRInt=RRInt(1:*) ; Vector with the RR intervals
    ;RRTop=RRTop(1:*) ; Vector with the values of the R waves of the ECG

    CoilChan=CoilChan-min(CoilChan) ; The first synchronization channel number is set to 0 in case it is not 0

    nlines=n_elements(VectorPosition) ; Number of lines of the file
    n_channels=max(coilchan)-min(coilchan)+1 ; The number of channels used is the difference between the maximum and minimum channel synchronization numbers

    Data=dcomplexarr(outputRes(0),outputRes(1),outputRes(2),outputRes(3)/10,1) ; Final output with the matrix size, the number of cardiac phases, the number of slices and the number of coils used. It is saved in a complex array
    ;TimeMask=fltarr(outputRes(0),outputRes(1),outputRes(2),1,outputRes(4)) ; Floating array of the same dimensions than "data" but only in one slice
    mask=fltarr(outputRes(0),outputRes(1),outputRes(2),1,outputRes(4)) ; Floating array of the same dimensions than "data" but only in one slice
    filter=exp(-((findgen(xdim)-xdim/2)/xdim)^2) ; Gaussian filter along dimensions of X. It is a 1D vector with length the number of slices

    openr, 2, DataData ; Open the data file, index=2
    ppppp=0
    filter=fltarr(xdim)
    filter=exp(-(findgen(xdim)-xdim/2.)^2./20000)
    filter(0:4)=0.0 ; First 5 elements of the filter, which will be used in the edges, will be 0
    filter(xdim-6:*)=0.0 ; Last 5 elements of the filter, which will be used in the edges, will be 0, too

    filter=convol(filter, FilterGen1D(30,30),/EDGE_t) ; Final filter is convolution of the original 1D filter with a Gaussian filter of 30 samples and width 30
    ;RRInterval=mean(rrtop)/outputRes(4) ; The RR interval is defined as the mean value of the distances between R waves in the ECG and the number of cardiac phases

    for i=0., nlines-1 do begin
       point_lun, 2, VectorPosition(i) ; Moves the pointer to the beginning of the file
       readu, 2, DataVector ; Reads each line the data file

       DataVector=DataVector*filter ; Filtered positions in the X dimension
       DataVector=shift(fft(shift(DataVector, -xdim/2),/inverse, /OVERWRITE ), xdim/2)
       ; Shift the filtered positions in the X dimension, compute the IFFT and undo the shifting
       DataPoint=DataVector(xdim/2-floor(outputRes(3)/2.)+slice:xdim/2-floor(outputRes(3)/2.)+outputRes(3)-1+slice)
       ; Data in a point is defined by the components of the vector with the data, "datavector" in half of the dimensions of X and half of the number of coils plus the slice number
       ; And that same value and the number of coils
;       RRInterval=RRTop(i)/outputRes(4)

       ;if (RRInterval gt RR_interval(0)) and (RRInterval lt RR_interval(1)) then begin ; Look if the current RR interval is longer than the initial but shorter than the following
           ; PhaseN=floor(RRInt(i)/RRinterval) ; The cardiac phase number is determined as the current RR interval over the whole RR interval
           ;if PhaseN eq outputRes(4) then PhaseN=outputRes(4)-1 ; If the previous cardiac phase number equals the cardiac phase number in OutputRes, the phase number is one lower than the one in the vector
           ;if PhaseN eq outputRes(4)+1 then begin ; If the previous cardiac phase number is one above the cardiac phase number in OutputRes, it is set to 0, as a new cycle starts
              ;PhaseN=0
              ;RRInt(i)= RRInt(i) -outputRes(4)*RRInterval ; In order to avoid acummulated values from other cycles, the following RR intervals are substracted the cardiac phase number in OutputRes times the whole RR interval
           ;endif
           ;if PhaseN lt outputRes(4) then begin ; Looks if the previous cardiac phase number is lower than the cardiac phase number in OutputRes
             ;if  (mask(KyPos(i), KzPos(i), CoilChan(i), 0, 0) lt 1) then begin ; If "mask" in the first slice is lower than 1
                Data(KyPos(i), KzPos(i), CoilChan(i),*,0) += DataPoint ; Data in each point are actualized to the previous data calculated for each point
                ;TimeMask(KyPos(i), KzPos(i), CoilChan(i), 0, PhaseN) += RRInt(i) ; Time is actualized to the new RR interval
                mask(KyPos(i), KzPos(i), CoilChan(i), 0, 0) +=1.0 ; Mask is set to 1
             ;endif
           ;endif

    endfor
    close,/all
    ind=where(mask gt 0, count) ; Indexes where mask is positive
    print, max(mask)   ; Print the maximum value of the mask
;    if count gt 0 then TimeMask(ind)=TimeMask(ind)/mask(ind)
;    maskData=rebin(mask, outputRes(0),outputRes(1),outputRes(2),outputRes(3),outputRes(4))
;    ind=where(maskData gt 0, count)
;    if count gt 0 then data(ind)=data(ind)/maskData(ind)
;    maskData=0
; All parameters are set to 0 again
    VectorPosition=0
    KzPos=0
    KyPos=0
    CoilChan=0
    RRInt=0
    RRTop=0

    ;TimeMask=reform(TimeMask(*,*,0,0,*)) ; Time mask dimensions are changed and set to the first coil and slice
  return, Data ; Data with the matrix size in the k-space, the number of coils and slices and the cardiac phase number are returned
end

Function MakeStatic, data, outputRes=outputRes,TimeMask=TimeMask,ECG_PPU=ECG_PPU
; Receives as input the data, outputRes, the time mask and the ECG data
  print, 'Make Static Dataset'
  ;     Ky          Kz            nCoil       Slices      Phases
  ;outputRes(0),outputRes(1),outputRes(2),outputRes(3),outputRes(4));
  static=reform(data(*,*,*,*,0)) ; Data is adapted to a new set of dimensions if the cardiac phase is just one

  return, static ; The function returns static data read
end

function read_dicom_info, filename
; Function that reads DICOM images
   if (FILE_INFO(filename)).exists then begin ; See if the file we want to read really exists
        info={$
              CDV: fltarr(3),       $
              pixelsize: fltarr(3),  $ ; Pixel size
              npixels: intarr(3)     $ ; Number of pixels
             }
         oImg = OBJ_NEW('IDLffDicomEx', filename) ; Image to read
         ; Obtain the number of pixels in the 3 directions: Matrix Size X Number of Slices
         info.npixels(0)= oImg->GetValue('0028,0010', COUNT=vCount)
         info.npixels(1)= oImg->GetValue('0028,0011', COUNT=vCount)
         info.npixels(2)= oImg->GetValue('2001,1018', COUNT=vCount)
         ; Obtain the pixel size in the 3 directions
         info.pixelsize(0:1)= oImg->GetValue('0028,0030', COUNT=vCount)
         info.pixelsize(2)= oImg->GetValue('0018,0088', COUNT=vCount)
         ; Obtain the FOV in the 3 directions
         info.CDV(*)=float(info.npixels)*float(info.pixelsize)
         OBJ_DESTROY, oImg
         ; Destroy the image to read and keep just the information
         return, info ; Return the information of each slice: Pixel size, number of pixels and FOV
   endif else return,0
end

pro ModDicomTAG, Tag, value, oImg
    VR=oImg->GetVR(Tag)
    oImg->SetValue,Tag,VR,value
end
pro write_dicom_volume, DynImage, dir
; Procedure receiving as inputs the dynamic image and the directory where to save it
print, "Generating DICOM images"
    if (FILE_INFO(dir+'DICOM\')).exists then begin ; See if the directory where to save the images exists
        filter='IM_*' ; Image naming
        recurrent_dir=''

        if (file_info(dir+'BrisaDicom')).exists then file_delete, dir+'BrisaDicom', /RECURSIVE ; If there are files in BrisaDicom directory, they will all be deleted

        file_mkdir, dir+'BrisaDicom' ; Create a new directory
        if (FILE_INFO(dir+'DICOM\00000001')).exists then begin ; See if there are files in the new directory
            file_mkdir, dir+'BrisaDicom\00000001' ; Create a new directory
             recurrent_dir=[recurrent_dir,'00000001\'] ; Save the image files in here
        endif
        ; Do the same again
        if (FILE_INFO(dir+'DICOM\00000002')).exists then begin
            file_mkdir, dir+'BrisaDicom\00000002'
            recurrent_dir=[recurrent_dir,'00000002\']
        endif

        for k=0, n_elements(recurrent_dir)-1 do begin ; Go through all the image files saved in the directories
          print, dir+'DICOM\'+recurrent_dir(k) ; Print the directory and the current image that is being read
            filename=file_search(dir+'DICOM\'+recurrent_dir(k)+filter,count=count) ; Search all files in the directory matching the name we are looking for
            if count gt 0 then begin ; See if there are files in the directory we are looking
              oImg = OBJ_NEW('IDLffDicomEx',filename(0)) ; Create an object related to each file
                n_slices= fix(oImg->GetValue('2001,1018', COUNT=vCount))-1 ; Obtain the number of slices
                slicethi=oImg->GetValue('0018,0088', COUNT=vCount) ; Obtain the slice thickness
              OBJ_DESTROY, oImg ; Eliminate the object and keep the information

              for i=0, count-1 do begin ; Go through all the files in the directory
                   name=strmid(filename(i),STRPOS(filename(i), 'IM_'),7) ; Name of the file is inside filename, for that reason we have to find where the name starts and then separate it from the rest of the string
                   print, name ; Print the name
                   oImg = OBJ_NEW('IDLffDicomEx', dir+'BrisaDicom\'+recurrent_dir(k)+name, CLONE=filename(i)) ; Create an object of the image that is being read
                     slice= oImg->GetValue('2001,100A', COUNT=vCount)
                     phase= oImg->GetValue('2001,1008', COUNT=vCount)
;                     print, slice,',', phase
                     ModDicomTAG, '0028,1050', '5000', oImg;WindowCenter
                     ModDicomTAG, '0028,1051', '10000', oImg;WindowWidth
                     ModDicomTAG, '0018,0050', slicethi, oImg;slice thickness
                     ModDicomTAG, '0018,1030', '3D-ESSOS', oImg;StoredBit
;                     ModDicomTAG, '0028,0102', '15', oImg;High Bit
                     ModDicomTAG, '0028,1052', '0.0', oImg;RescaleIntercept
                     ModDicomTAG, '0028,1053', '1.0', oImg;RescaleSlope
                     ModDicomTAG, '2005,100E', '1.0', oImg;ScaleSlope
                     imagen=uint(reform(DynImage(*,slice-1,*,phase-1))) ; Final image is unsigned integer datatype
                     imagen=reverse(imagen,2) ; Reverse the columns of the image
                     oImg->SetPixelData, imagen
                     oImg ->Commit
                   OBJ_DESTROY, oImg ; Remove the object and keep the information
;                  imagencita=read_dicom(dir+'BrisaDicom\'+name)
;                  tvscl, imagencita
              endfor
           endif
        endfor
    endif else print, 'Dicom files not saved' ; In case no files are found in the directory
end